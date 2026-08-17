//! This module provides a function to decode a byte array into a string representation of
//! annotations.
//!
//! # Shape of the decoder
//!
//! The encoded form is a nibble stream: three `,`-separated sections (EC, GO, InterPro), each a
//! `;`-separated list of annotations with their prefix stripped. Decoding therefore has to put the
//! prefixes back, and the obvious way to do that — materialise the nibble stream as a `String`,
//! `split(',')`, `split(';')`, and reassemble — costs two allocations and two passes over the data.
//!
//! This does it in one pass with no intermediate buffer, which matters because the server decodes
//! one of these per protein hit and a single request can have millions of them.
//!
//! The trick that makes a single pass possible: the old code appended `prefix + annotation + ';'`
//! for every annotation and then popped the trailing `;`. That final pop is unimplementable in a
//! streaming writer, where the byte may already have been flushed. But the output it produces is
//! exactly the `;`-join of `prefix ++ annotation`, so emitting the separator *before* each
//! annotation instead of after removes the pop entirely.

use std::fmt;

/// The prefixes for the different types of annotations.
static PREFIXES: [&[u8]; 3] = [b"EC:", b"GO:", b"IPR:IPR"];

/// Nibble value to the character it decodes to, in [`CharacterSet`] order.
///
/// The same mapping as [`CharacterSet::decode`], as a table rather than a `match` with a panicking
/// arm — every one of the 16 values is valid, so there is nothing to panic on. `decode_matches_the_character_set`
/// pins the two together.
const NIBBLE: [u8; 16] = *b"$0123456789-.n,;";

/// Byte to both of its characters, so the common case writes two bytes in one call.
const PAIR: [[u8; 2]; 256] = build_pairs();

/// Byte to whether either of its nibbles needs the slow path: `$` (0, padding), `,` (14, a section
/// break) or `;` (15, an annotation break).
///
/// Real annotation payloads are digits, `.`, `-` and `n`, so this is false for the overwhelming
/// majority of bytes and the hot loop is a table lookup and a two-byte write.
const SPECIAL: [bool; 256] = build_special();

const fn build_pairs() -> [[u8; 2]; 256] {
    let mut table = [[0u8; 2]; 256];
    let mut byte = 0;
    while byte < 256 {
        table[byte] = [NIBBLE[byte >> 4], NIBBLE[byte & 0b1111]];
        byte += 1;
    }
    table
}

const fn build_special() -> [bool; 256] {
    let mut table = [false; 256];
    let mut byte = 0;
    while byte < 256 {
        let (high, low) = (byte >> 4, byte & 0b1111);
        table[byte] = high == 0 || high == 14 || high == 15 || low == 0 || low == 14 || low == 15;
        byte += 1;
    }
    table
}

/// Where [`decode_to`] writes. Every slice it is handed is ASCII.
///
/// Deliberately private: it exists so the allocating and the streaming entry points can share one
/// decoder, not as something callers implement.
trait Sink {
    /// Appends `bytes`, which are always ASCII.
    fn write(&mut self, bytes: &[u8]);
    /// Whether writing has failed and the decoder should stop early. Only the formatter sink can
    /// fail; the `String` one never does.
    fn failed(&self) -> bool {
        false
    }
}

impl Sink for String {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        debug_assert!(bytes.is_ascii(), "the decoder only ever emits ASCII");
        // Every byte reaching here comes from `NIBBLE`/`PAIR` (a 16-entry ASCII table) or from
        // `PREFIXES` (ASCII literals), so this is always valid UTF-8. `from_utf8` on an all-ASCII
        // slice is a cheap vectorised scan, which is why this stays safe code.
        self.push_str(std::str::from_utf8(bytes).expect("the decoder only ever emits ASCII"));
    }
}

/// Buffers into a fixed array and flushes through [`fmt::Formatter::write_str`], so a typical
/// annotation set reaches the underlying writer in one call and nothing is heap-allocated.
struct FmtSink<'a, 'b> {
    formatter: &'a mut fmt::Formatter<'b>,
    buffer: [u8; 512],
    len: usize,
    error: Option<fmt::Error>
}

impl<'a, 'b> FmtSink<'a, 'b> {
    fn new(formatter: &'a mut fmt::Formatter<'b>) -> Self {
        Self {
            formatter,
            buffer: [0; 512],
            len: 0,
            error: None
        }
    }

    fn flush(&mut self) {
        if self.len == 0 || self.error.is_some() {
            return;
        }
        // Same ASCII argument as the `String` impl.
        let text = std::str::from_utf8(&self.buffer[..self.len]).expect("the decoder only ever emits ASCII");
        if let Err(err) = self.formatter.write_str(text) {
            self.error = Some(err);
        }
        self.len = 0;
    }

    fn finish(mut self) -> fmt::Result {
        self.flush();
        match self.error {
            Some(err) => Err(err),
            None => Ok(())
        }
    }
}

impl Sink for FmtSink<'_, '_> {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        if self.len + bytes.len() > self.buffer.len() {
            self.flush();
        }
        // A single write never exceeds the buffer: the longest is `PREFIXES[2]` at 7 bytes.
        self.buffer[self.len..self.len + bytes.len()].copy_from_slice(bytes);
        self.len += bytes.len();
    }

    fn failed(&self) -> bool {
        self.error.is_some()
    }
}

/// The decoder itself: one pass over the nibble stream, writing straight into `sink`.
///
/// Reproduces the two-pass `split`-based original exactly, including the cases well-formed
/// [`encode`](super::encode) output never produces:
///
/// * only the LOW nibble of a byte is skipped when it is `$`; a `$` in the high nibble is content,
///   so `decode(&[0])` is `"EC:$"`.
/// * `split(';')` keeps empty elements, so a doubled `;` emits a bare prefix.
/// * an empty section is dropped but still consumes its prefix slot, because the original filtered
///   *after* zipping with `PREFIXES` — `",b"` is `"GO:b"`, not `"EC:b"`.
/// * `zip(PREFIXES)` truncates, so a fourth and later section is silently dropped.
fn decode_to<S: Sink>(input: &[u8], sink: &mut S) {
    // `separator_owed`: something has been written, so the next annotation needs a `;` in front.
    // `open`: the current annotation's prefix has already been written.
    let mut separator_owed = false;
    let mut open = false;
    let mut section = 0usize;

    for &byte in input {
        if sink.failed() {
            return;
        }

        // Fast path: both nibbles are ordinary characters inside an already-open annotation, so
        // they go out as a pair with no decisions to make.
        if open && !SPECIAL[byte as usize] {
            sink.write(&PAIR[byte as usize]);
            continue;
        }

        let [high, low] = PAIR[byte as usize];
        // The high nibble is always content, even when it is `$`; the low one is skipped when it is.
        let characters: &[u8] = if low == b'$' { &[high] } else { &[high, low] };

        for &character in characters {
            match character {
                b',' => {
                    section += 1;
                    open = false;
                    // `zip(PREFIXES)` runs out here — nothing further can be emitted.
                    if section >= PREFIXES.len() {
                        return;
                    }
                }
                b';' => {
                    // A `;` as a section's first character means an empty leading element, which
                    // still gets its prefix.
                    if !open {
                        open_annotation(sink, section, &mut separator_owed);
                    }
                    open_annotation(sink, section, &mut separator_owed);
                    open = true;
                }
                _ => {
                    if !open {
                        open_annotation(sink, section, &mut separator_owed);
                        open = true;
                    }
                    sink.write(&[character]);
                }
            }
        }
    }
}

#[inline]
fn open_annotation<S: Sink>(sink: &mut S, section: usize, separator_owed: &mut bool) {
    if *separator_owed {
        sink.write(b";");
    }
    sink.write(PREFIXES[section]);
    *separator_owed = true;
}

/// Decodes a byte array into a string representation of annotations.
///
/// The input byte array is decoded by splitting each byte into two characters.
/// The decoded annotations are then reconstructed by appending the appropriate prefix
/// (e.g., "EC:", "GO:", "IPR:IPR") to each annotation.
///
/// # Arguments
///
/// * `input` - The byte array to decode.
///
/// # Returns
///
/// A string representation of the decoded annotations.
///
/// # Examples
///
/// ```
/// use fa_compression::algorithm1::decode;
///
/// let input = &[ 44, 44, 44, 190, 17, 26, 56, 174, 18, 116, 117 ];
/// let result = decode(input);
/// assert_eq!(result, "EC:1.1.1.-;GO:0009279;IPR:IPR016364");
/// ```
pub fn decode(input: &[u8]) -> String {
    // Each byte becomes at most two characters, and the prefixes roughly double that again. A hint,
    // not a bound — a `;`-heavy input can still grow the string, exactly as it could before.
    let mut result = String::with_capacity(input.len() * 3);
    decode_into(input, &mut result);
    result
}

/// Decodes a byte array into annotations, **appending** them to `out`.
///
/// The allocation-free variant of [`decode`] for callers that have a buffer to reuse. Nothing is
/// cleared first, and nothing is written for empty input.
///
/// # Arguments
///
/// * `input` - The byte array to decode.
/// * `out` - The string to append the decoded annotations to.
///
/// # Examples
///
/// ```
/// use fa_compression::algorithm1::decode_into;
///
/// let mut out = String::from("annotations: ");
/// decode_into(&[ 225, 17, 163, 138, 224 ], &mut out);
/// assert_eq!(out, "annotations: GO:0009279");
/// ```
pub fn decode_into(input: &[u8], out: &mut String) {
    decode_to(input, out);
}

/// Borrowed encoded annotations that render as their decoded text.
///
/// Writing this through [`fmt::Display`] decodes straight into the formatter with no intermediate
/// `String`, which is what lets a serialiser emit the annotations without materialising them —
/// `serde_json`'s `collect_str` streams a `Display` value directly into its output buffer.
///
/// # Note on formatting flags
///
/// The `Display` impl writes through [`fmt::Formatter::write_str`], which bypasses padding, so
/// width, fill and alignment are ignored: `format!("{:>20}", decoded(x))` is not padded. Formatting
/// it plainly, as `collect_str` and `to_string` do, is the supported use.
///
/// # Examples
///
/// ```
/// use fa_compression::algorithm1::decoded;
///
/// assert_eq!(decoded(&[ 238, 18, 116, 117 ]).to_string(), "IPR:IPR016364");
/// ```
#[derive(Clone, Copy)]
pub struct Decoded<'a>(&'a [u8]);

/// Views encoded annotations as their decoded text, without decoding them yet.
///
/// See [`Decoded`].
///
/// # Arguments
///
/// * `input` - The encoded byte array to view.
///
/// # Examples
///
/// ```
/// use fa_compression::algorithm1::{decode, decoded};
///
/// let input = &[ 44, 44, 44, 190, 17, 26, 56, 174, 18, 116, 117 ];
/// assert_eq!(decoded(input).to_string(), decode(input));
/// ```
pub fn decoded(input: &[u8]) -> Decoded<'_> {
    Decoded(input)
}

impl fmt::Display for Decoded<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut sink = FmtSink::new(formatter);
        decode_to(self.0, &mut sink);
        sink.finish()
    }
}

impl fmt::Debug for Decoded<'_> {
    /// Shows the decoded text rather than the encoded bytes, which is the only readable rendering
    /// of a nibble stream.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{:?}", self.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::{super::CharacterSet, super::Decode, *};

    #[test]
    fn test_decode_empty() {
        assert_eq!(decode(&[]), "")
    }

    #[test]
    fn test_decode_single_ec() {
        assert_eq!(decode(&[44, 44, 44, 190, 224]), "EC:1.1.1.-")
    }

    #[test]
    fn test_decode_single_go() {
        assert_eq!(decode(&[225, 17, 163, 138, 224]), "GO:0009279")
    }

    #[test]
    fn test_decode_single_ipr() {
        assert_eq!(decode(&[238, 18, 116, 117]), "IPR:IPR016364")
    }

    #[test]
    fn test_decode_no_ec() {
        assert_eq!(decode(&[225, 17, 163, 138, 225, 39, 71, 95, 17, 153, 39]), "GO:0009279;IPR:IPR016364;IPR:IPR008816")
    }

    #[test]
    fn test_decode_no_go() {
        assert_eq!(decode(&[44, 44, 44, 191, 44, 60, 44, 142, 225, 39, 71, 80]), "EC:1.1.1.-;EC:1.2.1.7;IPR:IPR016364")
    }

    #[test]
    fn test_decode_no_ipr() {
        assert_eq!(decode(&[44, 44, 44, 190, 17, 26, 56, 175, 17, 26, 56, 174]), "EC:1.1.1.-;GO:0009279;GO:0009279")
    }

    #[test]
    fn test_decode_all() {
        assert_eq!(
            decode(&[44, 44, 44, 190, 17, 26, 56, 174, 18, 116, 117, 241, 67, 116, 111, 17, 153, 39]),
            "EC:1.1.1.-;GO:0009279;IPR:IPR016364;IPR:IPR032635;IPR:IPR008816"
        )
    }

    // ── Byte-identity against the original implementation ─────────────────────
    //
    // The rewrite's whole contract is that it changes nothing about the output, so the old
    // implementation is kept here verbatim and the new one is checked against it. Every branch of
    // the state machine is reachable from arbitrary bytes, including combinations `encode` never
    // produces, so the corpus is exhaustive where it can afford to be.

    /// The pre-rewrite decoder, copied unchanged. Do not "simplify" this — its value is being the
    /// thing the current implementation is not.
    fn decode_reference(input: &[u8]) -> String {
        static REFERENCE_PREFIXES: [&str; 3] = ["EC:", "GO:", "IPR:IPR"];

        if input.is_empty() {
            return String::new();
        }

        let mut decoded = String::with_capacity(input.len() * 2);
        for &byte in input {
            let (c1, c2) = CharacterSet::decode_pair(byte);

            decoded.push(c1);
            if c2 != '$' {
                decoded.push(c2);
            }
        }

        let mut result = String::with_capacity(input.len() * 3);
        for (annotations, prefix) in decoded.split(',').zip(REFERENCE_PREFIXES).filter(|(s, _)| !s.is_empty()) {
            for annotation in annotations.split(';') {
                result.push_str(prefix);
                result.push_str(annotation);
                result.push(';');
            }
        }

        result.pop();

        result
    }

    #[test]
    fn matches_the_reference_on_every_one_byte_input() {
        for byte in 0..=u8::MAX {
            assert_eq!(decode(&[byte]), decode_reference(&[byte]), "input {byte:?}");
        }
    }

    /// Exhaustive over every nibble-pair interaction: `$`, `,` and `;` in either position, section
    /// transitions, padding, and the boundaries between them.
    #[test]
    fn matches_the_reference_on_every_two_byte_input() {
        for first in 0..=u8::MAX {
            for second in 0..=u8::MAX {
                let input = [first, second];
                assert_eq!(decode(&input), decode_reference(&input), "input {input:?}");
            }
        }
    }

    /// Three bytes is 16.7M cases — too slow for every run, but the place to look first if the
    /// two-byte sweep ever passes while something longer misbehaves.
    #[test]
    #[ignore = "16.7M cases; run explicitly with --ignored"]
    fn matches_the_reference_on_every_three_byte_input() {
        for first in 0..=u8::MAX {
            for second in 0..=u8::MAX {
                for third in 0..=u8::MAX {
                    let input = [first, second, third];
                    assert_eq!(decode(&input), decode_reference(&input), "input {input:?}");
                }
            }
        }
    }

    #[test]
    fn matches_the_reference_on_random_inputs() {
        // A fixed seed, so a failure is reproducible rather than a once-seen flake.
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };

        for _ in 0..100_000 {
            let length = (next() % 33) as usize;
            let input: Vec<u8> = (0..length).map(|_| next() as u8).collect();
            assert_eq!(decode(&input), decode_reference(&input), "input {input:?}");
        }
    }

    /// The cases worth naming, so a failure says which rule broke rather than printing bytes.
    #[test]
    fn reproduces_the_original_edge_cases() {
        let encoded = |text: &str| super::super::encode(text);

        // A `$` in the high nibble is content, not padding.
        assert_eq!(decode(&[0]), "EC:$");
        // Empty sections are dropped but still consume their prefix slot.
        assert_eq!(decode(&encoded("GO:0009279")), "GO:0009279");
        // Every section empty.
        assert_eq!(decode(&encoded("nonsense-with-no-prefixes")), "");
        // Round-trip through the real encoder.
        let text = "EC:1.1.1.-;GO:0009279;IPR:IPR016364";
        assert_eq!(decode(&encoded(text)), text);
    }

    // ── The lookup tables ─────────────────────────────────────────────────────

    #[test]
    fn decode_matches_the_character_set() {
        for value in 0..16u8 {
            assert_eq!(NIBBLE[value as usize] as char, CharacterSet::decode(value), "nibble {value}");
        }
    }

    #[test]
    fn pairs_match_the_character_set() {
        for byte in 0..=u8::MAX {
            let (first, second) = CharacterSet::decode_pair(byte);
            assert_eq!(PAIR[byte as usize], [first as u8, second as u8], "byte {byte}");
        }
    }

    #[test]
    fn special_flags_every_byte_the_fast_path_cannot_take() {
        for byte in 0..=u8::MAX {
            let [high, low] = PAIR[byte as usize];
            let expected = matches!(high, b'$' | b',' | b';') || matches!(low, b'$' | b',' | b';');
            assert_eq!(SPECIAL[byte as usize], expected, "byte {byte}");
        }
    }

    // ── The other two entry points ────────────────────────────────────────────

    #[test]
    fn decoded_renders_the_same_text_as_decode() {
        for first in 0..=u8::MAX {
            for second in 0..=u8::MAX {
                let input = [first, second];
                assert_eq!(decoded(&input).to_string(), decode(&input), "input {input:?}");
            }
        }
    }

    /// The `Display` sink buffers in 512-byte chunks, so an input that spans several flushes is the
    /// case that would catch a boundary bug.
    #[test]
    fn decoded_handles_inputs_larger_than_its_flush_buffer() {
        let text = (0..400).map(|i| format!("IPR:IPR{:06}", i)).collect::<Vec<_>>().join(";");
        let encoded = super::super::encode(&text);
        assert!(decode(&encoded).len() > 4 * 512, "the corpus must span several flushes");
        assert_eq!(decoded(&encoded).to_string(), decode(&encoded));
    }

    #[test]
    fn decode_into_appends_rather_than_replacing() {
        let mut out = String::from("annotations: ");
        decode_into(&[225, 17, 163, 138, 224], &mut out);
        assert_eq!(out, "annotations: GO:0009279");

        // Empty input leaves the buffer alone.
        decode_into(&[], &mut out);
        assert_eq!(out, "annotations: GO:0009279");
    }

    #[test]
    fn debug_shows_the_decoded_text() {
        assert_eq!(format!("{:?}", decoded(&[238, 18, 116, 117])), "\"IPR:IPR016364\"");
    }
}
