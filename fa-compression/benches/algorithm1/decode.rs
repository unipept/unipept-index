use std::hint::black_box;

use fa_compression::algorithm1::{decode, decoded, encode};

use super::util::generate_annotation;

/// Generate `count` random annotations, encoded.
fn generate_encoded_annotations(count: usize) -> Vec<u8> {
    let mut random = rand::rng();

    let mut annotations = String::new();
    for _ in 0..count {
        annotations.push_str(&generate_annotation(&mut random));
        annotations.push(';');
    }
    annotations.pop();

    encode(annotations.as_str())
}

/// The size a real protein's annotation set actually is.
///
/// Measured over a full Unipept protein dump (the `proteins.tsv` the index is built from, which is
/// not in this repo): the decoded text averages ~143 characters, with a median of 117 and a p90 of
/// 221 — roughly 10 annotations, not the 100 that `decode_algorithm1` below uses. The big case is
/// worth keeping (it is where per-call overhead disappears and the inner loop is all that is left),
/// but the server decodes one of these per protein hit and the small one is what it pays.
const TYPICAL_ANNOTATION_COUNT: usize = 10;

pub fn decode_benchmark(c: &mut criterion::Criterion) {
    c.bench_function("decode_algorithm1", |b| {
        b.iter_batched(
            || generate_encoded_annotations(100),
            |annotations| black_box(decode(annotations.as_slice())),
            criterion::BatchSize::SmallInput
        )
    });

    c.bench_function("decode_algorithm1_typical", |b| {
        b.iter_batched(
            || generate_encoded_annotations(TYPICAL_ANNOTATION_COUNT),
            |annotations| black_box(decode(annotations.as_slice())),
            criterion::BatchSize::SmallInput
        )
    });

    // The path `sa-server` takes: decoded straight into a writer, never materialised as a `String`.
    c.bench_function("decoded_algorithm1_typical_streamed", |b| {
        b.iter_batched(
            || generate_encoded_annotations(TYPICAL_ANNOTATION_COUNT),
            |annotations| {
                use std::fmt::Write;
                let mut out = String::with_capacity(512);
                write!(out, "{}", decoded(annotations.as_slice())).unwrap();
                black_box(out)
            },
            criterion::BatchSize::SmallInput
        )
    });
}
