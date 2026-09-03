# libsais64-rs

![CI](https://img.shields.io/github/actions/workflow/status/unipept/unipept-index/ci.yml?logo=github&label=ci)

Rust bindings to [`unipept/libsais-packed`](https://github.com/unipept/libsais-packed), the
suffix-array construction library `sa-builder` uses. One function does the work: `sais64(text,
sparseness)` returns the suffix array over `text`, or an error.

## Building

The build script (`builder.rs`) fetches the C library, compiles it with CMake, and generates the
FFI bindings with `bindgen`. So building this crate — and therefore `sa-builder`, and therefore a
plain `cargo build` of the workspace — needs:

* `git`, and **network access on the first build** (a shallow fetch of one commit from GitHub)
* `cmake` and `make`
* a C compiler, and `libclang` for `bindgen`

The checkout lands in `libsais64-rs/libsais-packed/`, which is gitignored: it is a build artefact,
not a vendored dependency or a submodule.

**The commit is pinned**, in `LIBSAIS_COMMIT` at the top of `builder.rs`, so that a given Rust
commit always builds the same binary and no change in the C repository reaches a build without a
commit here recording it. Bump it by editing that constant. A build whose checkout is already at
the pinned commit skips the fetch; cmake and make still run, and find their work already done.

## Sparseness and bit packing

A sparse suffix array indexes only every n-th text position. Rather than sampling after
construction, this crate packs n consecutive residues into one symbol and builds a dense suffix
array over the packed text — so libsais does the sampling for free. At 5 bits per residue that
selects one of three entry points by how wide the packed symbol has to be:

| sparseness | packed width | libsais entry point |
|---|---|---|
| 1 | 8 bits | `libsais64` |
| 2-3 | 16 bits | `libsais16x64` |
| 4-6 | 32 bits | `libsais32x64` |

Sparseness 1 puts a single residue in each symbol, so it is not a size reduction, but the text
still goes through the packer: that maps it onto the ranks `0..=27` rather than the scattered ASCII
bytes libsais would otherwise bucket over, and it rejects a byte outside the alphabet instead of
handing it to the algorithm.

Six is the widest this crate can pack, since `6 * 5 = 30` bits still fits a `u32`, and `sais64`
refuses anything outside `1..=6` rather than picking a branch with no packer behind it. `sa-builder`
stops at 5 and reaches any higher `--sparseness-factor` by sampling the result as well.

The returned positions are multiplied back by the sparseness before they are handed out, so callers
see positions in the original text.

`bitpacking` is public and does the packing. It rejects any byte outside the protein alphabet
(`$`, `-`, `A`-`Z`) rather than computing on it: a byte below `A` would underflow the rank
subtraction, and a byte at or above 95 — which includes every non-ASCII UTF-8 byte — would produce
a rank wide enough to spill into the neighbouring residue's field and silently corrupt a character
the caller never supplied. `sa-builder` cannot reach either case, since it packs text that has
already been through the 5-bit protein encoding.

## Where it sits

Depends on `bindgen` at build time and on the C library it clones. Used by `sa-builder` alone.

---

Part of the [Unipept Index](../README.md) workspace · full API docs with
`cargo doc -p libsais64-rs --open`
