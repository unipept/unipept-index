# bitarray

![CI](https://img.shields.io/github/actions/workflow/status/unipept/unipept-index/ci.yml?logo=github&label=ci)

Dense arrays of fixed-width values packed into `u64` words. The index stores hundreds of millions
of values that need far fewer than 64 bits each — a suffix array over a 300 M-residue text needs
29, and the protein text needs 5 for its residue alphabet. Packing them at their natural
width rather than rounding up to a byte is what keeps a compressed index small enough to be worth
loading.

## The layout is the contract

Values are packed **most-significant-bit first within each little-endian `u64` word**, and a value
may straddle a word boundary. Both implementations here and the mmap readers in `protein-text` and
`sa-index` depend on that layout matching exactly.

## Two implementations

| type | width fixed | use when |
|---|---|---|
| `BitArray<BITS>` | at compile time | the width is a property of the data — `protein-text` is always `BitArray<5>` |
| `DynBitArray` | at runtime | the width comes from a file header — the compressed suffix array picks its width from the text length at build time |

They are otherwise interchangeable, and `test_suite.rs` asserts they pack identically at every
width from 1 to 64. Nothing in the workspace benchmarks the two against each other; the const
generic buys const-folding, which is a mechanism, not a measurement.

## The rest of the crate

* `Binary` serialises a bit array's backing words — headerless, and read to EOF. Worth reading the
  contract before embedding one in a larger file.
* `data_to_writer` packs and writes a `Vec<i64>` in chunks, so building a file does not need a
  second full copy of the index in memory.

Both constructors advise transparent huge pages over the allocation **before** they zero it, via
`memory_hints::hugepages`. That ordering is the whole point and is argued there; callers do not
need to repeat it.

## Where it sits

Depends on `memory-hints`. Used by `protein-text` and `sa-index`.

---

Part of the [Unipept Index](../README.md) workspace · full API docs with
`cargo doc -p bitarray --open`
