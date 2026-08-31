# fa-compression

![CI](https://img.shields.io/github/actions/workflow/status/unipept/unipept-index/ci.yml?logo=github&label=ci)

Compression for the functional annotations attached to every protein — the `EC:`, `GO:` and
`IPR:IPR` strings the Unipept index builder emits. Annotations are stored encoded and decoded once
per reported search result, so how cheaply a single entry decodes matters more than the ratio.

## Two algorithms

| | `algorithm1` | `algorithm2` |
|---|---|---|
| Shared state | none | a `CompressionTable` over the whole database |
| Ratio | 68-71% (at least 50% even on one short annotation) | ~76% |
| Decode one entry in isolation | yes | no — needs the exact table that produced it |
| Encoding cost | one pass | a linear scan of the table, per annotation |
| Ceiling | none | 2^24 distinct annotations (3-byte indices) |

**`algorithm1` is what the index and the server use**, because a protein hit's annotations have to
be decodable on their own. It exploits the fact that Unipept annotations draw on only sixteen
characters: `encode` strips the prefixes, groups what is left into three `,`-separated sections,
and packs two characters per byte; `decode` puts the prefixes back.

`algorithm2` is the higher-ratio alternative for cases where a whole-database table is acceptable.

## Decoding without allocating

`algorithm1` has three entry points, all the same single pass over the input:

* `decode()` — allocates and returns a `String`
* `decode_into()` — appends to a buffer you already have
* `decoded()` — decodes lazily into a formatter, so a serialiser never materialises the string

The third exists for the search path, which has no consumer here yet: a large response can hold
millions of hits, and owning each one's annotations would cost an allocation per hit purely so
`serde` could copy out of it a moment later.

## Benchmarks

```bash
cargo bench -p fa-compression
```

Criterion benches for both algorithms, on the default `bench` profile (which inherits `release`).

## Where it sits

Depends on nothing at runtime. Used by `sa-mappings`.

---

Part of the [Unipept Index](../README.md) workspace · full API docs with
`cargo doc -p fa-compression --open`
