# binary-traits

![Test](https://img.shields.io/github/actions/workflow/status/unipept/unipept-index/test.yml?logo=github&label=test)
![Codecov](https://img.shields.io/codecov/c/github/unipept/unipept-index?token=IZ75A2FY98&flag=binary-traits&logo=codecov)

The I/O traits every on-disk index structure is written and read through. An index is built once by
`sa-builder` and read back through one of two storage backends, so most structures have one writer
and two readers — and these four traits are that contract.

## The four traits

| trait | direction | implemented by |
|---|---|---|
| `WriteBinary` | serialise | the preloaded type, once per structure |
| `ReadBinary` | deserialise into owned memory | the preloaded backend |
| `ReadBinaryMmap` | map the file and decode fields in place | the mmap backend |
| `LoadIndex` | load by whichever route this concrete type uses | both |

No type implements all four. A structure takes either the owned route (`WriteBinary` +
`ReadBinary` + `LoadIndex`) or the mapped one (`ReadBinaryMmap` + `LoadIndex`). `LoadIndex` is what
lets code that is generic over the backends load an index without knowing which route it took —
and what lets `sa-index`'s tests build all sixteen backend combinations without a single `#[cfg]`.

## Where the formats are documented

At their **writer**. The writer lives with the preloaded type but its output is also consumed by an
mmap reader in a different module, so the two are easy to drift apart; every format is therefore
described on the `impl WriteBinary` block, and each reader points back at it. `sa-builder` names
only the preloaded types for the same reason — one writer per structure, whichever backend reads it
later.

## Where it sits

Depends on nothing. `sa-index`, `protein-metadata`, `protein-text` and `sa-builder` all depend on
it directly, which is what lets them implement these traits for each other's types without a
dependency cycle. Each names them as `binary_traits::…`, so the crate a trait is imported from is
the crate that defines it.

---

Part of the [Unipept Index](../README.md) workspace · full API docs with
`cargo doc -p binary-traits --open`
