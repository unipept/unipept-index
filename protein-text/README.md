# protein-text

![CI](https://img.shields.io/github/actions/workflow/status/unipept/unipept-index/ci.yml?logo=github&label=ci)

The concatenated protein text, packed at 5 bits per residue. Every protein sequence in the database
is joined into one text, separated by `-` and terminated by `$`. The suffix array indexes positions
in *this* text, so essentially every search operation reads it — comparing a candidate suffix
against the query touches one residue per character compared. **It is the hottest data structure in
the index.**

The alphabet is 25 amino-acid letters plus the two delimiters, so a residue needs 5 bits rather
than 8: packed, the text is five eighths the size it would be at one byte per residue. The saving
scales with the database, and it is worth more than the unpacking costs.

## Two backends

Both are always compiled; callers pick by naming a type, and everything here that touches the text
is generic over `ProteinTextBackend`.

| type | where the bytes live |
|---|---|
| `InMemoryProteinText` | the text in owned RAM. Faster per access, but the process pays the full resident size |
| `MmapBackedProteinText` | decoded straight out of a memory mapping, so the kernel decides what stays resident |

`preloaded` owns the `WriteBinary` implementation that produces the file both backends read, which
is why `sa-builder` needs only that half.

Storing this structure differently from the rest of the index is worth a knob of its own —
`sa-server`'s `preloaded-text`. The text is the hottest structure in the index while the protein
metadata sharing its file is the one that grows most when preloaded, so the best place for one is
not the best place for the other. On the full UniProt index, `mmap,preloaded-text,preloaded-proteins`
is one of the two fastest configurations measured; see the [workspace README](../README.md).

## The alphabet is the on-disk format

`BIT5_TO_CHAR` maps a 5-bit code to an ASCII residue, and the inverse table is built from it at
load time — so that one array is the single definition of both the alphabet and the encoding, and
changing its order changes the file format.

It has 27 entries but is indexed by a 5-bit value, so codes 27..=31 are out of bounds. No encoder
can emit them, but a corrupt or truncated index file can, which would panic. A known limitation
with no tracking issue filed; padding the table to 32 entries would remove both the panic and the
bounds check from the hot path.

`bit_array_byte_size` is the other shared piece: every reader of the format goes through it to
relate a declared text length to a byte count. It is fallible because the length it is handed comes
straight out of an untrusted file header, and a validator must not fail on the very input it exists
to reject.

## Where it sits

Depends on `bitarray`, `binary-traits` and `memory-hints`. Used by `protein-metadata` (the text and
the metadata share one `proteins.bin`), `sa-index`, `sa-builder` and `sa-server`.

---

Part of the [Unipept Index](../README.md) workspace · full API docs with
`cargo doc -p protein-text --open`
