# Unipept Index

![Codecov](https://img.shields.io/codecov/c/github/unipept/unipept-index?token=IZ75A2FY98&logo=Codecov)

The unipept index written entirely in `Rust`. This repository consists of multiple different Rust projects that depend on each other. More information about each project can be found in their respective `README.md` file.

## Installation

> [!NOTE]  
> To build and use the Unipept Index, you need to have Rust installed. If you don't have Rust installed, you can get it from [rust-lang.org](https://www.rust-lang.org/).

Clone this repository by executing the following command:

```bash
git clone https://github.com/unipept/unipept-index.git
```

Next, build everything by executing the following command in the root of the repository.

```bash
cargo build --release
```

After a successful build, two executable binaries are available under `target/release/`:

- `sa-builder` — builds the index from a protein database file
- `sa-server` — serves peptide search queries over HTTP

---

## Utilities

### sa-builder

`sa-builder` constructs a suffix array index from a protein sequence database. It reads a tab-separated input file, builds a (sparse, optionally compressed) suffix array over the concatenated protein sequences, and writes three binary output files: the suffix array itself, the serialised protein metadata, and a suffix-to-protein mapping. These three files are the input required by `sa-server`.

#### Usage

```
sa-builder [OPTIONS] --database-file <DATABASE_FILE> --output-sa <OUTPUT_SA> --output-proteins <OUTPUT_PROTEINS> --output-mapping <OUTPUT_MAPPING>
```

#### Input file format

The database file must be a tab-separated file (TSV) with one protein per line and four columns:

```
<uniprot_accession>  <taxon_id>  <amino_acid_sequence>  <functional_annotations>
```

Functional annotations are a semicolon-separated list of terms such as `GO:0009279;IPR:IPR016364`.

#### Parameters

| Flag | Short | Required | Default | Description |
|------|-------|----------|---------|-------------|
| `--database-file` | `-d` | yes | — | Path to the input TSV protein database file. |
| `--output-sa` | — | yes | — | Output path for the binary suffix array file. |
| `--output-proteins` | — | yes | — | Output path for the binary proteins file. |
| `--output-mapping` | — | yes | — | Output path for the binary suffix-to-protein mapping file. |
| `--sparseness-factor` | `-s` | no | `1` | Sparseness factor *k* for the suffix array. Only every *k*-th suffix is stored, which reduces memory and file size by a factor of *k*. Peptides shorter than *k* amino acids cannot be searched in the resulting index. A value of `1` stores every suffix (no sparseness). Accepted range: 1–255. |
| `--construction-algorithm` | `-a` | no | `lib-sais` | Algorithm used to construct the suffix array. Accepted values: `lib-sais` (default), `lib-div-suf-sort`. See below for details. |
| `--compress-sa` | `-c` | no | `false` | Flag. When set, suffix array values are stored using the minimum number of bits required (derived from the total text length) rather than 64 bits per value. This reduces the output file size at no loss of information. |
| `--mapping-style` | — | no | `bit-vec` | Style of the suffix-to-protein mapping. Accepted values: `bit-vec` (default), `dense`, `sparse`. See below for details. |

#### Construction algorithm

- **`lib-sais`** (default): Uses the [libsais](https://github.com/IlyaGrebnov/libsais) algorithm. Supports native sparse suffix array construction, which is more memory-efficient than building a full array and sampling it afterwards. The maximum effective sparseness that libsais handles natively is 5; for larger sparseness factors the tool applies an additional sampling step on top of libsais output.
- **`lib-div-suf-sort`**: Uses the [libdivsufsort](https://github.com/y-256/libdivsufsort) algorithm. Always builds the full (non-sparse) suffix array first, then applies the sparseness sampling step afterwards. This requires more memory during construction when a sparseness factor greater than 1 is used.

#### Mapping style

The mapping file records which text position belongs to which protein. Three representations are available:

- **`bit-vec`** (default): Stores one bit per text position using a [Rank9](https://github.com/ssteinberg/succinct) rank/select data structure. Separator and terminator characters are marked with a 1-bit; protein characters with a 0-bit. Protein lookup is O(1) via a rank query. Memory usage is approximately 2 bits per character. This is the recommended option.
- **`dense`**: Stores a 32-bit protein index for every text position (4 bytes per character). Lookup is O(1) by direct array access. Higher memory and file size than `bit-vec`, but conceptually simpler.
- **`sparse`**: Stores only the starting text position of each protein (8 bytes per protein). Lookup is O(log m) where *m* is the number of proteins, via binary search. Smallest file size of the three options, but slower lookup.

#### Example

```bash
sa-builder \
  --database-file uniprot_proteins.tsv \
  --output-sa index.sa \
  --output-proteins index.proteins.bin \
  --output-mapping index.mapping \
  --sparseness-factor 3 \
  --construction-algorithm lib-sais \
  --compress-sa \
  --mapping-style bit-vec
```

---

### sa-server

`sa-server` loads the three binary files produced by `sa-builder` and exposes an HTTP API for peptide search. It accepts POST requests containing a list of peptide sequences and returns all matching proteins from the index along with their UniProt accession, taxon ID, and functional annotations.

The server binds to `0.0.0.0:3000` and is ready to accept requests once all index files are loaded.

#### Usage

```
sa-server [OPTIONS] --database-file <DATABASE_FILE> --index-file <INDEX_FILE> --mapping-file <MAPPING_FILE>
```

#### Parameters

| Flag | Short | Required | Default | Description |
|------|-------|----------|---------|-------------|
| `--database-file` | `-d` | yes | — | Path to the binary proteins file (`.proteins.bin`) produced by `sa-builder`. |
| `--index-file` | `-i` | yes | — | Path to the binary suffix array file produced by `sa-builder`. |
| `--mapping-file` | — | yes | — | Path to the binary suffix-to-protein mapping file produced by `sa-builder`. |
| `--mmap` | `-m` | no | `false` | Flag. When set, index files are loaded via memory-mapped I/O instead of being read fully into memory. This makes server startup near-instant because the OS pages data in on demand. The trade-off is that the first queries may be slower while the relevant pages are loaded from disk. Recommended for large indexes where startup time matters. |

#### HTTP API

**`POST /search`**

Accepts a JSON body and returns a JSON array of search results.

Request body fields:

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `peptides` | `string[]` | yes | — | List of peptide sequences to search. |
| `cutoff` | `number` | no | `10000` | Maximum number of suffix array matches processed per peptide. If more matches exist in the index, the result is flagged with `cutoff_used: true` and only the first `cutoff` matches are returned. |
| `equate_il` | `boolean` | no | `false` | When `true`, isoleucine (I) and leucine (L) are treated as equivalent during search. This reflects the fact that these two amino acids are indistinguishable by tryptic mass spectrometry. |
| `tryptic` | `boolean` | no | `false` | When `true`, only tryptic matches are returned. A match is considered tryptic if it starts at the beginning of a protein or directly after a lysine (K) or arginine (R) residue (not followed by proline), and ends at the end of a protein or at such a cleavage site. |

Response: a JSON array where each element corresponds to one peptide that had at least one match. Peptides with no matches are omitted from the response.

```json
[
  {
    "sequence": "MSKIAALLPSV",
    "cutoff_used": false,
    "proteins": [
      {
        "taxon": 9606,
        "uniprot_accession": "P12345",
        "functional_annotations": "GO:0005737;IPR:IPR016364"
      }
    ]
  }
]
```

#### Example

```bash
sa-server \
  --database-file index.proteins.bin \
  --index-file index.sa \
  --mapping-file index.mapping \
  --mmap
```

```bash
curl -X POST http://localhost:3000/search \
  -H "Content-Type: application/json" \
  -d '{
    "peptides": ["MSKIAALLPSV", "ACDEFGHIK"],
    "equate_il": true,
    "tryptic": false,
    "cutoff": 5000
  }'
```
