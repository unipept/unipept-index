Suffix Array Builder
====================

A rust implementation to build large generalized suffix arrays.

# Usage

```plain
Build a (sparse, compressed) suffix array from the given text

Usage: sa-builder [OPTIONS] --database-file <DATABASE_FILE> --taxonomy <TAXONOMY> --output <OUTPUT>

Options:
  -d, --database-file <DATABASE_FILE>
          File with the proteins used to build the suffix tree. All the proteins are expected to be concatenated using a hashtag `#`
  -t, --taxonomy <TAXONOMY>
          The taxonomy to be used as a tsv file. This is a preprocessed version of the NCBI taxonomy
  -o, --output <OUTPUT>
          Output location where to store the suffix array
  -s, --sparseness-factor <SPARSENESS_FACTOR>
          The sparseness_factor used on the suffix array (default value 1, which means every value in the SA is used) [default: 1]
  -a, --construction-algorithm <CONSTRUCTION_ALGORITHM>
          The algorithm used to construct the suffix array (default value LibSais) [default: lib-sais] [possible values: lib-div-suf-sort, lib-sais]
  -c, --compress-sa
          If the suffix array should be compressed (default value true)
  -h, --help
          Print help
```
