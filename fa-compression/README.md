# Functional Annotation Compression

![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/unipept/unipept-index/test.yml?logo=github)
![Codecov](https://img.shields.io/codecov/c/github/unipept/unipept-index?token=IZ75A2FY98&flag=fa-compression&logo=codecov)
![Static Badge](https://img.shields.io/badge/doc-rustdoc-blue)

The `fa-compression` library offers compression for Unipept's functional annotation strings. These strings follow a very specific 
format that the compression algorithm will use to achieve a guaranteed minimal compression of **50%** for both very large and very 
small input strings. The compression ratio will often situate around **60-70%**.

The compression algorithm never has to allocate extra memory to build an encoding table or something similar. We can encode each 
string separately. This is particullary useful when all strings have to be encoded/decoded on their own. There is no need to decode 
an entire database to only fetch a single entry.

## Example

```rust
use fa_compression;

fn main() {
    let encoded: Vec<u8> = fa_compression::encode(
        "IPR:IPR016364;EC:1.1.1.-;IPR:IPR032635;GO:0009279;IPR:IPR008816"
    );

    // [ 44, 44, 44, 189, 17, 26, 56, 173, 18, 116, 117, 225, 67, 116, 110, 17, 153, 39 ]
    println!("{:?}", encoded);

    let decoded: String = fa_compression::decode(&encoded);

    // "EC:1.1.1.-;GO:0009279;IPR:IPR016364;IPR:IPR032635;IPR:IPR008816"
    println!("{:?}", decoded);
}
```
