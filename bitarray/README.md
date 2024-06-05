# Bitarray

![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/unipept/unipept-index/test.yml?logo=github)
![Codecov](https://img.shields.io/codecov/c/github/unipept/unipept-index?token=IZ75A2FY98&flag=bitarray&logo=codecov)
![Static Badge](https://img.shields.io/badge/doc-rustdoc-blue)

The `bitarray` offers a special array where each item is represented by a specified amount of bits (smaller than 64). The bitarray uses a pre-alocated vector and allows you to `set` or `get` a value from the array.

## Example

```rust
use bitarray;

fn main() {
    let bitarray = BitArray::<40>::with_capacity(4);

    bitarray.set(0, 0b0001110011111010110001000111111100110010);

    assert_eq!(bitarray.get(0), 0b0001110011111010110001000111111100110010);
}
```
