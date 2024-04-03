use crate::{CharacterSet, Decode};

static PREFIXES: [&str; 3] = [ "EC:", "GO:", "IPR:IPR" ];

pub fn decode(input: &[u8]) -> String {
    if input.is_empty() {
        return String::new();
    }

    let mut decoded = String::with_capacity(input.len() * 2);

    for &byte in input {
        let (c1, c2) = CharacterSet::decode_pair(byte);

        decoded.push(c1.into());
        if c2 != '$' {
            decoded.push(c2.into());
        }
    }

    let mut result = String::new();
    for (i, annotations) in decoded.split(',').enumerate() {
        if !annotations.is_empty() {
            for annotation in annotations.split(';') {
                result.push_str(PREFIXES[i]);
                result.push_str(annotation);
                result.push(';');
            }
        }
    }

    // Remove the trailing semicolon
    result.pop();

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_empty() {
        assert_eq!(decode(&[]), "")
    }

    #[test]
    fn test_decode_single_ec() {
        assert_eq!(decode(&[ 44, 44, 44, 189, 208 ]), "EC:1.1.1.-")
    }

    #[test]
    fn test_decode_single_go() {
        assert_eq!(decode(&[ 209, 17, 163, 138, 208 ]), "GO:0009279")
    }

    #[test]
    fn test_decode_single_ipr() {
        assert_eq!(decode(&[ 221, 18, 116, 117 ]), "IPR:IPR016364")
    }

    #[test]
    fn test_decode_no_ec() {
        assert_eq!(decode(&[ 209, 17, 163, 138, 209, 39, 71, 94, 17, 153, 39 ]), "GO:0009279;IPR:IPR016364;IPR:IPR008816")
    }

    #[test]
    fn test_decode_no_go() {
        assert_eq!(decode(&[ 44, 44, 44, 190, 44, 60, 44, 141, 209, 39, 71, 80 ]), "EC:1.1.1.-;EC:1.2.1.7;IPR:IPR016364")
    }

    #[test]
    fn test_decode_no_ipr() {
        assert_eq!(decode(&[ 44, 44, 44, 189, 17, 26, 56, 174, 17, 26, 56, 173 ]), "EC:1.1.1.-;GO:0009279;GO:0009279")
    }

    #[test]
    fn test_decode_all() {
        assert_eq!(
            decode(&[ 44, 44, 44, 189, 17, 26, 56, 173, 18, 116, 117, 225, 67, 116, 110, 17, 153, 39 ]),
            "EC:1.1.1.-;GO:0009279;IPR:IPR016364;IPR:IPR032635;IPR:IPR008816"
        )
    }
}
