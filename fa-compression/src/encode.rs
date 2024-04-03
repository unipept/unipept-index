use crate::{CharacterSet, Encode};

pub fn encode(input: &str) -> Vec<u8> {
    if input.is_empty() {
        return Vec::new();
    }

    let mut interpros: Vec<&str> = Vec::new();
    let mut gos: Vec<&str> = Vec::new();
    let mut ecs: Vec<&str> = Vec::new();

    // If we can make sure the input is sorted, we can avoid the sorting step
    for annotation in input.split(';') {
        if annotation.starts_with("IPR") {
            interpros.push(&annotation[7..]);
        } else if annotation.starts_with("GO") {
            gos.push(&annotation[3..]);
        } else if annotation.starts_with("EC") {
            ecs.push(&annotation[3..]);
        }
    }

    let result = format!("{},{},{}", ecs.join(";"), gos.join(";"), interpros.join(";"));

    let mut encoded: Vec<u8> = Vec::new();

    let mut iter = result.as_bytes().chunks_exact(2);
    while let Some([ b1, b2 ]) = iter.next() {
        let c1 = CharacterSet::encode(*b1);
        let c2 = CharacterSet::encode(*b2);
        encoded.push(c1 | c2);
    }

    let remainder = iter.remainder();
    if !remainder.is_empty() {
        let c1 = CharacterSet::encode(remainder[0]);
        encoded.push(c1 | CharacterSet::EMPTY);
    }

    encoded
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_empty() {
        assert_eq!(encode(""), vec![])
    }

    #[test]
    fn test_encode_single_ec() {
        assert_eq!(encode("EC:1.1.1.-"), vec![ 44, 44, 44, 189, 208 ])
    }

    #[test]
    fn test_encode_single_go() {
        assert_eq!(encode("GO:0009279"), vec![ 209, 17, 163, 138, 208 ])
    }

    #[test]
    fn test_encode_single_ipr() {
        assert_eq!(encode("IPR:IPR016364"), vec![ 221, 18, 116, 117 ])
    }

    #[test]
    fn test_encode_no_ec() {
        assert_eq!(encode("IPR:IPR016364;GO:0009279;IPR:IPR008816"), vec![ 209, 17, 163, 138, 209, 39, 71, 94, 17, 153, 39 ])
    }

    #[test]
    fn test_encode_no_go() {
        assert_eq!(encode("IPR:IPR016364;EC:1.1.1.-;EC:1.2.1.7"), vec![ 44, 44, 44, 190, 44, 60, 44, 141, 209, 39, 71, 80 ])
    }

    #[test]
    fn test_encode_no_ipr() {
        assert_eq!(encode("EC:1.1.1.-;GO:0009279;GO:0009279"), vec![ 44, 44, 44, 189, 17, 26, 56, 174, 17, 26, 56, 173 ])
    }

    #[test]
    fn test_encode_all() {
        assert_eq!(
            encode("IPR:IPR016364;EC:1.1.1.-;IPR:IPR032635;GO:0009279;IPR:IPR008816"),
            vec![ 44, 44, 44, 189, 17, 26, 56, 173, 18, 116, 117, 225, 67, 116, 110, 17, 153, 39 ]
        )
    }
}
