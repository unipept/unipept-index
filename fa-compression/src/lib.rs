use std::ops::BitOr;

pub mod encode;
pub mod decode;

pub trait Encode {
    fn encode(value: u8) -> CharacterSet;
}

pub trait Decode {
    fn decode(value: u8) -> char;
    fn decode_pair(value: u8) -> (char, char) {
        ( Self::decode(value >> 4), Self::decode(value & 0b1111) )
    }
}

#[repr(u8)]
#[cfg_attr(test, derive(Clone, Copy))]
#[derive(PartialEq, Eq, Debug)]
pub enum CharacterSet {
    EMPTY,

    ZERO,
    ONE,
    TWO,
    THREE,
    FOUR,
    FIVE,
    SIX,
    SEVEN,
    EIGHT,
    NINE,

    DASH,
    POINT,
    COMMA,
    SEMICOLON,
}

impl Encode for CharacterSet {
    fn encode(value: u8) -> CharacterSet {
        match value {
            b'$' => CharacterSet::EMPTY,
            b'0' => CharacterSet::ZERO,
            b'1' => CharacterSet::ONE,
            b'2' => CharacterSet::TWO,
            b'3' => CharacterSet::THREE,
            b'4' => CharacterSet::FOUR,
            b'5' => CharacterSet::FIVE,
            b'6' => CharacterSet::SIX,
            b'7' => CharacterSet::SEVEN,
            b'8' => CharacterSet::EIGHT,
            b'9' => CharacterSet::NINE,
            b'-' => CharacterSet::DASH,
            b'.' => CharacterSet::POINT,
            b',' => CharacterSet::COMMA,
            b';' => CharacterSet::SEMICOLON,
            _ => panic!("Invalid character")
        }
    }
}

impl Decode for CharacterSet {
    fn decode(value: u8) -> char {
        match value {
            0 => '$',
            1 => '0',
            2 => '1',
            3 => '2',
            4 => '3',
            5 => '4',
            6 => '5',
            7 => '6',
            8 => '7',
            9 => '8',
            10 => '9',
            11 => '-',
            12 => '.',
            13 => ',',
            14 => ';',
            _ => panic!("Invalid character")
        }
    }
}

impl BitOr for CharacterSet {
    type Output = u8;

    fn bitor(self, rhs: Self) -> Self::Output {
        ((self as u8) << 4) | rhs as u8
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static CHARACTERS: [u8; 15] = [
        b'$', b'0', b'1', b'2', b'3', b'4', b'5', b'6', b'7', b'8', b'9', b'-', b'.', b',', b';'
    ];

    static CHARACTER_SETS: [CharacterSet; 15] = [
        CharacterSet::EMPTY,
        CharacterSet::ZERO,
        CharacterSet::ONE,
        CharacterSet::TWO,
        CharacterSet::THREE,
        CharacterSet::FOUR,
        CharacterSet::FIVE,
        CharacterSet::SIX,
        CharacterSet::SEVEN,
        CharacterSet::EIGHT,
        CharacterSet::NINE,
        CharacterSet::DASH,
        CharacterSet::POINT,
        CharacterSet::COMMA,
        CharacterSet::SEMICOLON
    ];

    #[test]
    fn test_or() {
        for i in 0..CHARACTERS.len() {
            for j in 0..CHARACTERS.len() {
                assert_eq!(CHARACTER_SETS[i] | CHARACTER_SETS[j], ((i as u8) << 4) | (j as u8));
            }
        }
    }

    #[test]
    fn test_encode() {
        for i in 0..CHARACTERS.len() {
            assert_eq!(CHARACTER_SETS[i], CharacterSet::encode(CHARACTERS[i]));
        }
    }

    #[test]
    fn test_decode() {
        for (i, c) in CHARACTERS.iter().enumerate() {
            assert_eq!(CharacterSet::decode(i as u8), *c as char);
        }
    }

    #[test]
    fn test_decode_pair() {
        for (i1, c1) in CHARACTERS.iter().enumerate() {
            for (i2, c2) in CHARACTERS.iter().enumerate() {
                let encoded = CharacterSet::encode(*c1) | CharacterSet::encode(*c2);
                assert_eq!(CharacterSet::decode_pair(encoded), (CharacterSet::decode(i1 as u8), CharacterSet::decode(i2 as u8)));
            }
        }
    }
}
