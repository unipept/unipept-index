use std::ops::{BitOr, Shl};

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
            _ => panic!("Invalid character"),   
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
            _ => panic!("Invalid character"),   
        }
    }
}

impl Into<char> for CharacterSet {
    fn into(self) -> char {
        match self {
            CharacterSet::EMPTY => '$',
            CharacterSet::ZERO => '0',
            CharacterSet::ONE => '1',
            CharacterSet::TWO => '2',
            CharacterSet::THREE => '3',
            CharacterSet::FOUR => '4',
            CharacterSet::FIVE => '5',
            CharacterSet::SIX => '6',
            CharacterSet::SEVEN => '7',
            CharacterSet::EIGHT => '8',
            CharacterSet::NINE => '9',
            CharacterSet::DASH => '-',
            CharacterSet::POINT => '.',
            CharacterSet::COMMA => ',',
            CharacterSet::SEMICOLON => ';',
        }
    }
}

impl BitOr for CharacterSet {
    type Output = u8;

    fn bitor(self, rhs: Self) -> Self::Output {
        (self << 4) | rhs as u8
    }
}

impl Shl<u8> for CharacterSet {
    type Output = u8;

    fn shl(self, rhs: u8) -> Self::Output {
        (self as u8) << rhs
    }
}
