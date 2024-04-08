use criterion::criterion_group;

use super::util;

mod encode;
mod decode;

criterion_group!(benches, encode::encode_benchmark, decode::decode_benchmark);
