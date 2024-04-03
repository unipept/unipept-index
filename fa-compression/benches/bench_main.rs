use criterion::{black_box, criterion_group, criterion_main};
use fa_compression::{decode::decode, encode::encode};

mod util;

fn encode_benchmark(c: &mut criterion::Criterion) {
    c.bench_function("encode", |b| {
        b.iter_batched(
            || util::generate_decoded_annotations(100), 
            |annotations| black_box(encode(annotations.as_str())), 
            criterion::BatchSize::SmallInput
        )
    });
}

fn decode_benchmark(c: &mut criterion::Criterion) {
    c.bench_function("decode", |b| {
        b.iter_batched(
            || util::generate_encoded_annotations(100), 
            |annotations| black_box(decode(annotations.as_slice())), 
            criterion::BatchSize::SmallInput
        )
    });
}

criterion_group!(benches, encode_benchmark, decode_benchmark);
criterion_main!(benches);
