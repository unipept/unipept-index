use criterion::criterion_main;

mod util;
mod algorithm1;
mod algorithm2;

criterion_main!(algorithm1::benches, algorithm2::benches);
