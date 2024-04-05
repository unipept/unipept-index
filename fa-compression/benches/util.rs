use fa_compression::algorithm1::encode;
use rand::{
    rngs::ThreadRng,
    Rng
};

/// Generate a random InterPro annotation.
pub fn generate_ipr(random: &mut ThreadRng) -> String {
    format!("IPR:IPR{:06}", random.gen_range(0 .. 999999))
}

/// Generate a random Gene Ontology annotation.
pub fn generate_go(random: &mut ThreadRng) -> String {
    format!("GO:{:07}", random.gen_range(0 .. 9999999))
}

/// Generate a random Enzyme Commission annotation.
pub fn generate_ec(random: &mut ThreadRng) -> String {
    format!(
        "EC:{}.{}.{}.{}",
        random.gen_range(0 .. 8),
        random.gen_range(0 .. 30),
        random.gen_range(0 .. 30),
        random.gen_range(0 .. 200)
    )
}

/// Generate a random annotation.
pub fn generate_annotation(random: &mut ThreadRng) -> String {
    match random.gen_range(0 .. 3) {
        0 => generate_ipr(random),
        1 => generate_go(random),
        2 => generate_ec(random),
        _ => unreachable!()
    }
}

/// Generate a random number of decoded annotations.
pub fn generate_decoded_annotations(count: usize) -> String {
    let mut random = rand::thread_rng();

    let mut annotations = String::new();
    for _ in 0 .. count {
        annotations.push_str(&generate_annotation(&mut random));
        annotations.push(';');
    }
    annotations.pop();
    annotations
}

/// Generate a random number of encoded annotations.
pub fn generate_encoded_annotations(count: usize) -> Vec<u8> {
    let mut random = rand::thread_rng();

    let mut annotations = String::new();
    for _ in 0 .. count {
        annotations.push_str(&generate_annotation(&mut random));
        annotations.push(';');
    }
    annotations.pop();

    encode(annotations.as_str())
}
