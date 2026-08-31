use rand::{RngExt, rngs::ThreadRng};

/// Generate a random InterPro annotation.
pub fn generate_ipr(random: &mut ThreadRng) -> String {
    format!("IPR:IPR{:06}", random.random_range(0..999999))
}

/// Generate a random Gene Ontology annotation.
pub fn generate_go(random: &mut ThreadRng) -> String {
    format!("GO:{:07}", random.random_range(0..9999999))
}

/// Generate a random Enzyme Commission annotation.
pub fn generate_ec(random: &mut ThreadRng) -> String {
    format!(
        "EC:{}.{}.{}.{}",
        random.random_range(0..8),
        random.random_range(0..30),
        random.random_range(0..30),
        random.random_range(0..200)
    )
}

/// Generate a random annotation.
pub fn generate_annotation(random: &mut ThreadRng) -> String {
    match random.random_range(0..3) {
        0 => generate_ipr(random),
        1 => generate_go(random),
        2 => generate_ec(random),
        _ => unreachable!()
    }
}
