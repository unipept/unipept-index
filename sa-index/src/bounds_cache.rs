pub struct BoundsCache {
    pub bounds: Vec<Option<(usize, usize)>>,
    pub base: usize,
    pub k: usize,

    ascii_array: [usize; 128],
    powers_array: [usize; 10],
    offsets_array: [usize; 10],
    alphabet: Vec<u8>
}

impl BoundsCache {
    pub fn new(alphabet: String, k: usize) -> BoundsCache {
        assert!(k < 10, "K must be less than 10");

        let alphabet = alphabet.to_uppercase().as_bytes().to_vec();
        let base = alphabet.len();

        let mut ascii_array: [usize; 128] = [0; 128];
        for (i, byte) in alphabet.iter().enumerate() {
            ascii_array[*byte as usize] = i;
        }
        //ascii_array[b'L' as usize] = ascii_array[b'I' as usize]; // I and L are treated as the same amino acid

        let mut powers_array = [0; 10];
        for i in 0..10 {
            powers_array[i] = base.pow(i as u32);
        }

        let mut offsets_array = [0; 10];
        for i in 2..10 {
            offsets_array[i] = offsets_array[i - 1] + powers_array[i - 1];
        }

        // 20^1 + 20^2 + 20^3 + ... + 20^(K) = (20^(K + 1) - 20) / 19
        let capacity = (base.pow(k as u32 + 1) - base) / (base - 1);

        Self {
            bounds: vec![None; capacity],
            ascii_array,
            powers_array,
            offsets_array,
            alphabet,
            base,
            k
        }
    }

    pub fn get_kmer(&self, kmer: &[u8]) -> Option<(usize, usize)> {
        self.bounds.get(self.kmer_to_index(kmer)).cloned()?
    }

    pub fn update_kmer(&mut self, kmer: &[u8], bounds: (usize, usize)) {
        let index = self.kmer_to_index(kmer);
        self.bounds[index] = Some(bounds);
    }

    pub fn index_to_kmer(&self, mut index: usize) -> Vec<u8> {
        if index < self.base {
            return vec![self.alphabet[index]];
        }

        // Calculate the length of the kmer
        let mut length = 2;
        while self.offsets_array[length + 1] <= index {
            length += 1;
        }

        // Calculate the offset of the kmer
        let offset = self.offsets_array[length];

        // Translate the index to be inside the bounds [0, 20^k)
        index -= offset;

        // Basic conversion from base 10 to base `length`
        let mut kmer = vec![0; length];
        for i in 0..length {
            kmer[length - i - 1] = self.alphabet[index % self.base];
            index /= self.base;
        }

        kmer
    }

    fn kmer_to_index(&self, kmer: &[u8]) -> usize {
        if kmer.len() == 1 {
            return self.ascii_array[kmer[0] as usize];
        }

        let mut result = 0;
        for i in 0..kmer.len() {
            result += (self.ascii_array[kmer[i] as usize] + 1) * self.powers_array[kmer.len() - i - 1];
        }

        result - 1
    }
}

#[cfg(test)]
mod tests {
    use std::str::from_utf8;
    use super::*;

    #[test]
    fn test_bounds_cache() {
        let kmer_cache = BoundsCache::new("ACDEFGHIKLMNPQRSTVWY".to_string(), 5);

        for i in 0..40 {
            let kmer = kmer_cache.index_to_kmer(i);
            eprintln!("{} -> {:?} -> {:?}", i, from_utf8(&kmer).unwrap(), kmer_cache.kmer_to_index(&kmer));
        }

        for i in 0..20_usize.pow(5) {
            let kmer = kmer_cache.index_to_kmer(i);

            if kmer.contains(&b'L') {
                let equated_kmer = kmer.iter().map(|&c| if c == b'L' { b'I' } else { c }).collect::<Vec<u8>>();
                assert_eq!(kmer_cache.kmer_to_index(&kmer), kmer_cache.kmer_to_index(&equated_kmer));
                continue;
            }

            assert_eq!(kmer_cache.kmer_to_index(&kmer), i);
        }
    }
}
