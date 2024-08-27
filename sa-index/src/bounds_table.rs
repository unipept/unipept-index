pub struct BoundsCache {
    pub bounds: Vec<Option<(usize, usize)>>,
    pub base: usize,
    pub k: usize,

    ascii_array: [usize; 128],
    powers_array: [usize; 10],
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

        let mut powers_array = [0; 10];
        for i in 0..10 {
            powers_array[i] = base.pow(i as u32);
        }

        // 20^1 + 20^2 + 20^3 + ... + 20^(K) = (20^(K + 1) - 20) / 19
        let capacity = (base.pow(k as u32 + 1) - base) / (base - 1);

        Self {
            bounds: vec![None; capacity],
            ascii_array,
            powers_array,
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

        let mut length = 2;
        let mut offset = self.base;
        while offset + self.powers_array[length] <= index {
            offset += self.powers_array[length];
            length += 1;
        }

        let mut index = index - offset;

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
    use super::*;

    #[test]
    fn test_bounds_cache() {
        let kmer_cache = BoundsCache::new("ACDEFGHIKLMNPQRSTVWY".to_string(), 5);

        for i in 0..20_usize.pow(5) {
            let kmer = kmer_cache.index_to_kmer(i);
            assert_eq!(kmer_cache.kmer_to_index(&kmer), i);
        }
    }
}
