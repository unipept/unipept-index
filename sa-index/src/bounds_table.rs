pub struct BoundsCache<const K: u32> {
    pub bounds: Vec<Option<(usize, usize)>>,

    ascii_array: [usize; 128],
    alphabet: Vec<u8>,
    base: usize
}

impl<const K: u32> BoundsCache<K> {
    pub fn new(alphabet: String) -> BoundsCache<K> {
        let alphabet = alphabet.to_uppercase().as_bytes().to_vec();
        let base = alphabet.len();

        let mut ascii_array: [usize; 128] = [0; 128];
        for (i, byte) in alphabet.iter().enumerate() {
            ascii_array[*byte as usize] = i;
        }

        // 20^1 + 20^2 + 20^3 + ... + 20^(K) = (20^(K + 1) - 20) / 19
        let capacity = (20_u32.pow(K + 1) - 20) / 19;

        BoundsCache {
            bounds: vec![None; capacity as usize],
            ascii_array,
            alphabet,
            base
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
        while offset + self.base.pow(length) <= index {
            offset += self.base.pow(length);
            length += 1;
        }

        let mut kmer = Vec::with_capacity(length as usize);

        let mut index = index - offset;

        for _ in 0..length {
            kmer.push(self.alphabet[index % self.base]);
            index /= self.base;
        }

        kmer.iter().rev().cloned().collect()
    }

    fn kmer_to_index(&self, kmer: &[u8]) -> usize {
        if kmer.len() == 1 {
            return self.ascii_array[kmer[0] as usize];
        }

        let a = kmer
            .iter()
            .rev()
            .enumerate()
            .map(|(i, n)| (self.ascii_array[*n as usize] + 1) * self.base.pow(i as u32))
            .sum::<usize>();

        let b = a - 1;
        b
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bounds_cache() {
        let kmer_cache = BoundsCache::<5>::new("ACDEFGHIKLMNPQRSTVWY".to_string());

        for i in 0..20_usize.pow(5) {
            let kmer = kmer_cache.index_to_kmer(i);
            assert_eq!(kmer_cache.kmer_to_index(&kmer), i);
        }
    }
}
