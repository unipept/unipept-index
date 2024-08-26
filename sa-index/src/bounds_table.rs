pub struct BoundsCache<const K: u32> {
    pub bounds: Vec<Option<(usize, usize)>>,

    ascii_array: [usize; 128],
    alphabet: Vec<u8>,
    base: usize
}

impl<const K: u32> BoundsCache<K> {
    pub fn new(alphabet: String) -> BoundsCache<K> {
        let alphabet = alphabet.to_uppercase().as_bytes().to_vec();
        let base = alphabet.len() + 1;

        let mut ascii_array: [usize; 128] = [0; 128];
        for (i, byte) in alphabet.iter().enumerate() {
            // Add 1 to the index, so we can reserve a different value for the 0 index
            ascii_array[*byte as usize] = i + 1;
        }

        BoundsCache {
            bounds: vec![None; base.pow(K)],
            ascii_array,
            alphabet,
            base
        }
    }

    pub fn get_kmer(&self, kmer: &[u8]) -> Option<(usize, usize)> {
        self.bounds.get(self.kmer_to_index(kmer)).cloned()?
    }

    pub fn update_all_kmers(&mut self, kmer: &[u8], bounds: (usize, usize)) {
        let index = self.kmer_to_index(kmer);
        self.bounds[index] = Some(bounds);
    }

    pub fn update_kmer(&mut self, kmer: &[u8], bounds: (usize, usize)) {
        let index = self.kmer_to_index(kmer);
        self.bounds[index] = Some(bounds);
    }

    pub fn index_to_kmer(&self, mut index: usize) -> Vec<u8> {
        let mut kmer = Vec::with_capacity(K as usize);

        for _ in 0..K {
            let modulo = index % self.base;
            if modulo == 0 {
                return kmer.iter().rev().cloned().collect();
            }
            kmer.push(self.alphabet[modulo - 1]);

            index /= self.base;
        }

        kmer.iter().rev().cloned().collect()
    }

    fn kmer_to_index(&self, kmer: &[u8]) -> usize {
        kmer
            .iter()
            .rev()
            .enumerate()
            .map(|(i, n)| self.ascii_array[*n as usize] * self.base.pow(i as u32))
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bounds_cache() {
        let kmer_cache = BoundsCache::<5>::new("ACDEFGHIKLMNPQRSTVWY".to_string());

        assert_eq!(kmer_cache.kmer_to_index("A".as_bytes()), 1);
        assert_eq!(kmer_cache.kmer_to_index("C".as_bytes()), 2);
        assert_eq!(kmer_cache.kmer_to_index("AA".as_bytes()), 22);
        assert_eq!(kmer_cache.kmer_to_index("DEMY".as_bytes()), 29798);
        assert_eq!(kmer_cache.kmer_to_index("DEMYS".as_bytes()), 625774);

        // assert_eq!(kmer_cache.index_to_kmer(1), "A".as_bytes());
        // assert_eq!(kmer_cache.index_to_kmer(2), "C".as_bytes());
        // assert_eq!(kmer_cache.index_to_kmer(22), "AA".as_bytes());
        assert_eq!(kmer_cache.index_to_kmer(29798), "DEMY".as_bytes());
        assert_eq!(kmer_cache.index_to_kmer(625774), "DEMYS".as_bytes());
    }
}
