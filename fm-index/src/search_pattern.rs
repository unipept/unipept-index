use std::error::Error;

pub struct SearchPattern {
    parts: Vec<Vec<u8>>,
    length: usize
}

impl SearchPattern {
    
    pub fn new(pattern: Vec<u8>, parts_amount: usize) -> Result<Self, Box<dyn Error + Send + Sync>> {

        if parts_amount == 0 {
            return Err("parts_amount must be greater than 0".into());
        }

        let total_len = pattern.len();

        if total_len < parts_amount {
            return Err("parts_amount is greater than pattern length".into());
        }

        let base_size = total_len / parts_amount;
        let remainder = total_len % parts_amount;

        let mut parts = Vec::with_capacity(parts_amount);
        let mut start = 0;

        for i in 0..parts_amount {
            let extra = if i < remainder { 1 } else { 0 };
            let end = start + base_size + extra;
            parts.push(pattern[start..end].to_vec());
            start = end;
        }

        let length = total_len;

        Ok(SearchPattern { parts, length })


    }

    pub fn get_part(&self, index: u8, direction_left: bool) -> Vec<u8> {
        let mut part = self.parts.get(index as usize).unwrap().clone();
        if direction_left { part.reverse() };
        part
    }

    pub fn get_part_len(&self, index: u8) -> usize {
        self.parts.get(index as usize).unwrap().len()
    }

    pub fn len(&self) -> usize {
        self.length
    }
}