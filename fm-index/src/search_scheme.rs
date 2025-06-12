use std::error::Error;
use std::fs;
use std::path::Path;

pub struct Search {
    pub order: Vec<u8>,
    pub min_mismatches: Vec<u8>,
    pub max_mismatches: Vec<u8>,
}

pub struct SearchIter<'a> {
    search: &'a Search,
    index: usize,
}

impl<'a> Iterator for SearchIter<'a> {
    type Item = (u8, u8, u8);

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.search.order.len() {
            let i = self.index;
            self.index += 1;
            Some((
                self.search.order[i],
                self.search.min_mismatches[i],
                self.search.max_mismatches[i],
            ))
        } else {
            None
        }
    }
}

impl Search {
    pub fn iter(&self) -> SearchIter {
        SearchIter {
            search: self,
            index: 0,
        }
    }

    pub fn get_direction_left(&self, idx: usize) -> bool {
        if idx == 0 {
            return self.order[1] < self.order[0];
        }

        self.order[idx] < self.order[idx-1]
    }

    pub fn get_upperbound(&self, idx: usize) -> u8 {
        self.max_mismatches[idx]
    }

    pub fn get_lowerbound(&self, idx: usize) -> u8 {
        self.min_mismatches[idx]
    }

    pub fn get_part(&self, idx: usize) -> u8 {
        self.order[idx]
    }

    pub fn len(&self) -> usize {
        self.order.len()
    }

    pub fn validate(&self) -> Result<(), Box<dyn Error>> {
        let len = self.order.len();
        if self.min_mismatches.len() != len || self.max_mismatches.len() != len {
            return Err("Length mismatch in Search fields".into());
        }

        let mut seen = (0, 0);
        let mut prev_min = 0;
        let mut prev_max = 0;

        for (i, &pos) in self.order.iter().enumerate() {
            // Check bounds
            let min = self.min_mismatches[i];
            let max = self.max_mismatches[i];

            if max < min {
                return Err(format!("max_mismatches[{}] < min_mismatches[{}]", i, i).into());
            }

            if i > 0 {
                if min < prev_min {
                    return Err(format!("min_mismatches decreased at step {}", i).into());
                }
                if max < prev_max {
                    return Err(format!("max_mismatches decreased at step {}", i).into());
                }
            }

            // Check contiguity
            if i == 0 {
                seen = (pos, pos);
            } else {
                let (begin, end) = seen;
                if begin == pos + 1 && pos < begin {
                    seen = (pos, end);
                } else if end + 1 == pos && pos < len as u8 {
                    seen = (begin, pos);
                } else {
                    return Err(format!("Part at pos {} does not border seen parts", pos).into());
                }
            }

            prev_min = min;
            prev_max = max;
        }

        Ok(())
    }
}

pub struct SearchScheme {
    searches: Vec<Search>
}

impl SearchScheme {

    pub fn from_file(path: &Path) -> Result<Self, Box<dyn Error>> {
        let content = fs::read_to_string(path)?;
        let mut searches = Vec::new();

        for line in content.lines() {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() != 3 {
                return Err(format!("Invalid line format: {}", line).into());
            }

            let order = Self::parse_braced_numbers(parts[0])?;
            let min_mismatches = Self::parse_braced_numbers(parts[1])?;
            let max_mismatches = Self::parse_braced_numbers(parts[2])?;
            
            searches.push(Search {
                order,
                min_mismatches,
                max_mismatches
            });
        }

        Ok(SearchScheme { searches })
    }

    pub fn validate(&self) -> Result<(), Box<dyn Error>> {
        for (i, search) in self.searches.iter().enumerate() {
            search.validate().map_err(|e| format!("Search {} invalid: {}", i, e))?;
        }

        Ok(())
    }

    fn parse_braced_numbers(s: &str) -> Result<Vec<u8>, Box<dyn Error>> {
        let s = s.trim();
        if !(s.starts_with('{') && s.ends_with('}')) {
            return Err(format!("Invalid format for list: {}", s).into());
        }

        let numbers = s[1..s.len() - 1]
            .split(',')
            .map(|num| num.trim().parse::<u8>())
            .collect::<Result<Vec<_>, _>>()?;

        Ok(numbers)
    }

    pub fn get_parts_amount(&self) -> u8 {
        self.searches[0].order.len() as u8
    }

    pub fn get_search(&self, index: usize) -> &Search {
        &self.searches[index]
    }

}

impl<'a> IntoIterator for &'a SearchScheme {
    type Item = &'a Search;
    type IntoIter = std::slice::Iter<'a, Search>;

    fn into_iter(self) -> Self::IntoIter {
        self.searches.iter()
    }
}

