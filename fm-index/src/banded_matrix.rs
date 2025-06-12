use std::ops::{Index, IndexMut};

pub struct BandedMatrix {
    matrix: Vec<u8>,
    width: u8,
    rows: usize,
    columns: usize,
    col_per_row: u8
}

impl Index<(usize, usize)> for BandedMatrix {
    type Output = u8;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (i, j) = index;
        &self.matrix[i * self.col_per_row as usize + j - i + self.width as usize]
    }
}

impl IndexMut<(usize, usize)> for BandedMatrix {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let (i, j) = index;
        &mut self.matrix[i * self.col_per_row as usize + j - i + self.width as usize]
    }
}

impl BandedMatrix {

    pub fn new(pattern_size: usize, width: u8, start_value: u8) -> Self {

        let columns = pattern_size + 1;
        let rows = pattern_size + width as usize + 1;
        let col_per_row = (2 * width + 1) + 2;
        let matrix = vec![0; rows * col_per_row as usize];

        let mut banded_matrix = Self { matrix, width, rows, columns, col_per_row };

        banded_matrix.initialize_matrix(start_value);

        banded_matrix
    }

    fn initialize_matrix(&mut self, start_value: u8) {

        let width = self.width;
        let width_usize = width as usize;

        // initialize the top row and leftmost column
        for i in 0..=(width+1) {
            let index = i as usize;
            self[(0, index)] = i + start_value;
            self[(index, 0)] = i + start_value;
        }

        // set max elements at sides
        // first the elements on rows [1, width]
        for i in 1..=width {
            let index = i as usize;
            // right of band
            self[(index, index + width_usize + 1)] = width + 1 + start_value;
        }

        // then the elements on rows [width + 1, x]
        if self.columns > width_usize {
            for i in (width_usize+1)..(self.columns-width_usize-1) {
                let index = i as usize;
                // right of band
                self[(index, index + width_usize + 1)] = width + 1 + start_value;
                // left of band
                self[(index, index - (width_usize + 1))] = width + 1 + start_value;
            } 
        }

        // finally the elements on the final rows
        let mut start = width_usize + 1;
        if self.columns > width_usize && self.columns - (width_usize + 1) > start {
            start = self.columns - (width_usize + 1);
        }
        for i in start..self.rows {
            // left of band
            self[(i, i - (width_usize + 1))] = width + 1 + start_value;
        }

    }

    
    fn get_first_column(&self, row: usize) -> usize {
        // leftmost cell of band
        if self.width as usize >= row {
            return 1;
        }
        row - self.width as usize
    }
    
    fn get_last_column(&self, row: usize) -> usize {
        // rightmost cell of band
        return std::cmp::min(self.columns - 1, self.width as usize + row);
    }

    pub fn update_matrix_cell(&mut self, not_match: bool, row: usize, column: usize) -> u8 {
        let match_score = if not_match { 1 } else { 0 };
        let diag = self[(row-1, column-1)] + match_score;
        let left = self[(row, column-1)] + 1;
        let top = self[(row-1, column)] + 1;
        let result = std::cmp::min(diag, std::cmp::min(left, top));
        
        self[(row, column)] = result;
        result
    }


    pub fn update_matrix_row(&mut self, pattern: &Vec<u8>, row: usize, c: u8) -> u8 {
        // Handle the case where row is not present in the matrix
        if row >= self.rows {
            // if the row exceeds the rows of the matrix return maximal value
            return u8::MAX;
        }

        let mut minimum: u8 = u8::MAX;
        for i in self.get_first_column(row)..=self.get_last_column(row) {
            let not_match: bool = pattern[i-1] != c;
            let result: u8 = self.update_matrix_cell(not_match, row, i);
            
            if result < minimum {
                minimum = result;
            }
        }

        minimum

    }

    pub fn is_final_column(&self, row: usize) -> bool {
        self.get_last_column(row) == self.columns - 1
    }

    pub fn get_value_in_final_column(&self, row: usize) -> u8 {
        self[(row, self.columns - 1)]
    }
}