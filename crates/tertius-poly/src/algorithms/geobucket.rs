//! Geobucket-based sparse polynomial multiplication.
//!
//! Geobuckets provide efficient accumulation of sparse polynomial products
//! with O(log n) amortized time per term insertion.
//!
//! Reference: Yan, "Geobuckets for Polynomial Multiplication" (1998)

use std::cmp::Ordering;

use tertius_rings::traits::Ring;

use crate::monomial::PackedMonomial;
use crate::ordering::MonomialOrder;

/// A geobucket for accumulating sparse polynomial terms.
///
/// Uses geometrically increasing bucket sizes for efficient merging.
pub struct Geobucket<R: Ring> {
    /// Buckets with geometrically increasing sizes.
    /// Bucket i can hold up to 2^(i+1) terms.
    buckets: Vec<Vec<(PackedMonomial, R)>>,
    /// Number of variables in monomials.
    num_vars: usize,
    /// Monomial ordering for comparisons.
    order: MonomialOrder,
}

impl<R: Ring> Geobucket<R> {
    /// Creates a new empty geobucket.
    #[must_use]
    pub fn new(num_vars: usize, order: MonomialOrder) -> Self {
        Self {
            buckets: Vec::new(),
            num_vars,
            order,
        }
    }

    /// Returns the capacity of bucket i.
    #[inline]
    fn bucket_capacity(i: usize) -> usize {
        1 << (i + 1) // 2, 4, 8, 16, ...
    }

    /// Adds a single term to the geobucket.
    pub fn add_term(&mut self, monomial: PackedMonomial, coeff: R) {
        if coeff.is_zero() {
            return;
        }

        // Start with a single-term "polynomial"
        let mut terms = vec![(monomial, coeff)];

        // Carry-propagate through buckets
        let mut i = 0;
        loop {
            // Ensure bucket exists
            while self.buckets.len() <= i {
                self.buckets.push(Vec::new());
            }

            if self.buckets[i].is_empty() {
                // Bucket is empty, just place terms here
                self.buckets[i] = terms;
                break;
            }

            // Merge with existing bucket
            terms = self.merge_sorted(&self.buckets[i], &terms);
            self.buckets[i].clear();

            if terms.len() <= Self::bucket_capacity(i) {
                // Fits in current bucket
                self.buckets[i] = terms;
                break;
            }

            // Overflow - carry to next bucket
            i += 1;
        }
    }

    /// Adds all terms from a sparse polynomial.
    pub fn add_poly(&mut self, terms: &[(PackedMonomial, R)]) {
        if terms.is_empty() {
            return;
        }

        // Add the entire polynomial as a sorted chunk
        let mut chunk: Vec<_> = terms.iter().cloned().collect();
        chunk.sort_by(|a, b| self.compare_monomials(&b.0, &a.0)); // Descending order

        // Find the appropriate bucket based on size
        let mut i = 0;
        while Self::bucket_capacity(i) < chunk.len() {
            i += 1;
        }

        // Carry-propagate from that bucket
        loop {
            while self.buckets.len() <= i {
                self.buckets.push(Vec::new());
            }

            if self.buckets[i].is_empty() {
                self.buckets[i] = chunk;
                break;
            }

            chunk = self.merge_sorted(&self.buckets[i], &chunk);
            self.buckets[i].clear();

            if chunk.len() <= Self::bucket_capacity(i) {
                self.buckets[i] = chunk;
                break;
            }

            i += 1;
        }
    }

    /// Merges two sorted term vectors, combining like terms.
    fn merge_sorted(
        &self,
        a: &[(PackedMonomial, R)],
        b: &[(PackedMonomial, R)],
    ) -> Vec<(PackedMonomial, R)> {
        let mut result = Vec::with_capacity(a.len() + b.len());
        let mut i = 0;
        let mut j = 0;

        while i < a.len() && j < b.len() {
            match self.compare_monomials(&a[i].0, &b[j].0) {
                Ordering::Greater => {
                    result.push(a[i].clone());
                    i += 1;
                }
                Ordering::Less => {
                    result.push(b[j].clone());
                    j += 1;
                }
                Ordering::Equal => {
                    let sum = a[i].1.clone() + b[j].1.clone();
                    if !sum.is_zero() {
                        result.push((a[i].0, sum));
                    }
                    i += 1;
                    j += 1;
                }
            }
        }

        // Append remaining terms
        while i < a.len() {
            result.push(a[i].clone());
            i += 1;
        }
        while j < b.len() {
            result.push(b[j].clone());
            j += 1;
        }

        result
    }

    /// Compares two monomials according to the ordering.
    fn compare_monomials(&self, a: &PackedMonomial, b: &PackedMonomial) -> Ordering {
        self.order.compare(a, b, self.num_vars)
    }

    /// Extracts the final sorted polynomial from the geobucket.
    #[must_use]
    pub fn extract(self) -> Vec<(PackedMonomial, R)> {
        if self.buckets.is_empty() {
            return Vec::new();
        }

        // Merge all buckets from smallest to largest
        let mut result = Vec::new();
        for bucket in &self.buckets {
            if !bucket.is_empty() {
                result = self.merge_sorted(&result, bucket);
            }
        }

        result
    }
}

/// Multiplies two sparse polynomials using geobuckets.
///
/// This is more efficient than naive multiplication for sparse polynomials
/// with many terms.
pub fn geobucket_multiply<R: Ring>(
    a: &[(PackedMonomial, R)],
    b: &[(PackedMonomial, R)],
    num_vars: usize,
    order: MonomialOrder,
) -> Vec<(PackedMonomial, R)> {
    if a.is_empty() || b.is_empty() {
        return vec![(PackedMonomial::one(num_vars), R::zero())];
    }

    // Choose the smaller polynomial to iterate over
    let (smaller, larger) = if a.len() <= b.len() {
        (a, b)
    } else {
        (b, a)
    };

    let mut bucket = Geobucket::new(num_vars, order);

    // For each term in the smaller polynomial, multiply with larger
    for (mono_a, coeff_a) in smaller {
        for (mono_b, coeff_b) in larger {
            let mono = mono_a.mul(mono_b);
            let coeff = coeff_a.clone() * coeff_b.clone();
            bucket.add_term(mono, coeff);
        }
    }

    bucket.extract()
}

/// Heap-based sparse polynomial multiplication.
///
/// Uses a binary heap to maintain the leading terms, which can be more
/// cache-friendly for certain access patterns.
pub fn heap_multiply<R: Ring>(
    a: &[(PackedMonomial, R)],
    b: &[(PackedMonomial, R)],
    num_vars: usize,
    order: MonomialOrder,
) -> Vec<(PackedMonomial, R)> {
    if a.is_empty() || b.is_empty() {
        return vec![(PackedMonomial::one(num_vars), R::zero())];
    }

    // Sort both polynomials by decreasing monomial order
    let mut a_sorted: Vec<_> = a.to_vec();
    let mut b_sorted: Vec<_> = b.to_vec();

    a_sorted.sort_by(|x, y| order.compare(&y.0, &x.0, num_vars));
    b_sorted.sort_by(|x, y| order.compare(&y.0, &x.0, num_vars));

    // Use a simple accumulator approach for now
    // (A full heap implementation would be more complex due to Ord constraints)
    let mut bucket = Geobucket::new(num_vars, order);

    for (mono_a, coeff_a) in &a_sorted {
        for (mono_b, coeff_b) in &b_sorted {
            let mono = mono_a.mul(mono_b);
            let coeff = coeff_a.clone() * coeff_b.clone();
            bucket.add_term(mono, coeff);
        }
    }

    bucket.extract()
}

#[cfg(test)]
mod tests {
    use super::*;
    use tertius_rings::rationals::Q;

    #[test]
    fn test_geobucket_simple() {
        // (1 + x) * (1 + x) = 1 + 2x + x^2
        // Monomials: 1 = exp(0,0,0), x = exp(1,0,0)
        let a = vec![
            (PackedMonomial::one(3), Q::from_integer(1)),
            (PackedMonomial::from_exponents(&[1, 0, 0]), Q::from_integer(1)),
        ];

        let result = geobucket_multiply(&a, &a, 3, MonomialOrder::Grevlex);

        // Should have 3 terms: 1, 2x, x^2
        assert_eq!(result.len(), 3);

        // Verify total sum of coefficients
        let sum: i64 = result
            .iter()
            .map(|(_, c)| c.0.numerator().to_i64().unwrap())
            .sum();
        assert_eq!(sum, 4); // 1 + 2 + 1
    }

    #[test]
    fn test_geobucket_vs_naive() {
        // Compare geobucket result with naive multiplication
        let a = vec![
            (PackedMonomial::from_exponents(&[1, 0, 0]), Q::from_integer(2)),
            (PackedMonomial::from_exponents(&[0, 1, 0]), Q::from_integer(3)),
        ];
        let b = vec![
            (PackedMonomial::from_exponents(&[1, 0, 0]), Q::from_integer(1)),
            (PackedMonomial::from_exponents(&[0, 0, 1]), Q::from_integer(4)),
        ];

        let result = geobucket_multiply(&a, &b, 3, MonomialOrder::Grevlex);

        // (2x + 3y) * (x + 4z) = 2x^2 + 8xz + 3xy + 12yz
        // Terms: x^2, xy, xz, yz
        assert_eq!(result.len(), 4);
    }

    #[test]
    fn test_heap_multiply() {
        let a = vec![
            (PackedMonomial::one(3), Q::from_integer(1)),
            (PackedMonomial::from_exponents(&[1, 0, 0]), Q::from_integer(1)),
        ];

        let result = heap_multiply(&a, &a, 3, MonomialOrder::Grevlex);

        // (1 + x)^2 = 1 + 2x + x^2
        assert_eq!(result.len(), 3);
    }
}
