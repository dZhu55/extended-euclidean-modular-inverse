//! Modular multiplicative inverses and linear congruences via the extended Euclidean algorithm.
//!
//! All arithmetic is performed on `u128`. Inputs wider than 128 bits are not supported;
//! the standard library has no generic integer trait that would let us write this once for
//! every unsigned width without pulling in a third-party crate. We chose `u128` because it
//! is the widest unsigned integer in the standard library and covers essentially every
//! practical modulus (up to 2^128 - 1).

/// Returns `(g, x, y)` such that `a * x + b * y == g` and `g == gcd(a, b)`.
///
/// The coefficients are returned as signed values because the intermediate
/// Bezout coefficients are naturally signed. `g` is unsigned: the gcd is
/// always non-negative.
///
/// The algorithm is iterative rather than recursive. Recursion depth is
/// proportional to the number of Euclidean steps, which for Fibonacci-like
/// inputs is logarithmic but still unbounded in principle; iteration keeps
/// the stack flat.
pub fn extended_gcd(a: u128, b: u128) -> (u128, i128, i128) {
    // We track the Bezout coefficients for `a` and `b` as signed integers.
    // `old_r`/`r` follow the Euclidean remainders; `old_s`/`s` and `old_t`/`t`
    // follow the corresponding coefficients.
    let mut old_r: u128 = a;
    let mut r: u128 = b;
    let mut old_s: i128 = 1;
    let mut s: i128 = 0;
    let mut old_t: i128 = 0;
    let mut t: i128 = 1;

    while r != 0 {
        let q = old_r / r;

        // The quotient can exceed i128::MAX when old_r is near u128::MAX and
        // r is 1, so we compute the coefficient updates in i128 with wrapping
        // checked arithmetic. In practice the coefficients stay small because
        // they shrink alongside the remainders; we use checked math to fail
        // loudly rather than silently wrap.
        let q_i = i128::try_from(q).expect("quotient does not fit in i128");

        let new_r = old_r - q * r;
        old_r = r;
        r = new_r;

        let new_s = old_s
            .checked_sub(q_i.checked_mul(s).expect("Bezout coefficient overflow"))
            .expect("Bezout coefficient overflow");
        old_s = s;
        s = new_s;

        let new_t = old_t
            .checked_sub(q_i.checked_mul(t).expect("Bezout coefficient overflow"))
            .expect("Bezout coefficient overflow");
        old_t = t;
        t = new_t;
    }

    (old_r, old_s, old_t)
}

/// Computes the modular multiplicative inverse of `a` modulo `m`.
///
/// Returns `Some(x)` with `0 <= x < m` and `(a * x) % m == 1` when the
/// inverse exists. Returns `None` when `gcd(a, m) != 1`, including the
/// cases `m == 0` and `m == 1` (mod 1 every value is congruent to 0, so no
/// inverse exists).
///
/// `a` is reduced modulo `m` internally, so `a` larger than `m` is fine.
pub fn modular_inverse(a: u128, m: u128) -> Option<u128> {
    if m <= 1 {
        return None;
    }
    let a = a % m;
    let (g, x, _) = extended_gcd(a, m);
    if g != 1 {
        return None;
    }
    // `x` may be negative; bring it into [0, m).
    let x = x.rem_euclid(m as i128);
    Some(x as u128)
}

/// A single linear congruence `a * x ≡ b (mod m)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LinearCongruence {
    pub a: u128,
    pub b: u128,
    pub m: u128,
}

/// Solves `a * x ≡ b (mod m)`.
///
/// Returns `Some((x0, period))` where `x0` is the smallest non-negative
/// solution and `period` is the modulus of the solution set, i.e. every
/// solution is `x0 + k * period` for integer `k`. `period == m / g` where
/// `g = gcd(a, m)`.
///
/// Returns `None` when `m == 0` (no modular structure) or when
/// `b % g != 0` (the congruence is inconsistent).
///
/// Reduction of `a` and `b` modulo `m` happens inside, so oversized inputs
/// are accepted.
pub fn solve_linear_congruence(a: u128, b: u128, m: u128) -> Option<(u128, u128)> {
    if m == 0 {
        return None;
    }
    let a = a % m;
    let b = b % m;
    let (g, x, _) = extended_gcd(a, m);
    if b % g != 0 {
        return None;
    }
    let m_g = m / g;
    // x is a Bezout coefficient for a, so a*x ≡ g (mod m).
    // Multiply by b/g to get a*(x * b/g) ≡ b (mod m).
    let factor = (b / g) % m_g;
    let x_mod = x.rem_euclid(m as i128) as u128;
    let x0 = (x_mod % m_g).wrapping_mul(factor) % m_g;
    Some((x0, m_g))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extended_gcd_basic() {
        assert_eq!(extended_gcd(240, 46), (2, -9, 47));
        assert_eq!(extended_gcd(46, 240), (2, 47, -9));
        assert_eq!(extended_gcd(0, 5), (5, 0, 1));
        assert_eq!(extended_gcd(5, 0), (5, 1, 0));
        assert_eq!(extended_gcd(0, 0), (0, 1, 0));
    }

    #[test]
    fn extended_gcd_identity_holds() {
        let cases = [
            (240u128, 46u128),
            (46, 240),
            (17, 13),
            (1, 1),
            (0, 7),
            (7, 0),
            (123456789, 987654321),
        ];
        for (a, b) in cases {
            let (g, x, y) = extended_gcd(a, b);
            assert_eq!(g, gcd_u128(a, b), "gcd mismatch for {a}, {b}");
            let lhs = (a as i128)
                .checked_mul(x)
                .and_then(|v| v.checked_add((b as i128).checked_mul(y)?));
            assert_eq!(lhs, Some(g as i128), "Bezout identity failed for {a}, {b}");
        }
    }

    fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
        while b != 0 {
            let r = a % b;
            a = b;
            b = r;
        }
        a
    }

    #[test]
    fn modular_inverse_small() {
        assert_eq!(modular_inverse(3, 11), Some(4));
        assert_eq!(modular_inverse(10, 17), Some(12));
        assert_eq!(modular_inverse(1, 7), Some(1));
        assert_eq!(modular_inverse(6, 7), Some(6));
    }

    #[test]
    fn modular_inverse_reduces_a() {
        // 14 mod 11 == 3, whose inverse is 4.
        assert_eq!(modular_inverse(14, 11), Some(4));
        assert_eq!(modular_inverse(25, 11), Some(4));
    }

    #[test]
    fn modular_inverse_none_when_not_coprime() {
        assert_eq!(modular_inverse(2, 4), None);
        assert_eq!(modular_inverse(0, 5), None);
        assert_eq!(modular_inverse(5, 5), None);
    }

    #[test]
    fn modular_inverse_none_for_trivial_moduli() {
        assert_eq!(modular_inverse(1, 0), None);
        assert_eq!(modular_inverse(1, 1), None);
        assert_eq!(modular_inverse(0, 0), None);
        assert_eq!(modular_inverse(0, 1), None);
    }

    #[test]
    fn modular_inverse_roundtrip() {
        let cases = [(3u128, 11u128), (7, 40), (13, 100), (2, 101), (1, 2)];
        for (a, m) in cases {
            let inv = modular_inverse(a, m).expect("inverse should exist");
            assert!(inv < m, "inverse out of range for {a} mod {m}");
            assert_eq!((a * inv) % m, 1, "roundtrip failed for {a} mod {m}");
        }
    }

    #[test]
    fn modular_inverse_large_modulus() {
        // A large prime modulus where the inverse is known to exist.
        let m: u128 = 170141183460469231731687303715884105727; // 2^127 - 1, a Mersenne prime.
        let a: u128 = 3;
        let inv = modular_inverse(a, m).expect("inverse should exist modulo a prime");
        assert!(inv < m);
        assert_eq!((a * inv) % m, 1);
    }

    #[test]
    fn solve_congruence_simple() {
        // 3x ≡ 1 (mod 11) -> x ≡ 4 (mod 11).
        assert_eq!(solve_linear_congruence(3, 1, 11), Some((4, 11)));
        // 2x ≡ 4 (mod 6) -> x ≡ 2 (mod 3).
        assert_eq!(solve_linear_congruence(2, 4, 6), Some((2, 3)));
    }

    #[test]
    fn solve_congruence_inconsistent() {
        // 2x ≡ 3 (mod 6): gcd(2, 6) = 2, 3 % 2 != 0.
        assert_eq!(solve_linear_congruence(2, 3, 6), None);
        // 4x ≡ 1 (mod 8): gcd(4, 8) = 4.
        assert_eq!(solve_linear_congruence(4, 1, 8), None);
    }

    #[test]
    fn solve_congruence_zero_modulus() {
        assert_eq!(solve_linear_congruence(1, 1, 0), None);
        assert_eq!(solve_linear_congruence(0, 0, 0), None);
    }

    #[test]
    fn solve_congruence_reduces_inputs() {
        // 14x ≡ 8 (mod 6) reduces to 2x ≡ 2 (mod 6) -> x ≡ 1 (mod 3).
        assert_eq!(solve_linear_congruence(14, 8, 6), Some((1, 3)));
    }

    #[test]
    fn solve_congruence_solution_set_verifies() {
        let cases = [
            (3u128, 1u128, 11u128),
            (2, 4, 6),
            (10, 15, 35),
            (7, 0, 13),
            (0, 0, 5),
        ];
        for (a, b, m) in cases {
            if let Some((x0, period)) = solve_linear_congruence(a, b, m) {
                assert!(x0 < period, "x0 out of range for {a}x={b} mod {m}");
                assert_eq!(period, m / gcd_u128(a % m, m));
                // Verify the base solution.
                assert_eq!((a % m) * (x0 % m) % m, b % m, "base solution wrong for {a}x={b} mod {m}");
                // Verify the next representative in the set.
                let x1 = x0 + period;
                assert_eq!((a % m) * (x1 % m) % m, b % m, "shifted solution wrong for {a}x={b} mod {m}");
            }
        }
    }

    #[test]
    fn solve_congruence_identity_when_coprime() {
        // When gcd(a, m) == 1 the solution is unique mod m.
        let cases = [(3u128, 7u128, 11u128), (5, 1, 17), (2, 9, 101)];
        for (a, b, m) in cases {
            let (x0, period) = solve_linear_congruence(a, b, m).expect("coprime should solve");
            assert_eq!(period, m);
            assert_eq!((a * x0) % m, b % m);
        }
    }
}
