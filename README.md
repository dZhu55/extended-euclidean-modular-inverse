# extended-euclidean-modular-inverse

Computes modular multiplicative inverses and solves linear congruences over `u128` via the extended Euclidean algorithm.

```rust
use extended_euclidean_modular_inverse::{extended_gcd, modular_inverse, solve_linear_congruence};

// Bezout coefficients: 240 * (-9) + 46 * 47 == 2.
assert_eq!(extended_gcd(240, 46), (2, -9, 47));

// 3 * 4 == 12 == 1 (mod 11).
assert_eq!(modular_inverse(3, 11), Some(4));

// 2x ≡ 4 (mod 6) has solutions x ≡ 2 (mod 3).
assert_eq!(solve_linear_congruence(2, 4, 6), Some((2, 3)));

// No inverse when gcd(a, m) != 1.
assert_eq!(modular_inverse(2, 4), None);
```

## Why this exists

Modular inverses show up constantly in number theory code, hash-table sizing, and
cryptographic primitives. The standard library has no integer type that is generic over
unsigned widths, so a library that wants to be dependency-free has to pick one width. We
chose `u128`: it is the widest unsigned integer in the standard library, it covers every
modulus a caller is likely to hand us, and it keeps the implementation to a single concrete
monomorphisation instead of a tower of macros.

The trade-off is explicit: inputs wider than 128 bits are not representable, and the Bezout
coefficients are returned as `i128` because the extended Euclidean algorithm naturally
produces signed values. Callers who need arbitrary precision should reach for a bignum
crate; this library deliberately does not.

## The awkward edge

The extended Euclidean algorithm returns signed Bezout coefficients, but modular inverses
and congruence solutions are conventionally reported in `[0, m)`. The conversion is not
just `x % m` — in Rust, `%` on a negative `i128` yields a negative result. We use
`rem_euclid` so the returned inverse is always in range, and we test that path explicitly
with inputs whose raw coefficient is negative (`modular_inverse(3, 11)` is the smallest such
case).

A second edge is `m == 1`. Modulo 1 every integer is congruent to 0, so `1` and `0` are
the same residue and no element has a multiplicative inverse. Both `modular_inverse` and
`solve_linear_congruence` treat `m == 1` as having no solution rather than returning a
misleading `0`.
