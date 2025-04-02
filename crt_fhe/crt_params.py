from math import gcd


def public_params(lam):
    n = int(input("Enter degree of polynomial (in powers of 2): "))
    # n = 4
    q = 2 ** lam
    p1 = 16
    p2 = 3

    # Standard deviations for error bounds
    std1 = 1
    std2 = 2
    error_bound = (std1, std2)

    # Check if p1 and p2 are relatively prime
    assert gcd(p1, p2) == 1, "p1 and p2 should be pairwise relatively prime"
    # Ensure p1 * p2 is less than q
    assert p1 * p2 < q, "p1 * p2 should be less than q"

    pp = (n, q, error_bound, p1, p2)

    return pp
