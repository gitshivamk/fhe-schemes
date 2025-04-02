from utils.poly_operations import poly_add, poly_mult
from utils.poly_operations import coeffs_mod


def decrypt(c, sk, pp):
    n, q, error_bound, p1, p2 = pp
    c0, c1 = c
    c0 = c0.coeffs
    c1 = c1.coeffs
    sk = sk.coeffs

    c1_s = poly_mult(c1, sk, n)
    m = poly_add(c0, c1_s, q)

    m = coeffs_mod(m, p1)

    return m