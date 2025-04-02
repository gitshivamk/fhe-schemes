from util.polynomial import Polynomial
from utils.random_samples import gaussian_distribution, ternary_distribution
from utils.crt_function import crt_encode_e, crt_encode_m
from utils.poly_operations import poly_add, poly_mult


def encrypt(m, pp, pk):
    n, q, error_bound, p1, p2 = pp
    u = ternary_distribution(n)
    e1 = gaussian_distribution(n, 0, error_bound[0], p2)
    e2 = gaussian_distribution(n, 0, error_bound[1], p2)

    # CRT encoding
    crt_m = crt_encode_m(pp, m)
    crt_e1 = crt_encode_e(pp, e1)
    crt_e2 = crt_encode_e(pp, e2)
    b, a = pk[0].coeffs, pk[1].coeffs

    au = poly_mult(a, u, n)
    bu = poly_mult(b, u, n)
    bu_e = poly_add(bu, crt_e1)

    c0 = poly_add(bu_e, crt_m, q)
    c1 = poly_add(au, crt_e2, q)

    c0 = Polynomial(n, c0)
    c1 = Polynomial(n, c1)

    return c0, c1
