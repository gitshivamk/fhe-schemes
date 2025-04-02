from util.polynomial import Polynomial
from utils.random_samples import gaussian_distribution, uniform_distribution, ternary_distribution
from utils.crt_function import crt_encode_e
from utils.poly_operations import poly_add, poly_mult


def keygen(pp):
    n, q, error_bound, p1, p2 = pp
    s = ternary_distribution(n)
    a = uniform_distribution(0, q, n)
    e = gaussian_distribution(n, 0, error_bound[0], p2)
    crt_e = crt_encode_e(pp, e)

    a_s = poly_mult(a, s, n, q)
    b = poly_add(a_s, crt_e, q)
    a = [-1 * coeffs for coeffs in a]

    # Convert to Polynomial object for FFT operations
    s = Polynomial(n, s)
    b = Polynomial(n, b)
    a = Polynomial(n, a)

    return (b, a), s
