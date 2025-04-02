from util.polynomial import Polynomial
from utils.random_samples import gaussian_distribution, uniform_distribution, ternary_distribution

def keygen(params):
    n, q, t, error_bound, delta = params

    s_coeffs = ternary_distribution(n)
    s_poly = Polynomial(n, s_coeffs)

    a_coeffs = uniform_distribution(0, q, n)
    a_poly = Polynomial(n, a_coeffs)

    e_coeffs = gaussian_distribution(n, 0, error_bound[0])
    e_poly = Polynomial(n, e_coeffs)

    a_s_poly = a_poly.multiply(s_poly, q)
    b_poly = a_s_poly.add(e_poly, coeff_modulus=q)
    b_poly.coeffs = [-coeff for coeff in b_poly.coeffs]

    return (b_poly, a_poly), s_poly
