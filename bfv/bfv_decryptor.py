from utils.poly_operations import poly_add, poly_mult
from utils.round_precisely import divide_and_round
from utils.poly_operations import coeffs_mod
from util.polynomial import Polynomial

def decrypt(c, sk, params):
    n, q, t, error_bound, delta = params
    c0, c1 = c

    c1_s_coeffs = poly_mult(c1.coeffs, sk.coeffs, n)
    c1_s_poly = Polynomial(n, c1_s_coeffs)

    c0_poly = Polynomial(n, c0.coeffs)
    m_dash_poly = c0_poly.add(c1_s_poly, coeff_modulus=q)

    m_dash = m_dash_poly.coeffs

    m = divide_and_round(m_dash, delta)
    m = coeffs_mod(m, t)

    return m
