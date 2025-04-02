from utils.random_samples import gaussian_distribution, ternary_distribution
from utils.poly_operations import poly_add, poly_mult
from util.polynomial import Polynomial

def encrypt(m, params, pk):
    n, q, t, error_bound, delta = params

    # Generate random and error polynomials
    u_coeffs = ternary_distribution(n)
    e1_coeffs = gaussian_distribution(n, 0, error_bound[0])
    e2_coeffs = gaussian_distribution(n, 0, error_bound[1])

    b_poly, a_poly = pk
    delta_m_coeffs = [coeff * delta for coeff in m]

    # Perform polynomial operations using Polynomial class
    u_poly = Polynomial(n, u_coeffs)
    delta_m_poly = Polynomial(n, delta_m_coeffs)
    e1_poly = Polynomial(n, e1_coeffs)
    e2_poly = Polynomial(n, e2_coeffs)

    au_poly = a_poly.multiply(u_poly, coeff_modulus=q)
    bu_poly = b_poly.multiply(u_poly, coeff_modulus=q)

    c0_poly = delta_m_poly.add(bu_poly, coeff_modulus=q).add(e1_poly, coeff_modulus=q)
    c1_poly = au_poly.add(e2_poly, coeff_modulus=q)

    return c0_poly, c1_poly
