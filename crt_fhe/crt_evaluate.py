from util.polynomial import Polynomial
from utils.poly_operations import poly_add, poly_mult
from math import ceil, sqrt, floor, log
from utils.poly_operations import base_decompose
from utils import round_precisely


def eval_add(ct1, ct2, q):
    c0_1, c1_1 = ct1
    c0_2, c1_2 = ct2

    c0_coeffs = poly_add(c0_1.coeffs, c0_2.coeffs, q)
    c1_coeffs = poly_add(c1_1.coeffs, c1_2.coeffs, q)

    c0 = Polynomial(len(c0_coeffs), c0_coeffs)
    c1 = Polynomial(len(c1_coeffs), c1_coeffs)

    return c0, c1

def mult_coeff(ct1, ct2, pp):
    n, q, error_bound, p1, p2 = pp

    c0_1, c1_1 = ct1
    c0_2, c1_2 = ct2

    c0 = c0_1.multiply_fft(c0_2)
    c0 = c0.mod(q)

    c1 = c0_1.multiply_fft(c1_2).add(c1_1.multiply_fft(c0_2))
    c1 = c1.mod(q)

    c2 = c1_1.multiply_fft(c1_2)
    c2 = c2.mod(q)

    return c0, c1, c2

def relinearize_v1(c0, c1, c2, relin_key, pp):
    n, q, error_bound, p1, p2 = pp
    num_levels = len(relin_key)

    base = ceil(sqrt(q))
    c2_decomposed = c2.base_decompose(base, num_levels)

    new_c0 = c0
    new_c1 = c1

    for i in range(num_levels):
        new_c0 = new_c0.add(relin_key[i][0].multiply(c2_decomposed[i], q), q)
        new_c1 = new_c1.add(relin_key[i][1].multiply(c2_decomposed[i], q), q)

    return new_c0, new_c1


def relinearize_v2(c_0, c_1, c_2, pp, relin_key, p):
    rlk_0, rlk_1 = relin_key

    c2_k0 = poly_mult(c_2, rlk_0, pp[0])
    c2_k1 = poly_mult(c_2, rlk_1, pp[0])

    c2_0 = round_precisely.precise_round_with_mod(c2_k0, p, pp[1])
    c2_1 = round_precisely.precise_round_with_mod(c2_k1, p, pp[1])

    new_c0 = poly_add(c_0, c2_0, pp[1])
    new_c1 = poly_add(c_1, c2_1, pp[1])

    return new_c0, new_c1


def eval_mult_v1(ct1, ct2, relin_key, pp):
    c0, c1, c2 = mult_coeff(ct1, ct2, pp)

    return relinearize_v1(c0, c1, c2, relin_key, pp)


def eval_mult_v2(ct1, ct2, pp, relin_key, p):
    c_0 = poly_mult(ct1[0], ct2[0], pp[0], pp[1])

    c_1_0 = poly_mult(ct1[0], ct2[1], pp[0], pp[1])
    c_1_1 = poly_mult(ct1[1], ct2[0], pp[0], pp[1])
    c_1 = poly_add(c_1_0, c_1_1, pp[1])

    c_2 = poly_mult(ct1[1], ct2[1], pp[0], pp[1])

    return relinearize_v2(c_0, c_1, c_2, pp, relin_key, p)
