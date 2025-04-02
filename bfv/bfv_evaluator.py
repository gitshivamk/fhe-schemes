from utils.poly_operations import poly_add
from util.polynomial import Polynomial
from math import ceil, sqrt


def eval_add(ct1, ct2, q):
    c0_1, c1_1 = ct1
    c0_2, c1_2 = ct2

    c0_coeffs = poly_add(c0_1.coeffs, c0_2.coeffs, q)
    c1_coeffs = poly_add(c1_1.coeffs, c1_2.coeffs, q)

    c0 = Polynomial(len(c0_coeffs), c0_coeffs)
    c1 = Polynomial(len(c1_coeffs), c1_coeffs)

    return c0, c1


def mult_coeff(ct1, ct2, params):
    n, q, t, error_bound, delta = params

    c0_1, c1_1 = ct1
    c0_2, c1_2 = ct2

    c0 = c0_1.multiply_fft(c0_2)
    c0 = c0.scalar_multiply(1 / delta)
    c0 = c0.round().mod(q)

    c1 = c0_1.multiply_fft(c1_2).add(c1_1.multiply_fft(c0_2))
    c1 = c1.scalar_multiply(1 / delta)
    c1 = c1.round().mod(q)

    c2 = c1_1.multiply_fft(c1_2)
    c2 = c2.scalar_multiply(1 / delta)
    c2 = c2.round().mod(q)

    return c0, c1, c2


def relinearize_v1(c0, c1, c2, relin_key, q):
    num_levels = len(relin_key)

    base = ceil(sqrt(q))
    c2_decomposed = c2.base_decompose(base, num_levels)

    new_c0 = c0
    new_c1 = c1

    for i in range(num_levels):
        new_c0 = new_c0.add(relin_key[i][0].multiply(c2_decomposed[i], q), q)
        new_c1 = new_c1.add(relin_key[i][1].multiply(c2_decomposed[i], q), q)

    return new_c0, new_c1


def relinearize_v2(c0, c1, c2, relin_key, q, p):
    rlk_0, rlk_1 = relin_key

    c2_k0 = c2.multiply_fft(rlk_0)
    c2_k1 = c2.multiply_fft(rlk_1)

    # c2_0 = c2_k0.divide_round(p)
    c2_0 = c2_k0.scalar_multiply(1 / p)
    c2_0 = c2_0.mod(q)

    # c2_1 = c2_k1.divide_round(p)
    c2_1 = c2_k1.scalar_multiply(1 / p)
    c2_1 = c2_1.round().mod(q)

    new_c0 = c0.add(c2_0).mod(q)
    new_c1 = c1.add(c2_1).mod(q)

    return new_c0, new_c1


def eval_mult_v1(ct1, ct2, relin_key_v1, params):
    n, q, t, error_bound, delta = params
    c0, c1, c2 = mult_coeff(ct1, ct2, params)
    return relinearize_v1(c0, c1, c2, relin_key_v1, q)


def eval_mult_v2(ct1, ct2, relin_key_v2, params, p):
    n, q, t, error_bound, delta = params
    c0, c1, c2 = mult_coeff(ct1, ct2, params)
    return relinearize_v2(c0, c1, c2, relin_key_v2, q, p)
