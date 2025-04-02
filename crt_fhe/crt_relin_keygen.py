from math import ceil, sqrt, floor, log
import gmpy2
from utils.random_samples import gaussian_distribution, uniform_distribution
from utils.poly_operations import poly_add, poly_mult
from utils.crt_function import crt_encode_e
from util.polynomial import Polynomial


def relin_keygen_v1(secret_key, pp):
    base = ceil(sqrt(pp[1]))
    num_levels = floor(log(pp[1], base)) + 1

    relin_keys = []
    power = 1
    secret_key = secret_key.coeffs

    sk_squared = poly_mult(secret_key, secret_key, pp[0], pp[1])

    for i in range(num_levels):
        ai = uniform_distribution(0, pp[1], pp[0])
        ei = gaussian_distribution(pp[0], 0, 1, pp[4])
        crt_e = crt_encode_e(pp, ei)
        t_i = base ** i

        sk_times_ai = poly_mult(secret_key, ai, pp[0])
        t_sk2 = [t_i * coeffs for coeffs in sk_squared]

        k0_temp = poly_add(sk_times_ai, crt_e, pp[1])
        k0 = poly_add(k0_temp, t_sk2, pp[1])

        k1 = [(-1) * coeff for coeff in ai]

        k0 = Polynomial(pp[0], k0)
        k1 = Polynomial(pp[0], k1)

        relin_keys.append((k0, k1))
        power = (power * base) % pp[1]

    return relin_keys


def relin_keygen_v2(secret_key, pp):
    n, q, error_bound, p1, p2 = pp
    p = q ** 3 + 1
    new_modulus = (q * p)
    secret_key = secret_key.coeffs

    ai = uniform_distribution(0, new_modulus, n)
    ei = gaussian_distribution(pp[0], 0, 1, p2)
    crt_e = crt_encode_e(pp, ei)

    sk_squared = poly_mult(secret_key, secret_key, n)
    secret_part = [p * coeffs for coeffs in sk_squared]

    a_s = poly_mult(ai, secret_key, n)
    as_e = poly_add(a_s, crt_e)
    b = poly_add(as_e, secret_part, new_modulus)

    ek0 = b
    ek1 = [(-1) * coeffs for coeffs in ai]

    ek0 = Polynomial(pp[0], ek0)
    ek1 = Polynomial(pp[0], ek1)

    relin_keys = ek0, ek1

    return relin_keys, p

