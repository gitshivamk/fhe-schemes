from math import ceil, sqrt, floor, log
from utils.random_samples import gaussian_distribution, uniform_distribution
from utils.poly_operations import poly_add, poly_mult
from util.polynomial import Polynomial


def relin_keygen_v1(secret_key, params):
    n, q, t, error_bound, delta = params

    base = ceil(sqrt(q))
    num_levels = floor(log(q, base)) + 1

    relin_keys = []
    power = 1

    sk_squared = secret_key.multiply(secret_key, q)

    for i in range(num_levels):

        ai_coeffs = uniform_distribution(0, q, n)
        ai = Polynomial(n, ai_coeffs)

        ei_coeffs = gaussian_distribution(n, 0, error_bound[0])
        ei = Polynomial(n, ei_coeffs)

        t_i = base ** i
        t_sk2 = sk_squared.scalar_multiply(t_i, q)

        sk_times_ai = secret_key.multiply(ai, q)
        k0_temp = sk_times_ai.add(ei, coeff_modulus=q).scalar_multiply(-1, q)
        k0 = k0_temp.add(t_sk2, coeff_modulus=q)

        k1 = ai
        relin_keys.append((k0, k1))

        power = (power * base) % q

    return relin_keys

def relin_keygen_v2(secret_key, params):
    n, q, t, error_bound, delta = params
    p = q ** 3 + 1
    new_modulus = (q * p)

    a = uniform_distribution(0, new_modulus, n)
    e = gaussian_distribution(n, 0, error_bound[0])

    sk_coeffs = secret_key.coeffs
    sk_squared = poly_mult(sk_coeffs, sk_coeffs, n)
    secret_part = [p * coeffs for coeffs in sk_squared]

    a_s = poly_mult(a, sk_coeffs, n)
    as_e = poly_add(a_s, e)
    as_e = [-1 * coeff for coeff in as_e]
    b = poly_add(as_e, secret_part, new_modulus)

    # Convert to polynomial
    b = Polynomial(n, b)
    a = Polynomial(n, a)

    ek0 = b
    ek1 = a
    relin_keys = ek0, ek1

    return relin_keys, p




