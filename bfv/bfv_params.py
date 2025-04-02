import math
import gmpy2
import random

def gen_modulus(n, lam):
    while True:
        q = random.randint(2 ** (lam - 1), 2 ** lam - 1)
        q = gmpy2.next_prime(q)
        if (q - 1) % n == 0:
            return int(q)

def bfv_params(lam):
    n = int(input("Enter degree of polynomial (in powers of 2): "))
    # n = 4
    q = gen_modulus(n, lam)  # Generate modulus q using the new function
    t = 17

    # Standard deviations for error bounds
    std1 = 1
    std2 = 2
    error_bound = (std1, std2)

    # Scaling factor
    delta = math.floor(q // t)

    params = (n, q, t, error_bound, delta)

    return params
