from math import sqrt


def mod_exp(base, exp, mod):
    result = 1
    while exp > 0:
        if exp % 2 == 1:
            result = (result * base) % mod
        base = (base * base) % mod
        exp //= 2
    return result


def find_primitive_root(n, q):
    factors = set()
    phi = q - 1
    temp = phi
    for i in range(2, int(sqrt(temp)) + 1):
        while temp % i == 0:
            factors.add(i)
            temp //= i
    if temp > 1:
        factors.add(temp)

    for g in range(2, q):
        if all(mod_exp(g, phi // factor, q) != 1 for factor in factors):
            return g
    return None


def ntt(a, q, omega):
    n = len(a)
    A = [0] * n
    for k in range(n):
        A[k] = sum(a[j] * mod_exp(omega, j * k, q) for j in range(n)) % q
    return A


def inverse_ntt(A, q, omega):
    n = len(A)
    omega_inv = mod_exp(omega, q - 2, q)
    a = [0] * n
    for k in range(n):
        a[k] = sum(A[j] * mod_exp(omega_inv, j * k, q) for j in range(n)) % q
        a[k] = (a[k] * mod_exp(n, q - 2, q)) % q
    return a
