
def coeffs_mod(y, modulus):
    return [(coeff % modulus + modulus) % modulus for coeff in y]

def poly_mult(coeffs1, coeffs2, ring_degree, coeff_modulus=None):
    poly_prod = [0] * ring_degree

    for d in range(2 * ring_degree - 1):
        index = d % ring_degree
        sign = int(d < ring_degree) * 2 - 1

        coeff = 0
        for i in range(ring_degree):
            if 0 <= d - i < ring_degree:
                coeff += coeffs1[i] * coeffs2[d - i]
        poly_prod[index] += sign * coeff

        if coeff_modulus:
            poly_prod[index] %= coeff_modulus

    return poly_prod


def poly_add(x, y, coeff_modulus=None):
    poly_sum = [0] * len(x)

    for i in range(len(x)):
        poly_sum[i] = (x[i] + y[i]) % coeff_modulus if coeff_modulus else x[i] + y[i]
        if coeff_modulus:
            poly_sum[i] = (poly_sum[i] % coeff_modulus + coeff_modulus) % coeff_modulus

    return poly_sum


def base_decompose(poly, base, num_levels):
    decomposed = [[0] * len(poly) for _ in range(num_levels)]
    current_poly = poly[:]

    for i in range(num_levels):
        decomposed[i] = [coeff % base for coeff in current_poly]
        current_poly = [coeff // base for coeff in current_poly]

    return decomposed


def pad_with_zeros(lst, n):
    if len(lst) < n:
        lst.extend([0] * (n - len(lst)))
    return lst


def correctness(a, b):
    if a == b:
        return True
    else:
        print("--- It's an error! ---")
        print("Expected value:", a)
        print("Decrypted value: ", b)
        return False
