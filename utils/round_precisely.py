import gmpy2


def precise_round_with_mod(coef_list, p, mod):
    rounded_list = []
    half_p = gmpy2.mpz(p) // 2

    for coef in coef_list:
        # Exact integer rounding with half-addition to ensure consistency
        rounded_coef = gmpy2.divexact(coef + half_p, p)
        rounded_mod_coef = gmpy2.mpz(rounded_coef % mod)
        rounded_list.append(rounded_mod_coef)

    return rounded_list


def divide_and_round(poly, p):
    rounded_coeff = [(coeff + p // 2) // p for coeff in poly]
    return rounded_coeff
