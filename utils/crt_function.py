import gmpy2


def crt_encode_e(pp, e):
    p1, p2 = pp[3], pp[4]
    delta2 = int(p1 * gmpy2.invert(p1, p2))
    crt_e = [(delta2 * coeff) % (p1 * p2) for coeff in e]
    return crt_e


def crt_encode_m(pp, m):
    p1, p2 = pp[3], pp[4]
    delta1 = int(p2 * gmpy2.invert(p2, p1))
    crt_m = [(delta1 * coeff) % (p1 * p2) for coeff in m]
    return crt_m
