import time

from crt_params import public_params
from crt_keygen import keygen
from crt_relin_keygen import relin_keygen_v1, relin_keygen_v2
from crt_encryptor import encrypt
from crt_decryptor import decrypt
from utils.poly_operations import pad_with_zeros, poly_mult, poly_add
from crt_evaluate import eval_add, eval_mult_v1, eval_mult_v2
from utils.poly_operations import correctness

# Define security parameter
lam = 20

# Get BFV parameters
params = public_params(lam)
n, q, error_bound, p1, p2 = params
print("-"*50)
print("Parameters are: ")
print(f"n = {n}, q = {q}, p1 = {p1}, p2 = {p2}")
print("-"*50)

# Generate keys
public_key, secret_key = keygen(params)

# Print the keys
# print("Public Key:")
# print("b:", public_key[0].coeffs)
# print("a:", public_key[1].coeffs)
#
# print("\nSecret Key:")
# print("s:", secret_key.coeffs)

# # Generate relinearisation keys
relin_key_v1 = relin_keygen_v1(secret_key, params)
relin_key_v2, p = relin_keygen_v2(secret_key, params)

# print("Relinearisation Keys V1:")
# for idx, (k0, k1) in enumerate(relin_key_v1):
#     print(f"\nKey Level {idx + 1}:")
#     print(f"k0 coefficients: {k0.coeffs}")
#     print(f"k1 coefficients: {k1.coeffs}")
#
# print("Relinearisation Key V2:")
# print("rlk[0]: ", relin_key_v2[0].coeffs)
# print("rlk[1]: ", relin_key_v2[1].coeffs)

m = [1, 2, 3, 4]

# m = list(map(int, input("Enter message polynomial: ").split()))
m = pad_with_zeros(m, n)

c = encrypt(m, params, public_key)
d = decrypt(c, secret_key, params)

# c0_poly, c1_poly = encrypt(m, params, public_key)

# print("Enc(m):")
# print("c0:", c[0].coeffs)
# print("c1:", c[1].coeffs)
#
# print(f"Dec(c): {d}")

if correctness(m, d):
    print("CRT FHE Enc-Dec successful!")

print("-"*50)
print("Executing Homomorphic evaluations using FFT")
print("-"*50)
c_sum = eval_add(c, c, q)
# print("c+c =")
# print("c0_sum:", c_sum[0].coeffs)
# print("c1_sum:", c_sum[1].coeffs)

d_sum = decrypt(c_sum, secret_key, params)
# print("Dec (c + c) = ", d_sum)
if correctness(d_sum, (poly_add(m, m, p1))):
    print("CRT_FHE.add executed successfully!")

c_mult_v1 = eval_mult_v1(c, c, relin_key_v1, params)
# c_mult_v2 = eval_mult_v2(c, c, relin_key_v2, params, p)
#
d_mult1 = decrypt(c_mult_v1, secret_key, params)
# print("Dec (c * c) V1 = ", d_mult1)
#
# d_mult2 = decrypt(c_mult_v2, secret_key, params)
# # print("Dec (c * c) V2 = ", d_mult2)
#
m12 = poly_mult(m, m, n, p1)
# print(f"(m * m) = {m12}")
# #
if correctness(m12, d_mult1):
    print("CRT_FHE.mult_v1 executed successfully!")

# if correctness(m12, d_mult2):
#     print("CRT_FHE.mult_v2 executed successfully!")

print("-"*50)

