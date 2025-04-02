from util.ntt import NTTContext


class EncoderDecoder:
    def __init__(self, degree, modulus):
        self.ntt_context = NTTContext(degree, modulus)
        self.modulus = modulus

    def encode(self, plaintext):
        encoded_poly = self.ntt_context.ftt_inv(plaintext)
        return encoded_poly

    def decode(self, encoded_poly):
        decoded_values = self.ntt_context.ftt_fwd(encoded_poly)
        decoded_values = [val % self.modulus for val in decoded_values]
        return decoded_values
