import random


def gaussian_distribution(degree, mu, sigma, p2=None):
    if p2 is not None:
        e = [int(random.gauss(mu, sigma)) % p2 for _ in range(degree)]
    else:
        e = [int(random.gauss(mu, sigma)) for _ in range(degree)]
    return e


def ternary_distribution(degree):
    s = random.choices([-1, 0, 1], weights=[0.25, 0.5, 0.25], k=degree)

    if all(value == 0 for value in s):
        random_index = random.randint(0, degree - 1)
        s[random_index] = random.choice([-1, 1])

    return s


def uniform_distribution(min_value, max_value, degree):
    if degree == 1:
        return random.randrange(min_value, max_value)
    else:
        return [random.randrange(min_value, max_value) for _ in range(degree)]
