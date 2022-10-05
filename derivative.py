import math

def derivative(f, a, is_good):
    dx = 0.1
    prev = 0
    xn = a
    cnt = 0
    err = 1
    while not is_good(prev, xn, err, cnt):
        lim = (f(a + dx) - f(a))/dx
        err = abs(lim - prev)
        prev = lim
        dx *= 0.1
        cnt += 1

    return prev

def is_good(prev, xn, err, cnt):
    return err <= 10**(-15) or cnt > 500


print(derivative(lambda x: math.pow(math.e, x), 2, is_good))
print(derivative(lambda x: x**2, 2, is_good))
