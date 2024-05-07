import numpy as np

def summands(x_0=0.0, terms=1):
    summands = [x_0]
    for n in range(1, terms+1):
        if n % 2 == 0 or n<2:
            continue
        factor = -(x_0*x_0)/(n*(n-1))
        summands.append(summands[-1]*factor)
    return summands

def sum(arr, reverse=False):
    sum = 0
    if reverse:
        arr = arr[::-1]
    for num in arr:
        sum += num
    return sum

def sum_v2(arr, reverse=False):
    if reverse:
        arr = arr[::-1]
    sum_positive = sum(arr[::2])
    sum_negative = sum(arr[1::2])
    return sum_positive + sum_negative

def sin_interval(x, terms=1):
    sinarr = summands(terms=terms, x_0=x)
    result = sum(sinarr, reverse=True)
    return result

def reduce_sin(x, terms=1):
    if x < 0:
        return -reduce_sin(-x, terms=terms)
    if x <= np.pi:
        return sin_interval(x, terms=terms)
    elif x <= 2 * np.pi:
        return -sin_interval(np.abs(np.pi - x), terms=terms)
    else:
        k = np.floor(x / (2 * np.pi))
        r = x - 2 * k * np.pi
        return reduce_sin(r, terms=terms)

def main():
    a = 3.1
    sinarr = summands(x_0=a, terms=15)
    print(sum(sinarr))
    print(reduce_sin(a, 15))
    print(np.sin(a))
    print((reduce_sin(a, 25)-np.sin(a))/np.sin(a))


if __name__ == "__main__":
    main()
