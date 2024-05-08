import numpy as np
import matplotlib.pyplot as plt
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

def rel_error(test=0.0, function=reduce_sin, terms=1):
    true_value = np.longdouble(np.sin(test))
    test_value = function(test, terms)
    if true_value != 0:
        return (np.abs((test_value-true_value)/true_value))
    else:
        return 0

def abs_error(test=0.0, function=reduce_sin, terms=1):
    true_value = np.longdouble(np.sin(test))
    test_value = function(test, terms)
    return np.abs((test_value-true_value))

def terms_needed_forRelErr(test=0.0, function=reduce_sin, terms=1, accuracy = 1e-13, limit=100):
    while rel_error(test=test, function=function, terms=terms) > accuracy and terms<limit:
        terms += 1
    return terms

def terms_needed_forAbsErr(test=0.0, function=reduce_sin, terms=1, accuracy = 1e-13, limit=100):
    while abs_error(test=test, function=function, terms=terms) > accuracy and terms<limit:
        terms += 1
    return terms

def plot_terms_needed_rel(x_arr=None, function=reduce_sin, accuracy=1e-13):
    if x_arr is None:
        x_arr = [1.0]
    rel_arr=[]
    for x in x_arr:
        rel_arr.append(terms_needed_forRelErr(test=x, function=function, accuracy=accuracy))
    plt.plot(x_arr, rel_arr)
    plt.show()
    return

def plot_terms_needed_abs(x_arr=None, function=reduce_sin, accuracy=1e-13):
    if x_arr is None:
        x_arr = [1.0]
    abs_arr=[]
    for x in x_arr:
        abs_arr.append(terms_needed_forAbsErr(test=x, function=function, accuracy=accuracy))
    plt.plot(x_arr, abs_arr)
    plt.show()
    return

def plot_rel_error_interval(start, end, num_points, functions=None, terms=1):
    if functions is None:
        functions = [reduce_sin]
    x_arr = np.linspace(start, end, num_points)
    for func in functions:
        rel_arr = []
        for x in x_arr:
            rel_arr.append(rel_error(test=x, function=func, terms=terms))
        plt.plot(x_arr, rel_arr, label=func.__name__)
    plt.xlabel('x')
    plt.ylabel('Relative Error')
    plt.title('Relative Error for Interval [{}, {}]'.format(start, end))
    plt.yscale('log')
    plt.legend()
    plt.show()

def main():
    a = 0.0
    print(np.sin(a))
    print(sin_interval(a, 15))
    print(rel_error(test=a, function=sin_interval, terms=15))
    print(abs_error(test=a, function=sin_interval, terms=15))
    print(terms_needed_forRelErr(a, sin_interval))
    print(terms_needed_forAbsErr(a, sin_interval, accuracy=1e-15))
    print(reduce_sin(a, 15))
    print(rel_error(test=a, function=reduce_sin, terms=15))
    print(abs_error(test=a, function=reduce_sin, terms=15))
    print(terms_needed_forRelErr(a, reduce_sin))
    print(terms_needed_forAbsErr(a, reduce_sin, accuracy=1e-15))
    plot_rel_error_interval(-5, 5, 1001, functions=[reduce_sin, sin_interval], terms=10)





if __name__ == "__main__":
    main()