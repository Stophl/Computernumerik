import numpy as np
import matplotlib.pyplot as plt


def summands(x_0=0.0, terms=1):
    summands = [x_0]
    for n in range(1, terms + 1):
        if n % 2 == 0 or n < 2:
            continue
        factor = -(x_0 * x_0) / (n * (n - 1))
        summands.append(summands[-1] * factor)
    return summands


def summing(arr, reverse=False):
    sum = 0
    if reverse:
        arr = arr[::-1]
    for num in arr:
        sum += num
    return sum


def summing_v2(arr, reverse=False):
    if reverse:
        arr = arr[::-1]
    sum_positive = summing(arr[::2])
    sum_negative = summing(arr[1::2])
    return sum_positive + sum_negative


def sin_standard_taylor(x, terms=1, reverse=True):
    sinarr = summands(terms=terms, x_0=x)
    result = summing(sinarr, reverse=reverse)
    return result


def sin_reduction_taylor(x, terms=1):
    if x < 0:
        return -sin_reduction_taylor(-x, terms=terms)
    if x <= np.pi:
        return sin_standard_taylor(x, terms=terms)
    elif x <= 2 * np.pi:
        return -sin_standard_taylor(np.abs(np.pi - x), terms=terms)
    else:
        k = np.floor(x / (2 * np.pi))
        r = x - 2 * k * np.pi
        return sin_reduction_taylor(r, terms=terms)


def error(*reverse, test=0.0, function=sin_reduction_taylor, terms=1, absolute=False):
    true_value = np.longdouble(np.sin(test))
    test_value = function(test, terms, *reverse)
    if absolute:
        return np.abs((test_value - true_value))
    else:
        if true_value != 0:
            return np.abs((test_value - true_value) / true_value)
        else:
            return 0


def terms_needed(test=0.0, function=sin_reduction_taylor, terms=1, accuracy=1e-13, limit=100, absolute=False):
    while error(test=test, function=function, terms=terms, absolute=absolute) > accuracy and terms < limit:
        terms += 1
    return terms


def plot_terms_interval(start, end, num_points, functions=None, accuracy=1e-13, limit=500, absolute=False):
    if functions is None:
        functions = [sin_reduction_taylor]
    x_arr = np.linspace(start, end, num_points)
    for func in functions:
        term_arr = []
        for x in x_arr:
            term_arr.append(terms_needed(test=x, function=func, accuracy=accuracy, limit=limit, absolute=absolute))
        plt.plot(x_arr, term_arr, label=func.__name__)
    plt.xlabel('x')
    if absolute:
        err_type = "absolute"
    else:
        err_type = "relative"
    plt.ylabel('Terms')
    plt.title(f'Terms needed to have {accuracy} accuracy in {err_type} error \n Interval [{start}, {end}]')
    plt.legend()
    plt.tight_layout()
    plt.show()
    return


def plot_error_interval(start, end, num_points, functions=None, terms=1, absolute=False):
    if functions is None:
        functions = [sin_reduction_taylor]
    x_arr = np.linspace(start, end, num_points)
    for func in functions:
        error_arr = []
        for x in x_arr:
            error_arr.append(error(test=x, function=func, terms=terms, absolute=absolute))
        plt.plot(x_arr, error_arr, label=func.__name__)
    plt.xlabel('x')
    if absolute:
        err_type = "absolute"
    else:
        err_type = "relative"

    plt.ylabel(f'{err_type} Error')
    plt.title(f'{err_type} Error for Interval [{start}, {end}]')
    # plt.yscale('log')
    plt.legend()
    plt.tight_layout()
    plt.show()
    return


def txt_write_error(term_values=None, x_0_values=None):
    if term_values is None:
        term_values = [1]
    if x_0_values is None:
        x_0_values = [0]
    txt_file_path = 'error_analysis_table.txt'

    with open(txt_file_path, 'w') as txtfile:
        txtfile.write(
            "{:<8} {:<8} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}\n".format("terms", "x_0",
                                                                                                  "True Value",
                                                                                                  "reduction (x)",
                                                                                                  "Abs Err",
                                                                                                  "Rel Err",
                                                                                                  "Forward Sum",
                                                                                                  "Abs Err",
                                                                                                  "Rel Err",
                                                                                                  "Backward Sum",
                                                                                                  "Abs Err",
                                                                                                  "Rel Err",
                                                                                                  ))

        for x_0 in x_0_values:
            for terms in term_values:
                true_value = np.longdouble(np.sin(x_0))
                reduction = sin_reduction_taylor(x_0, terms)
                reduction_rel = error(test=x_0, function=sin_reduction_taylor, terms=terms)
                reduction_abs = error(test=x_0, function=sin_reduction_taylor, terms=terms, absolute=True)
                forward_sum = sin_standard_taylor(x_0, terms, reverse=False)
                forward_rel_error = error(False, test=x_0, function=sin_standard_taylor, terms=terms, absolute=False)
                forward_abs_error = error(False, test=x_0, function=sin_standard_taylor, terms=terms, absolute=True)
                backward_sum = sin_standard_taylor(x_0, terms, reverse=False)
                backward_rel_error = error(test=x_0, function=sin_standard_taylor, terms=terms, absolute=False)
                backward_abs_error = error(test=x_0, function=sin_standard_taylor, terms=terms, absolute=True)

                txtfile.write(
                    "{:<8} {:<8} {:<15.3e} {:<15.3e} {:<15.3e} {:<15.3e} {:<15.3e} {:<15.3e} {:<15.3e} {:<15.3e} {:<15.3e} {:<15.3e}\n".format(
                        terms, x_0, true_value, reduction, reduction_abs, reduction_abs,
                        forward_sum, forward_abs_error, forward_rel_error, backward_sum,
                        backward_abs_error, backward_rel_error))

    print("Data has been written to", txt_file_path)
    return


def main():
    plot_terms_interval(0, 7, 701, functions=[sin_standard_taylor, sin_reduction_taylor],
                        accuracy=1e-13, limit=100)
    plot_terms_interval(0, 7, 701, functions=[sin_standard_taylor, sin_reduction_taylor],
                        accuracy=1e-15, limit=100)
    plot_terms_interval(0, 7, 701, functions=[sin_standard_taylor, sin_reduction_taylor],
                        accuracy=1e-13, limit=100, absolute=True)
    plot_terms_interval(0, 7, 701, functions=[sin_standard_taylor, sin_reduction_taylor],
                        accuracy=1e-15, limit=100, absolute=True)
    plot_error_interval(0, 7, 701, functions=[sin_standard_taylor, sin_reduction_taylor],
                        terms=25)
    plot_error_interval(0, 7, 701, functions=[sin_standard_taylor, sin_reduction_taylor],
                        terms=25, absolute=True)
    txt_write_error([10,15,20,25,30], [2,4,8,16,32,64,128])


if __name__ == "__main__":
    main()
