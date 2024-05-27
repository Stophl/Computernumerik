import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal, getcontext


def summands(x_0=0.0, terms=1):
    summands = [x_0]
    for n in range(1, 2 * terms + 1):
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


def sin_standard_taylor(x, terms=1, reverse=True):
    sinarr = summands(terms=terms, x_0=x)
    result = summing(sinarr, reverse=reverse)
    return result


def sin_reduction_taylor(x, terms=1, high_precision=False):
    if x < 0:
        return -sin_reduction_taylor(-x, terms=terms)
    if x <= np.pi:
        return sin_standard_taylor(x, terms=terms)
    elif x <= 2 * np.pi:
        return -sin_standard_taylor(np.abs(np.pi - x), terms=terms)
    else:
        if high_precision:
            getcontext().prec = 32
            x_decimal = Decimal(x)
            pi = Decimal(np.pi)
            k = Decimal(np.floor(x / (2 * np.pi)))
            r = x_decimal - Decimal(2 * k * pi)
            return sin_reduction_taylor(float(r), terms=terms)
        else:
            k = np.floor(x / (2 * np.pi))
            r = x - 2 * k * np.pi
            return sin_reduction_taylor(r, terms=terms)


def error(*reverse_high, test=0.0, function=sin_reduction_taylor, terms=1, absolute=False):
    true_value = np.longdouble(np.sin(test))
    test_value = function(test, terms, *reverse_high)
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


def txt_write_error(term_values=None, x_0_values=None, name="error_analysis_table.txt"):
    if term_values is None:
        term_values = [1]
    if x_0_values is None:
        x_0_values = [0]
    txt_file_path = name

    with (open(txt_file_path, 'w') as txtfile):
        for x_0 in x_0_values:
            txtfile.write("{:<15}".format(f"\nx_0 = {x_0}"))
            first_term = True
            for terms in term_values:
                true_value = np.longdouble(np.sin(x_0))
                reduction = sin_reduction_taylor(x_0, terms)
                reduction_rel = error(test=x_0, function=sin_reduction_taylor, terms=terms)
                reduction_abs = error(test=x_0, function=sin_reduction_taylor, terms=terms, absolute=True)
                reduction_high = sin_reduction_taylor(x_0, terms, high_precision=True)
                reduction_rel_high = error(True, test=x_0, function=sin_reduction_taylor, terms=terms)
                reduction_abs_high = error(True, test=x_0, function=sin_reduction_taylor, terms=terms, absolute=True)
                forward_sum = sin_standard_taylor(x_0, terms, reverse=False)
                forward_rel_error = error(False, test=x_0, function=sin_standard_taylor, terms=terms, absolute=False)
                forward_abs_error = error(False, test=x_0, function=sin_standard_taylor, terms=terms, absolute=True)
                backward_sum = sin_standard_taylor(x_0, terms)
                backward_rel_error = error(test=x_0, function=sin_standard_taylor, terms=terms, absolute=False)
                backward_abs_error = error(test=x_0, function=sin_standard_taylor, terms=terms, absolute=True)
                if first_term:
                    txtfile.write(f"True Value = {true_value:.3e}\n")

                    txtfile.write(
                        "{:<7} {:<10} {:<10} {:<12} {:<10} {:<10} {:<12} {:<10} {:<10} {:<12} {:<10} {:<10} {:<12} {:<12}\n".format
                        ("terms",
                         "Taylor",
                         "Abs Err",
                         "Rel Err",
                         "Tay-High",
                         "Abs Err",
                         "Rel Err",
                         "Forward",
                         "Abs Err",
                         "Rel Err",
                         "Backward",
                         "Abs Err",
                         "Rel Err",
                         "F-B"))
                    first_term = False

                txtfile.write(
                    "{:<7} {:<10.3e} {:<10.3e} {:<12.3e} {:<10.3e} {:<10.3e} {:<12.3e} {:<10.3e} {:<10.3e} {:<12.3e} {:<10.3e} {:<10.3e} "
                    "{:<12.3e} {:<12.3e}\n".format(
                        terms,
                        reduction, reduction_abs, reduction_rel,
                        reduction_high, reduction_abs_high, reduction_rel_high,
                        forward_sum, forward_abs_error, forward_rel_error,
                        backward_sum, backward_abs_error, backward_rel_error,
                        forward_abs_error - backward_abs_error))

    print("Data has been written to", txt_file_path)
    return


def neville(datax, datay, x):
    n = len(datax)
    p = n * [0]
    for k in range(n):
        for i in range(n - k):
            if k == 0:
                p[i] = datay[i]
            else:
                p[i] = ((x - datax[i + k]) * p[i] + (datax[i] - x) * p[i + 1]) / \
                       (datax[i] - datax[i + k])
    return p[0]


def generate_interpolation_data(interval_start, interval_end, num_points, chebyshev=False):
    if chebyshev:
        chebyshev_nodes = [(interval_start + interval_end) / 2 +
                           (interval_end - interval_start) / 2 *
                           np.cos((2 * i + 1) * np.pi / (2 * num_points))
                           for i in range(num_points)]
        datax = np.array(chebyshev_nodes)
        datay = np.sin(datax)
        return datax, datay
    else:
        datax = np.linspace(interval_start, interval_end, num_points)
        datay = np.sin(datax)
        return datax, datay


def error_neville(test=0, datax=None, datay=None, absolute=False):
    if datax is None:
        datax = [0]
    if datay is None:
        datay = [np.sin(0)]
    true_value = np.longdouble(np.sin(test))
    test_value = neville(datax, datay, test)

    if absolute:
        return np.abs((test_value - true_value))
    else:
        if true_value != 0:
            return np.abs((test_value - true_value) / true_value)
        else:
            return 0


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
    plt.title(f'Terms needed to have {accuracy:.2e} accuracy in {err_type} error \n Interval [{start}, {end}]')
    plt.legend()
    plt.tight_layout()
    plt.show()
    return


def plot_error_interval(start, end, num_points, functions=None, terms=1, absolute=False, log=False,
                        chebyshev=False, both_knots=False, stuetz=50):
    if functions is None:
        functions = [sin_reduction_taylor]

    x_arr = np.linspace(start, end, num_points)

    for func in functions:
        if func == neville:
            looper = [False, True] if both_knots else [chebyshev]

            for item in looper:
                error_arr = []
                for x in x_arr:
                    datax, datay = generate_interpolation_data(start, end, stuetz, chebyshev=item)
                    error_arr.append(error_neville(test=x, datax=datax, datay=datay, absolute=absolute))
                label = f'{func.__name__} ({"Chebyshev" if item else "Equidistant"})'
                plt.plot(x_arr, error_arr, label=label)
        else:
            error_arr = []
            for x in x_arr:
                error_arr.append(error(test=x, function=func, terms=terms, absolute=absolute))
            plt.plot(x_arr, error_arr, label=func.__name__)

    plt.xlabel('x')
    err_type = "absolute" if absolute else "relative"
    plt.ylabel(f'{err_type} Error')
    plt.title(f'{err_type} Error for Interval [{start}, {end}]')

    if log:
        plt.yscale('log')

    machine_precision = np.finfo(float).eps
    horizontal_line_value = 100 * machine_precision
    plt.axhline(y=horizontal_line_value, color='r', linestyle='--',
                label=f'100eps ({horizontal_line_value:.2e})')

    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_error_terms(xs, functions=None, max_terms=30, absolute=False):
    if functions is None:
        functions = [sin_reduction_taylor]

    if not isinstance(xs, (list, np.ndarray)):
        xs = [xs]

    terms = np.arange(1, max_terms + 1)

    for func in functions:
        for x in xs:
            errors = []
            for term in terms:
                err = error(test=x, function=func, terms=term, absolute=absolute)
                errors.append(err)
            plt.plot(terms, errors, label=f'{func.__name__} at x={x}')

    plt.xlabel('Number of Terms')
    err_type = "Absolute" if absolute else "Relative"
    plt.ylabel(f'{err_type} Error')
    plt.title(f'{err_type} Error vs. Number of Terms')
    plt.yscale('log')
    plt.axhline(y=100 * np.finfo(float).eps, color='r', linestyle='--',
                label=f'100eps ({100 * np.finfo(float).eps:.2e})')
    plt.legend()
    plt.tight_layout()
    plt.show()


def main():
    accuracy = 100 * np.finfo(np.double).eps
    terms = 13
    stuetz = 17
    plot_terms_interval(0, 7, 701, functions=[sin_standard_taylor, sin_reduction_taylor],
                        accuracy=accuracy, limit=100)
    plot_terms_interval(0, 7, 701, functions=[sin_standard_taylor, sin_reduction_taylor],
                        accuracy=accuracy, limit=100, absolute=True)
    plot_error_interval(0, 4, 201, functions=[sin_standard_taylor, sin_reduction_taylor],
                        terms=terms, log=True)
    plot_error_interval(0, 4, 201, functions=[sin_standard_taylor, sin_reduction_taylor],
                        terms=terms, absolute=True, log=True)
    plot_error_interval(0, 4, 201, functions=[sin_standard_taylor, sin_reduction_taylor, neville],
                        terms=terms, log=True, both_knots=True, stuetz=stuetz)
    plot_error_interval(0, 4, 201, functions=[sin_standard_taylor, sin_reduction_taylor, neville],
                        terms=terms, absolute=True, log=True, both_knots=True, stuetz=stuetz)
    plot_error_interval(707, 711, 401, functions=[sin_reduction_taylor, neville],
                        terms=terms, log=True, both_knots=True, stuetz=stuetz)
    plot_error_terms([2, 4, 8, 16, 32, 64, 128, 710],max_terms=50)
    plot_error_terms([2, 4, 8, 16, 32, 64, 128], functions=[sin_standard_taylor], max_terms=50)


    txt_write_error([2, 4, 6, 8, 10, 12, 13, 14, 15, 20, 25, 30], [0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 710])


if __name__ == "__main__":
    main()
