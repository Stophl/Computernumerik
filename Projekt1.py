import numpy as np
import matplotlib.pyplot as plt


def summands_for_exp(terms=1, x_0=0.0):
    invert = False
    if x_0 < 0:
        invert = True
    x_0 = np.abs(x_0)
    summands = [1, x_0]

    for n in range(terms + 1):
        if n < 2:
            continue
        factor = x_0 / n
        summands.append(summands[-1] * factor)
    summands = np.array(summands)

    return invert, summands


def taylorExponential_allsums(terms=1, x_0=0.0):
    invert, summands = summands_for_exp(terms=terms, x_0=x_0)

    forward = 0
    backward = 0
    for n in range(len(summands)):
        forward = forward + summands[n]
        reverse = len(summands) - 1 - n
        backward = backward + summands[reverse]
    check_sum = np.sum(summands)

    if invert:
        forward = 1 / forward
        backward = 1 / backward
        check_sum = 1 / check_sum

    return check_sum, forward, backward


def taylorExponential_back(terms=1, x_0=0.0):
    invert, summands = summands_for_exp(terms=terms, x_0=x_0)

    backward = 0
    for n in range(len(summands)):
        reverse = len(summands) - 1 - n
        backward = backward + summands[reverse]

    if invert:
        backward = 1.0 / backward

    return backward


def taylor_exp(x):
    ln2 = np.log(2.0)
    if np.abs(x) > ln2:
        invert = False
        if x < 0:
            invert = True
        x = np.abs(x)

        u = int(x / ln2)
        r = x - u * ln2
        result = (np.power(2.0, u) * taylorExponential_back(terms=10, x_0=r))
        if invert:
            result = 1.0 / result
    else:
        result = taylorExponential_back(terms=10, x_0=x)
    return result


def allErrors(terms=1, x_0=0.0, output=False, long=False):
    np_sum, forward_sum, backward_sum = taylorExponential_allsums(terms=terms, x_0=x_0)
    if long:
        true_value = np.exp(np.longdouble(x_0))
    else:
        true_value = np.exp(x_0)

    taylor_value = taylor_exp(x_0)

    taylor_value_abs_error = np.abs(taylor_value - true_value)
    taylor_value_rel_error = taylor_value_abs_error / true_value

    forward_abs_error = np.abs(forward_sum - true_value)
    forward_rel_error = forward_abs_error / true_value

    backward_abs_error = np.abs(backward_sum - true_value)
    backward_rel_error = backward_abs_error / true_value

    check_abs_error = np.abs(np_sum - true_value)
    check_rel_error = check_abs_error / true_value

    if output:
        print('"True" Value:', f"{true_value:.3E}")
        print("Taylor Value:", f"{taylor_value:.3E}")
        print("Taylor Absolute Error:", f"{taylor_value_abs_error:.3E}")
        print("Taylor Relative Error:", f"{taylor_value_rel_error:.3E}")
        print("Forward Sum Value:", f"{forward_sum:.3E}")
        print("Forward Absolute Error:", f"{forward_abs_error:.3E}")
        print("Forward Relative Error:", f"{forward_rel_error:.3E}")
        print("Backward Sum Value:", f"{backward_sum:.3E}")
        print("Backward Absolute Error:", f"{backward_abs_error:.3E}")
        print("Backward Relative Error:", f"{backward_rel_error:.3E}")
        print("CheckSum Value:", f"{np_sum:.3E}")
        print("CheckSum Absolute Error:", f"{check_abs_error:.3E}")
        print("CheckSum Relative Error:", f"{check_rel_error:.3E}")

    return (true_value,
            taylor_value, taylor_value_abs_error, taylor_value_rel_error,
            forward_sum, forward_abs_error, forward_rel_error,
            backward_sum, backward_abs_error, backward_rel_error, np_sum,
            check_abs_error, check_rel_error)


def txt_write_error(term_values=None, x_0_values=None, long=False, name="error_analysis_table.txt"):
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
                true_value, taylor_value, taylor_value_abs_error, taylor_value_rel_error, \
                    forward_sum, forward_abs_error, forward_rel_error, \
                    backward_sum, backward_abs_error, backward_rel_error, \
                    check_sum, check_abs_error, check_rel_error = allErrors(terms=terms, x_0=x_0, long=long)
                if first_term:
                    txtfile.write(f"True Value = {true_value:.3e}\n")

                    txtfile.write(
                        "{:<7} {:<10} {:<10} {:<12} {:<10} {:<10} {:<12} {:<10} {:<10} {:<12} "
                        "{:<10} {:<10} {:10}\n".format
                        ("terms",
                         "Taylor",
                         "Abs Err",
                         "Rel Err",
                         "Forward",
                         "Abs Err",
                         "Rel Err",
                         "Backward",
                         "Abs Err",
                         "Rel Err",
                         "CheckSum",
                         "Abs Err",
                         "Rel Err"))
                    first_term = False

                txtfile.write(
                    "{:<7} {:<10.3e} {:<10.3e} {:<12.3e} {:<10.3e} {:<10.3e} {:<12.3e} {:<10.3e} {:<10.3e} "
                    "{:<12.3e} {:<10.3e} {:<10.3e} {:<10.3e}\n".format(
                        terms,
                        taylor_value, taylor_value_abs_error, taylor_value_rel_error,
                        forward_sum, forward_abs_error, forward_rel_error,
                        backward_sum, backward_abs_error, backward_rel_error,
                        check_sum, check_abs_error, check_rel_error))

    print("Data has been written to", txt_file_path)


def plotting_delta(x_values, function, *args, title="Plot of δ(x) - default", \
                   logscale=False, helping_line=True, long=False):
    def delta_long(x):
        return np.abs((np.exp(np.longdouble(x)) - function(*args, x)) / np.exp(np.longdouble(x)))

    def delta(x):
        return np.abs(np.exp(x) - function(*args, x)) / np.exp(x)

    save = title.split("-", 1)[1].strip()


    y_values = []
    for item in x_values:
        y_values.append(delta(item))
    plt.plot(x_values, y_values, label='Standard Precision')

    if long:
        y_values_long = []
        for item in x_values:
            y_values_long.append(delta_long(item))
        plt.plot(x_values, y_values_long, label='Extended Precision')
        save = save +"_long"

    if logscale:
        plt.yscale('log')
        save = save + "_log"

    plt.xlabel('x')
    plt.ylabel('δ(x)')
    plt.title(title)
    plt.grid(True)

    lower = min(x_values)
    upper = max(x_values)
    save = save + f"_{round(lower, 2)}_{round(upper, 2)}"
    if helping_line:

        multiples_of_log2 = np.arange(np.ceil(lower / np.log(2)), np.floor(upper / np.log(2)) + 1)
        for multiple in multiples_of_log2:
            plt.axvline(x=multiple * np.log(2), color='purple', linestyle='--', alpha=0.5)
        plt.plot([], [], color='purple', linestyle='--', label="Multiples of ln(2)")
        save = save + "_helping"

    plt.legend()

    plt.savefig(f"{save}.pdf", format="pdf", dpi=300)
    plt.show()


def main():
    error = True
    if error:
        txt_write_error(term_values=[10, 50], x_0_values=[1.0, 10.0, 20.0, -10.0, -20.0], long=False,
                        name="error_analysis_table_normal.txt")
        txt_write_error(term_values=[10, 50], x_0_values=[1.0, 10.0, 20.0, -10.0, -20.0], long=True,
                        name="error_analysis_table_long.txt")

    plot = True
    if plot:
        x_values = np.linspace(0, np.log(2), 100)
        plotting_delta(x_values, taylorExponential_back, 10, title="Plot of δ(x) - Backward Summation",
                       logscale=False, helping_line=False)
        plotting_delta(x_values, taylor_exp, title="Plot of δ(x) - Taylor Algorithm",
                       logscale=False, helping_line=False)

        x_values = np.linspace(0, 5, 1000)
        plotting_delta(x_values, taylorExponential_back, 10, title="Plot of δ(x) - Backward Summation",
                       logscale=False, helping_line=True)
        plotting_delta(x_values, taylor_exp, title="Plot of δ(x) - Taylor Algorithm",
                       logscale=False, helping_line=True)

        x_values = np.linspace(0, 1, 300)
        plotting_delta(x_values, taylorExponential_back, 10, title="Plot of δ(x) - Backward Summation",
                       logscale=True, helping_line=False, long=True)
        plotting_delta(x_values, taylor_exp, title="Plot of δ(x) - Taylor Algorithm",
                       logscale=True, helping_line=False, long=True)


if __name__ == "__main__":
    main()
