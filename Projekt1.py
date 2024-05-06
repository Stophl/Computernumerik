import numpy as np
import matplotlib.pyplot as plt


def sums_for_exp(terms=1, x_0=0.0):
    summands = [1, x_0]

    for n in range(terms + 1):
        if n < 2:
            continue
        factor = x_0 / n
        summands.append(summands[-1] * factor)

    forward = 0
    backward = 0
    for n in range(len(summands)):
        forward = forward + summands[n]
        reverse = len(summands) - 1 - n
        backward = backward + summands[reverse]

    return forward, backward


def taylor_exp(terms=10, x_0=0.0):
    ln2 = np.log(2.0)

    invert = False
    if x_0 < 0:
        invert = True
    x = np.abs(x_0)

    if np.abs(x) > ln2:
        u = int(x / ln2)
        r = x - u * ln2
        result = (np.power(2.0, u) * sums_for_exp(terms=terms, x_0=r)[1])
    else:
        result = sums_for_exp(terms=terms, x_0=x)[1]

    if invert:
        result = 1.0 / result

    return result


def allErrors(terms=1, x_0=0.0, output=False):
    forward_sum, backward_sum = sums_for_exp(terms=terms, x_0=x_0)

    true_value = np.exp(np.longdouble(x_0))

    taylor_value = taylor_exp(terms=terms, x_0=x_0)

    taylor_value_abs_error = np.abs(taylor_value - true_value)
    taylor_value_rel_error = taylor_value_abs_error / true_value

    forward_abs_error = np.abs(forward_sum - true_value)
    forward_rel_error = forward_abs_error / true_value

    backward_abs_error = np.abs(backward_sum - true_value)
    backward_rel_error = backward_abs_error / true_value

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

    return (true_value,
            taylor_value, taylor_value_abs_error, taylor_value_rel_error,
            forward_sum, forward_abs_error, forward_rel_error,
            backward_sum, backward_abs_error, backward_rel_error)


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
                    backward_sum, backward_abs_error, backward_rel_error = allErrors(terms=terms, x_0=x_0)
                if first_term:
                    txtfile.write(f"True Value = {true_value:.3e}\n")

                    txtfile.write(
                        "{:<7} {:<10} {:<10} {:<12} {:<10} {:<10} {:<12} {:<10} {:<10} {:<12} \n".format
                        ("terms",
                         "Taylor",
                         "Abs Err",
                         "Rel Err",
                         "Forward",
                         "Abs Err",
                         "Rel Err",
                         "Backward",
                         "Abs Err",
                         "Rel Err"))
                    first_term = False

                txtfile.write(
                    "{:<7} {:<10.3e} {:<10.3e} {:<12.3e} {:<10.3e} {:<10.3e} {:<12.3e} {:<10.3e} {:<10.3e} "
                    "{:<12.3e}\n".format(
                        terms,
                        taylor_value, taylor_value_abs_error, taylor_value_rel_error,
                        forward_sum, forward_abs_error, forward_rel_error,
                        backward_sum, backward_abs_error, backward_rel_error))

    print("Data has been written to", txt_file_path)


def plotting_delta(x_values, function, terms=10, index=None, title="Plot of δ(x) - default", \
                   logscale=False, helping_line=True, fontsize=20):
    def delta(x):
        if index:
            return np.abs((np.exp(np.longdouble(x)) - function(terms, x)[index]) / np.exp(np.longdouble(x)))
        else:
            return np.abs((np.exp(np.longdouble(x)) - function(terms, x)) / np.exp(np.longdouble(x)))

    save = title.split("-", 1)[1].strip() + "_terms" + str(terms)

    y_values = []
    for item in x_values:
        y_values.append(delta(item))
    plt.plot(x_values, y_values, label="relativer Fehler")

    if logscale:
        plt.yscale('log')
        save = save + "_log"

    plt.xlabel('x', fontsize=fontsize)
    plt.ylabel('δ(x)', fontsize=fontsize)
    #plt.title(title, fontsize=fontsize)
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

    plt.legend(fontsize=fontsize - 5)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.gca().yaxis.get_offset_text().set_fontsize(fontsize - 3)
    plt.gca().ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.tight_layout()
    plt.savefig(f"{save}.pdf", format="pdf", dpi=300)
    plt.show()


def main():
    error = False
    if error:
        txt_write_error(term_values=[10, 50], x_0_values=[1.0, 10.0, 20.0, -10.0, -20.0],
                        name="error_analysis_table.txt")

    plot = True
    if plot:
        for item in [10, 11, 12, 13, 14, 15, 20]:
            x_values = np.linspace(0, 3, 1000)
            plotting_delta(x_values, sums_for_exp, item, index=1, title="Plot of δ(x) - Backward Summation",
                           logscale=False, helping_line=True)
            plotting_delta(x_values, taylor_exp, item, title="Plot of δ(x) - Taylor Algorithm",
                           logscale=False, helping_line=True)


if __name__ == "__main__":
    main()
