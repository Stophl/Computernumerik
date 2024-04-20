
import numpy as np
import csv
import matplotlib.pyplot as plt


def taylorExponential_allsums(terms=1, x_0=0):
    invert = False
    if x_0 < 0:
        invert = True
    x_0 = np.abs(x_0)
    summands = [1, x_0]

    for n in range(terms+1):
        if n < 2:
            continue
        factor = x_0/n
        summands.append(summands[-1]*factor)
    summands = np.array(summands)

    forward = 0
    backward = 0
    for n in range(len(summands)):
        forward = forward + summands[n]
        reverse = len(summands) - 1 - n
        backward = backward + summands[reverse]
    check_sum = np.sum(summands)

    if invert:
        forward = 1/forward
        backward = 1/backward
        check_sum = 1/check_sum

    return summands, check_sum, forward, backward

def allErrors(terms=1, x_0=0, output=False):
    summands, np_sum, forward_sum, backward_sum = taylorExponential_allsums(terms=terms, x_0=x_0)
    true_value = np.exp(x_0)

    forward_abs_error = np.abs(forward_sum - true_value)
    forward_rel_error = forward_abs_error / true_value

    backward_abs_error = np.abs(backward_sum - true_value)
    backward_rel_error = backward_abs_error / true_value

    check_abs_error = np.abs(np_sum - true_value)
    check_rel_error = check_abs_error / true_value

    if output:
        print('"True" Value:', true_value)
        print("Forward Sum Value:", forward_sum)
        print("Forward Absolute Error:", forward_abs_error)
        print("Forward Relative Error:", forward_rel_error)
        print("Backward Sum Value:", backward_sum)
        print("Backward Absolute Error:", backward_abs_error)
        print("Backward Relative Error:", backward_rel_error)
        print("CheckSum Absolute Error:", check_abs_error)
        print("CheckSum Relative Error:", check_rel_error)

    return true_value, forward_sum, forward_abs_error, forward_rel_error, backward_sum, backward_abs_error, backward_rel_error, check_abs_error, check_rel_error

def txt_write_error(term_values=[1], x_0_values=[0]):
    txt_file_path = 'error_analysis_table.txt'

    with open(txt_file_path, 'w') as txtfile:
        txtfile.write(
            "{:<8} {:<8} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}\n".format("terms", "x_0", "True Value",
                                                                                           "Forward Sum",
                                                                                           "Forward Abs Err",
                                                                                           "Forward Rel Err",
                                                                                           "Backward Sum",
                                                                                           "Backward Abs Err",
                                                                                           "Backward Rel Err",
                                                                                           "CheckSum Abs Err"))

        for x_0 in x_0_values:
            for terms in term_values:
                true_value, forward_sum, forward_abs_error, forward_rel_error, backward_sum, backward_abs_error, backward_rel_error, check_abs_error, check_rel_error = allErrors(
                    terms=terms, x_0=x_0)

                txtfile.write(
                    "{:<8} {:<8} {:<15.10e} {:<15.10e} {:<15.10e} {:<15.10e} {:<15.10e} {:<15.10e} {:<15.10e} {:<15.10e}\n".format(
                        terms, x_0, true_value, forward_sum, forward_abs_error, forward_rel_error, backward_sum,
                        backward_abs_error, backward_rel_error, check_abs_error))

    print("Data has been written to", txt_file_path)

def csv_write_error(term_values=[1], x_0_values=[0]):
    csv_file_path = 'error_analysis_table.csv'

    with open(csv_file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')

        writer.writerow(
            ["terms", "x_0", "True Value", "Forward Sum", "Forward Abs Err", "Forward Rel Err", "Backward Sum",
             "Backward Abs Err", "Backward Rel Err"])

        for x_0 in x_0_values:
            for terms in term_values:
                true_value, forward_sum, forward_abs_error, forward_rel_error, backward_sum, backward_abs_error, backward_rel_error, check_abs_error, check_rel_error = allErrors(
                    terms=terms, x_0=x_0)

                writer.writerow(
                    [terms, x_0, true_value, forward_sum, forward_abs_error, forward_rel_error, backward_sum,
                     backward_abs_error, backward_rel_error])

    print("Data has been written to", csv_file_path)

def taylorExponential_back(terms=1, x_0=0):
    invert = False
    if x_0 < 0:
        invert = True
    x_0 = np.abs(x_0)
    summands = [1, x_0]

    for n in range(terms+1):
        if n < 2:
            continue
        factor = x_0/n
        summands.append(summands[-1]*factor)

    backward = 0
    for n in range(len(summands)):
        reverse = len(summands) - 1 - n
        backward = backward + summands[reverse]


    if invert:
        backward = 1 / backward

    return backward

def plotting_delta(x_values, function, *args):
    def delta(x):
        return (function(*args, x) - np.exp(x)) / np.exp(x)

    y_values = []
    for item in x_values:
        y_values.append(delta(item))

    # Plot the function
    plt.plot(x_values, y_values)
    plt.xlabel('x')
    plt.ylabel('δ(x)')
    plt.title("Plot of δ(x)")
    plt.grid(True)
    plt.show()

def taylor_exp(x):
    ln2 = np.log(2)
    if np.abs(x) > ln2:
        invert = False
        if x < 0:
            invert = True
        x = np.abs(x)

        u = int(x / ln2)
        r = x - u * ln2
        result = (np.power(2, u) * taylorExponential_back(terms=10, x_0=r))
        if invert:
            result = 1/result
    else:
        result = taylorExponential_back(terms=10, x_0=x)
    return result

def main():
    error = False
    if error:
        txt_write_error(term_values=[10, 50], x_0_values=[1, 10, 20, -10, -20])
        csv_write_error(term_values=[10, 50], x_0_values=[1, 10, 20, -10, -20])

    plot = True
    if plot:
        x_values = np.linspace(0, np.log(2), 100)
        plotting_delta(x_values, taylorExponential_back, 10)




if __name__ == "__main__":
    main()








