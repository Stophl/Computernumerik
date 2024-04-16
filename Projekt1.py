
import numpy as np


def taylorExponential(terms=1, x_0=0):
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

def errorFinding(terms=1, x_0=0):
    summands, np_sum, forward_sum, backward_sum = taylorExponential(terms=terms, x_0=x_0)
    true_value = np.exp(x_0)

    forward_abs_error = np.abs(forward_sum - true_value)
    forward_rel_error = forward_abs_error / true_value

    backward_abs_error = np.abs(backward_sum - true_value)
    backward_rel_error = backward_abs_error / true_value

    check_abs_error = np.abs(np_sum - true_value)
    check_rel_error = check_abs_error / true_value

    print('"True" Value:', f"{true_value:.3}")
    print("Forward Sum Value:", f"{forward_sum:.3}")
    print("Forward Absolute Error:", f"{forward_abs_error:.3}")
    print("Forward Relative Error:", f"{forward_rel_error:.3}")
    print("Backward Sum Value:", f"{backward_sum:.3}")
    print("Backward Absolute Error:", f"{backward_abs_error:.3}")
    print("Backward Relative Error:", f"{backward_rel_error:.3}")
    print("CheckSum Absolute Error:", f"{check_abs_error:.3}")
    print("CheckSum Relative Error:", f"{check_rel_error:.3}")

term_list = [10, 50]
x0_list = [1, 10, 20, -10, -20]
#errorFinding(10, 20)

for A in x0_list:
    for B in term_list:
        print("***** x0 = ",A,", n = ",B," *****")
        errorFinding(B, A)
        print("***** ***** ***** ***** *****")
        print("")







