
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

    print('"True" Value:', true_value)
    print("Forward Sum Value:", forward_sum)
    print("Forward Absolute Error:", forward_abs_error)
    print("Forward Relative Error:", forward_rel_error)
    print("Backward Sum Value:", backward_sum)
    print("Backward Absolute Error:", backward_abs_error)
    print("Backward Relative Error:", backward_rel_error)
    print("CheckSum Absolute Error:", check_abs_error)
    print("CheckSum Relative Error:", check_rel_error)

errorFinding(10, 20)







