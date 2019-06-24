__author__ = "Daniel Winklehner"
__doc__ = "Find out if a number is a power of two"


def power_of_two(number):
    """
    Function that checks if the input value (data) is a power of 2
    (i.e. 2, 4, 8, 16, 32, ...)
    """
    res = 0

    while res == 0:

        res = number % 2
        number /= 2.0
        print("res: {}, data: {}".format(res, number))

        if number == 1 and res == 0:

            return True

    return False
