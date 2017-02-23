#!/usr/bin/env python3
import numpy as np


def find_rowcolumn_of_LTri(number):
    """
    Return the row number and column number of Nth entry  of a lower triangle matrix.
    NUMBER, ROW, COLUMN are counted from ZERO!!!
    Example:
    0
    1  2
    3  4  5
    6  7  8  9
    10 11 12 13 14
    15 16 17 18 19 20

    >>> findrowcolumn(18)
    (5,3)
    # 18(19th) is at row 5(6th) and column 3(4th)
    """
    number += 1
    y = int((np.sqrt(1 + 8 * number) - 1) / 2)
    b = int(number - (y**2 + y) / 2)
    if b == 0:
        return (y - 1, y - 1)
    else:
        return (y, b - 1)


def find_num_of_LTri(i, j):
    if i < j:
        i, j = j, i
    num = i * (i - 1) / 2 + j
    num = int(num)
    return num - 1
