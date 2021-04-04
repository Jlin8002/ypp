import numpy as np
from PIL import Image

from typing import Tuple


def invert(data: np.ndarray) -> np.ndarray:
    if data.dtype is np.int8:
        bits = 8
    elif data.dtype is np.int16:
        bits = 16
    else:
        raise TypeError("Image data is not 8-bit or 16-bit.")

    return ((1 << bits) - 1) - data


def crop(data: np.ndarray, center: Tuple(int, int), radius: int) -> np.ndarray:
    max_row, max_col = data.shape
    cen_row, cen_col = center
    row_lower_bound = max(cen_row - radius, 0)
    row_upper_bound = min(cen_row + radius, max_row)
    col_lower_bound = max(cen_col - radius, 0)
    col_upper_bound = min(cen_col + radius, max_col)

    return data[row_lower_bound:row_upper_bound, col_lower_bound:col_upper_bound]
