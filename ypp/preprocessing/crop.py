import numpy as np
from PIL import Image

from ypp.pipeline.step import PipeStep

class CropStep(PipeStep):
    def __init__(self, directory=None):
        super().__init__(directory=directory)
        super().set_input_exts('fits', 'crop')


def invert(data: np.ndarray) -> np.ndarray:
    dt = data.dtype.type
    if (dt is np.int8) or (dt is np.uint8):
        bits = 8
    elif (dt is np.int16) or (dt is np.uint16):
        bits = 16
    else:
        raise TypeError("Image data is not 8-bit or 16-bit.")

    return ((1 << bits) - 1) - data


def crop(data: np.ndarray, center: tuple((int, int)), radius: tuple((int, int))) -> np.ndarray:
    max_row, max_col = data.shape
    cen_row, cen_col = center
    rad_row, rad_col = radius
    row_lower_bound = max(cen_row - rad_row, 0)
    row_upper_bound = min(cen_row + rad_row, max_row)
    col_lower_bound = max(cen_col - rad_col, 0)
    col_upper_bound = min(cen_col + rad_col, max_col)

    return data[row_lower_bound:row_upper_bound, col_lower_bound:col_upper_bound]
