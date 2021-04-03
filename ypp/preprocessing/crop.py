import numpy as np
from PIL import Image

from typing import Tuple


def invert(data: np.ndarray) -> np.ndarray:
    datatype = data.dtype

    data = ((1 << 16) - 1) - data


def crop(data: np.ndarray, center: Tuple(int, int), radius: float) -> np.ndarray:
    x = center[0]
    y = center[1]
