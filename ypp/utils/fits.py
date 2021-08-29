from dataclasses import dataclass
import numpy as np
from astropy.io import fits


@dataclass
class Fits:
    name: str
    data: np.ndarray
    header: fits.Header


def ndarray_to_fits(data: np.ndarray, filename: str):
    hdu = fits.PrimaryHDU(data=data)
    hdu.writeto(filename, overwrite=True)

