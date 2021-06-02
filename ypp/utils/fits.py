import numpy as np
from astropy.io import fits


def ndarray_to_fits(data: np.ndarray, filename: str):
    hdu = fits.PrimaryHDU(data=data)
    hdu.writeto(filename, overwrite=True)
