import numpy as np
from PIL import Image
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

from ypp.pipeline.step import PipeStep


class CropStep(PipeStep):
    def __init__(self, directory=None):
        super().__init__(directory=directory)
        super().set_input_exts("fits", "crop")


def invert(data: np.ndarray) -> np.ndarray:
    dt = data.dtype.type
    if (dt is np.int8) or (dt is np.uint8):
        bits = 8
    elif (dt is np.int16) or (dt is np.uint16):
        bits = 16
    else:
        raise TypeError("Image data is not 8-bit or 16-bit.")

    return ((1 << bits) - 1) - data


def crop_ndarray(
    data: np.ndarray, position: tuple((int, int)), size: tuple((int, int))
) -> np.ndarray:
    max_row, max_col = data.shape
    cen_row, cen_col = position
    rad_row, rad_col = size
    row_lower_bound = max(cen_row - rad_row, 0)
    row_upper_bound = min(cen_row + rad_row, max_row)
    col_lower_bound = max(cen_col - rad_col, 0)
    col_upper_bound = min(cen_col + rad_col, max_col)

    return data[row_lower_bound:row_upper_bound, col_lower_bound:col_upper_bound]


def crop_fits(
    filename: str,
    position: tuple((int, int)),
    size: tuple((int, int)),
    new_filename: str,
) -> None:

    # Load the image and the WCS
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

    # Make the cutout, including the WCS
    cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)

    # Put the cutout image in the FITS HDU
    hdu.data = cutout.data

    # Update the FITS header with the cutout WCS
    hdu.header.update(cutout.wcs.to_header())

    # Write the cutout to a new FITS file
    cutout_filename = new_filename
    hdu.writeto(cutout_filename, overwrite=True)
