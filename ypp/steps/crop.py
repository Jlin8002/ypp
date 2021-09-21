import numpy as np
from PIL import Image
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

from ypp.pipeline.step import Step


class CropStep(Step):
    """A pipeline step used to crop and invert an image before processing.

    Config parameters:
    ----------
    name : str
        The name of the step to be used in the config file.
    interactive : bool
        If True, the user will be prompted for input.
    crop : bool
        If True, the image will be cropped.
    invert : bool
        If True, the image will be inverted.
    """

    DEFAULT_CONFIG = {
        "name": "crop",
        "interactive": False,
        "crop": False,
        "invert": False,
    }

    def __init__(self, directory=None, config=None, name=None, verbose=False):
        super().__init__(directory=directory, config=config, name=name, verbose=verbose)
        if config is None:
            config = CropStep.DEFAULT_CONFIG
        self.set_input_ext("fits")
        self.set_output_ext("fits")
        self.set_output_tag("crop")


def invert(data: np.ndarray) -> np.ndarray:
    """Invert the pixel values of an image by subtracting from the maximum possible value.

    Parameters
    ----------
    data : np.ndarray
        The input image.

    Returns
    -------
    np.ndarray
        The inverted image.

    Raises
    ------
    TypeError
        If `data` is not an 8-bit or 16-bit integer array.
    """

    dt = data.dtype.type
    if (dt is np.int8) or (dt is np.uint8):
        bits = 8
    elif (dt is np.int16) or (dt is np.uint16):
        bits = 16
    else:
        raise TypeError("Image data is not 8-bit or 16-bit.")

    return ((1 << bits) - 1) - data


def crop_ndarray(
    data: np.ndarray, position: tuple((int, int)), size: tuple((int, int),)
) -> np.ndarray:
    """Crop a 2D numpy array.

    Parameters
    ----------
    data : np.ndarray
        The numpy array to crop.
    position : tuple
        The center pixel of the crop in (x, y) order.
    size : tuple
        The number of pixels in each direction that the crop should extend from the center pixel.

    Returns
    -------
    np.ndarray
        The cropped numpy array.
    """

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
    """Crop a FITS image.

    Parameters
    ----------
    filename : str
        The name of the FITS file to crop.
    position : tuple
        The center pixel of the crop in (x, y) order.
    size : tuple
        The number of pixels in each direction that the crop should extend from the center pixel.
    new_filename : str
        The name to be given to the cropped FITS file.
    """

    # Load the image and the WCS
    with fits.open(filename) as hdul:
        hdu = hdul[0]
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
