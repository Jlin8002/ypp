import numpy as np
from astropy.io import fits
from astroquery.astrometry_net import AstrometryNet


def solve_image(filename, api_key):
    ast = AstrometryNet()
    ast.api_key = api_key

    wcs_header = ast.solve_from_image(filename)

    return wcs_header
