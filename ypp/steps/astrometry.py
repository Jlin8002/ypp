import numpy as np
from astropy.io import fits
from astroquery.astrometry_net import AstrometryNet

from pipeline.step import Step
from utils.fits import Fits

TEMP_FILE_NAME = "temp.fits"


class AstrometryStep(Step):
    def __init__(self, directory=None, config=None, name=None, verbose=False):
        super().__init__(directory=directory)
        self.set_input_ext("fits")
        self.set_output_ext("fits")
        self.set_output_tag("astrometry")
        self.api_key = None

    def set_api_key(self, api_key):
        self.api_key = api_key

    def process_data(self):
        """
        Run astrometry.net on the input data.
        """
        self.write_data(TEMP_FILE_NAME)
        header = solve_image(TEMP_FILE_NAME, self.api_key)
        self.output_data = Fits(
            name=self.input_data.name, data=np.copy(self.input_data.data), header=header
        )


def solve_image(filename, api_key):
    ast = AstrometryNet()
    ast.api_key = api_key

    wcs_header = ast.solve_from_image(filename)

    return wcs_header
