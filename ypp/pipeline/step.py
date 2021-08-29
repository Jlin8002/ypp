import os
import numpy as np
from astropy.io import fits
import pandas as pd
from PIL import Image
import json

from utils.fits import Fits


class Step(object):
    """
    A base pipeline step class that can be inherited from.
    """

    def __init__(self, directory=None, config=None):
        self.directory = directory
        self.input_ext = None
        self.input_data = None
        self.output_ext = None
        self.output_tag = None
        self.output_data = None

        if type(config) == str:
            with open(config, "r") as f:
                self.config = json.load(f)
        elif type(config) == dict:
            self.config = config

    def set_input_ext(self, file_ext):
        self.input_file_ext = file_ext

    def set_output_ext(self, file_ext):
        self.output_file_ext = file_ext

    def set_output_tag(self, tag):
        self.output_file_tag = tag

    def set_directory(self, directory):
        self.directory = directory

    def read_fits_file(self, filename):
        filepath = os.path.join(self.directory, filename)
        data = fits.getdata(filepath)
        header = fits.getheader(filepath)
        self.input_data = Fits(filepath, data, header)

    def write_fits_file(self, filename):
        filepath = os.path.join(self.directory, filename)
        hdu = fits.PrimaryHDU(data=self.input_data.data, header=self.input_data.header)
        hdu.writeto(filepath, overwrite=True)

    def read_tiff_file(self, filename):
        filepath = os.path.join(self.directory, filename)
        with Image.open(filepath) as img:
            img = Image.open(filepath)
            self.input_data = np.array(img)

    def write_tiff_file(self, filename):
        filepath = os.path.join(self.directory, filename)
        img = Image.fromarray(self.output_data)
        img.save(filepath)

    def read_csv_file(self, filename):
        filepath = os.path.join(self.directory, filename)
        self.input_data = pd.read_csv(filepath)

    def write_csv_file(self, filename):
        filepath = os.path.join(self.directory, filename)
        self.output_data.to_csv(filepath)

    def read_data(self, filename):
        if self.input_file_ext == "fits" or self.input_file_ext == "fit":
            self.read_fits_file(filename)
        elif self.input_file_ext == "tif" or self.input_file_ext == "tiff":
            self.read_tiff_file(filename)
        elif self.input_file_ext == "csv":
            self.read_csv_file(filename)
        else:
            raise NameError("Unrecognized file extension.")

    def write_data(self, filename):
        if self.output_file_ext == "fits" or self.output_file_ext == "fit":
            self.write_fits_file()
        elif self.output_file_ext == "tif" or self.output_file_ext == "tiff":
            self.write_tiff_file()
        elif self.output_file_ext == "csv":
            self.write_csv_file()
        else:
            raise NameError("Unrecognized file extension.")

    def process_data(self):
        pass

    def run(self):
        self.process_data()

