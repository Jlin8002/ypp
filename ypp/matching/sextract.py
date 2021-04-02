class Plate:
    def __init__(self, path):
        self.path = path
        self.hdul = fits.open(path)
        if (
            self.hdul[0].header["NAXIS"] == 3
        ):  # todo: find a way to deal with RGB images
            self.hdul[0].data = self.hdul[0].data[0]
            self.hdul[0].header["NAXIS"] = 2
        self.data_original = np.copy(self.hdul[0].data)
        self.wcs = WCS(self.hdul[0].header)
        self.bits = self.hdul[0].header["BITPIX"]
        self.data = np.copy(self.data_original)
        self.datain = np.copy(self.data_original)
        self.dataout = None

        self.crop_center = None
        self.crop_radius = None
        self.inverted = False
        self.center = (np.shape(self.data)[1] / 2, np.shape(self.data)[0] / 2)
        # todo: look at SIP distortion to figure out center
        self.wcs_bounds = self.wcs.all_pix2world(
            np.asarray(
                [
                    (0, 0),
                    (len(self.data[0]), 0),
                    (0, len(self.data)),
                    (len(self.data[0]), len(self.data)),
                ]
            ),
            1,
        )
        self.wcs_center = self.wcs.all_pix2world(self.center, 1)
        self.pix_per_deg = len(self.data) / np.hypot(
            self.wcs_bounds[0][0] - self.wcs_bounds[1][0],
            self.wcs_bounds[0][1] - self.wcs_bounds[1][1],
        )

        self.bkg = None
        self.objects = None
        self.catalog = None

        self.fluxes = None
        self.nan_mask = None
        self.mag_plate = None
        self.table = None
