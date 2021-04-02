def invert(self, display=True):
    self.data = ((1 << self.bits) - 1) - self.data


def crop(self, center, radius, display=True):
    x = center[0]
    y = center[1]
    self.crop_center = center
    self.crop_radius = radius
    cutout = Cutout2D(
        self.data_original, center, 2 * radius, WCS(self.hdul[0].header), copy=True
    )
    self.data = cutout.data
    self.wcs = cutout.wcs
    self.center = (x - self.center[0], y - self.center[1])
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
