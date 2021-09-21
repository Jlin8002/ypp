import numpy as np
import matplotlib.pyplot as plt
import sep
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

from ypp.pipeline.step import Step


class SextractStep(Step):
    """A pipeline step used to crop and invert an image before processing.

    Config parameters:
    ----------
    name : str
        The name of the step to be used in the config file.
    interactive : bool
        If True, the user will be prompted for input.
    units : str
        The units to be used for the photometry arguments. Can be "pixel" or "arcsecond".
    aperture : float
        The radius of the object photometry aperture in the units specified above.
    buffer : float
        The radius of the object photometry buffer in the units specified above.
    background : float
        The radius of the object photometry background in the units specified above.
    """

    DEFAULT_CONFIG = {
        "name": "match",
        "interactive": False,
        "aperture": None,
        "buffer": None,
        "background": None,
    }

    def __init__(self, directory=None, config=None, name=None, verbose=False):
        super().__init__(directory=directory, config=config, name=name, verbose=verbose)
        if config is None:
            config = SextractStep.DEFAULT_CONFIG
        self.set_input_ext("fits")
        self.set_output_ext("csv")
        self.set_output_tag("sextract")

    def process_data(self):
        wcs = WCS(self.get_input_data().header)


def sextract_bkg(data):
    data = data.astype(float)
    bkg = sep.Background(data)
    data = data - bkg
    return data, bkg


def sextract_objects(data: np.ndarray, sigma: float, bkg=None):
    if bkg is None:
        data, bkg = sextract_bkg(data)
    objects = sep.extract(data, sigma, err=bkg.globalrms)
    return objects


def pix_to_wcs(objects, wcs: WCS):
    pix_locs = np.array([objects["x"], objects["y"]]).T
    # object_coords = wcs.pixel_to_world_values(pix_locs)
    object_coords = wcs.all_pix2world(pix_locs, 0)
    return SkyCoord(object_coords, unit=u.deg, frame="icrs")


def sextract_world_locs(data, sigma, wcs):
    objects = sextract_objects(data, sigma)
    object_coords = pix_to_wcs(objects, wcs)
    return object_coords


def sextract_magnitudes(
    data,
    aperture_radius,
    buffer_radius,
    background_radius,
    bkg=None,
    objects=None,
    sigma=3.0,
):
    aperture_area = np.pi * (aperture_radius ** 2)
    buffer_area = np.pi * (buffer_radius ** 2)
    background_area = (np.pi * (background_radius ** 2)) - buffer_area

    if bkg is None:
        data, bkg = sextract_bkg(data)
    if objects is None:
        objects = sextract_objects(data, sigma, bkg)

    aperture_list, _, _ = sep.sum_circle(
        data, objects["x"], objects["y"], aperture_radius, err=bkg.globalrms, gain=1.0
    )
    background_list, _, _ = sep.sum_circann(
        data,
        objects["x"],
        objects["y"],
        buffer_radius,
        background_radius,
        err=bkg.globalrms,
        gain=1.0,
    )
    background_densities = [
        annulus_sum / background_area for annulus_sum in background_list
    ]
    fluxes = np.array(
        [
            aperture_list[i] - (aperture_area * background_densities[i])
            for i in range(len(aperture_list))
        ]
    )
    return fluxes


# sep_table = Table()
# sep_table["x"] = objects["x"]
# sep_table["y"] = objects["y"]
# sep_table["ra"] = object_ra
# sep_table["dec"] = object_dec
# sep_table.write(
#     datapath + "/data_tables/" + extensionless_name + "_sep.csv",
#     format="csv",
#     overwrite=True,
# )

