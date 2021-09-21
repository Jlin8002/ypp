import numpy as np
from ypp.utils import query

from ypp.pipeline.step import Step


class MatchStep(Step):
    """A pipeline step used to crop and invert an image before processing.

    Config parameters:
    ----------
    name : str
        The name of the step to be used in the config file.
    interactive : bool
        If True, the user will be prompted for input.
    objects : str
        A file path to a csv file containing object coordinates. If None, the input will be piped in from the previous step.
    catalog : str
        A file path to a csv file containing an astronomical catalog. If None, SDSS will be queried.
    """

    DEFAULT_CONFIG = {
        "name": "match",
        "interactive": False,
        "objects": None,
        "catalog": None,
    }

    def __init__(self, directory=None, config=None, name=None, verbose=False):
        super().__init__(directory=directory, config=config, name=name, verbose=verbose)
        if config is None:
            config = MatchStep.DEFAULT_CONFIG
        self.set_input_ext("csv")
        self.set_output_ext("csv")
        self.set_output_tag("match")


def match(objects, catalog, wcs, bin_size=0.1, tolerance=0.5, display=True):
    object_coords = np.array(
        wcs.all_pix2world(
            [[objects["x"][i], objects["y"][i]] for i in range(len(objects))], 1,
        )
    )
    catalog_coords = [
        (catalog["ra"][i], catalog["dec"][i]) for i in range(len(catalog))
    ]
    match_indices = quick_match(
        object_coords, catalog_coords, bin_size=0.1, tolerance=0.5
    )
    table = catalog[match_indices]
    table["flux"] = fluxes
    table["mag_plate"] = mag_plate
    table["plate_ra"] = object_coords[:, 0]
    table["plate_dec"] = object_coords[:, 1]
    table["plate_x"] = objects["x"]
    table["plate_y"] = objects["y"]
    nan_mask = fluxes < 0
    table.remove_rows(np.arange(0, len(nan_mask))[nan_mask])
    clean_mask = table["clean"] == 0
    table.remove_rows(np.arange(0, len(clean_mask))[clean_mask])
