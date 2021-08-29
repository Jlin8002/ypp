import numpy as np
from ypp.utils import query


def match(self, bin_size=0.1, tolerance=0.5, display=True):
    object_coords = np.array(
        self.wcs.all_pix2world(
            [
                [self.objects["x"][i], self.objects["y"][i]]
                for i in range(len(self.objects))
            ],
            1,
        )
    )
    catalog_coords = [
        (self.catalog["ra"][i], self.catalog["dec"][i])
        for i in range(len(self.catalog))
    ]
    match_indices = quick_match(
        object_coords, catalog_coords, bin_size=0.1, tolerance=0.5
    )
    self.table = self.catalog[match_indices]
    self.table["flux"] = self.fluxes
    self.table["mag_plate"] = self.mag_plate
    self.table["plate_ra"] = object_coords[:, 0]
    self.table["plate_dec"] = object_coords[:, 1]
    self.table["plate_x"] = self.objects["x"]
    self.table["plate_y"] = self.objects["y"]
    nan_mask = self.fluxes < 0
    self.table.remove_rows(np.arange(0, len(nan_mask))[nan_mask])
    clean_mask = self.table["clean"] == 0
    self.table.remove_rows(np.arange(0, len(clean_mask))[clean_mask])
