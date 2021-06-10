import numpy as np
from scipy.spatial import distance as dt


def quick_match(sextracted_coords, catalog_coords, bin_size, tolerance):
    """
    quick_match: implementation of TOPCAT matching (bin-based approach)
        NOTE: ensure that all inputs use the same units (pixels or degrees)

        sextracted_coords: array/list of tuples ('x' and 'y' columns from sep.extract in [x, y] form)
        catalog_coords: array/list of tuples ('x' and 'y' columns from a star catalog in [x, y] form)
        bin_size: float (width of each bin)
        tolerance: float (search radius around an object)

        return: array of ints (elements are indexes of closest objects from catalog_coords, array indexed by sextracted_coords)
    """
    minima = np.amin(sextracted_coords, axis=0)
    maxima = np.amax(sextracted_coords, axis=0)
    #    print(minima, maxima)
    width_x = maxima[0] - minima[0]
    width_y = maxima[1] - minima[1]
    #    print(width_x, width_y)
    bin_num_x = np.ceil(width_x / bin_size).astype(int)
    bin_num_y = np.ceil(width_y / bin_size).astype(int)
    #    print(bin_num_x, bin_num_y)
    bins = [[[] for j in range(bin_num_y + 2)] for i in range(bin_num_x + 2)]
    bins_indices = [[[] for j in range(bin_num_y + 2)] for i in range(bin_num_x + 2)]
    sextracted_x, sextracted_y = (
        [loc[0] for loc in sextracted_coords],
        [loc[1] for loc in sextracted_coords],
    )
    catalog_x, catalog_y = (
        [loc[0] for loc in catalog_coords],
        [loc[1] for loc in catalog_coords],
    )
    bin_range_x = np.linspace(minima[0], maxima[0], bin_num_x + 1)
    bin_range_y = np.linspace(minima[1], maxima[1], bin_num_y + 1)
    #    print(bin_range_x, bin_range_y)
    sextracted_x_bins = np.digitize(sextracted_x, bin_range_x)
    sextracted_y_bins = np.digitize(sextracted_y, bin_range_y)
    catalog_x_bins = np.digitize(catalog_x, bin_range_x)
    catalog_y_bins = np.digitize(catalog_y, bin_range_y)
    for i in range(len(catalog_coords)):  # binning catalog stars
        #        print(i, catalog_x_bins[i], catalog_y_bins[i])
        bins[catalog_x_bins[i]][catalog_y_bins[i]].append(catalog_coords[i])
        bins_indices[catalog_x_bins[i]][catalog_y_bins[i]].append(i)
    bin_x = -1
    bin_y = -1

    # recursive_match: returns catalog index of closest match to a given set of coordinates
    def recursive_match(current_index, cm_distance, cm_index, center_distance):
        # check borders of square with side length ((center_distance * 2) + 1) centered on bin_x, bin_y
        closest_match_distance = cm_distance
        closest_match_index = cm_index
        first_match_found = 0
        if closest_match_distance != -1:
            first_match_found = 1
        for j in range(
            np.max([0, bin_x - center_distance]),
            np.min([bin_num_x + 2, bin_x + center_distance + 1]),
        ):
            for k in range(
                np.max([0, bin_y - center_distance]),
                np.min([bin_num_y + 2, bin_y + center_distance + 1]),
            ):
                if (
                    (np.abs(bin_x - j) == center_distance)
                    or (np.abs(bin_y - k) == center_distance)
                ) and (len(bins[j][k]) > 0):
                    #                    print(current_index, j, k)
                    #                    print(bins[j][k])
                    match_distances = dt.cdist(
                        bins[j][k], [sextracted_coords[current_index]]
                    )
                    #                    print(match_distances)
                    match_index_bin = np.argmin(match_distances)
                    if (match_distances[match_index_bin] < closest_match_distance) or (
                        closest_match_distance == -1
                    ):
                        closest_match_distance = match_distances[match_index_bin]
                        closest_match_index = bins_indices[j][k][match_index_bin]
        #                        print("UPDATE: " + str(closest_match_distance))
        #                        print(catalog_coords[closest_match_index])
        if closest_match_distance == -1:  # if no match found
            if ((center_distance + 1) * bin_size) > tolerance:  # if outside tolerance
                return -1
            else:  # if inside tolerance
                return recursive_match(
                    current_index,
                    closest_match_distance,
                    closest_match_index,
                    center_distance + 1,
                )
        else:  # if match found
            if (first_match_found == 0) or (
                closest_match_distance > (bin_size * (center_distance))
            ):  # check adjacent bins after first match found
                return recursive_match(
                    current_index,
                    closest_match_distance,
                    closest_match_index,
                    center_distance + 1,
                )
            elif (
                closest_match_distance > tolerance
            ):  # check if closest match is within tolerance
                return -1
        return closest_match_index

    index_list = []
    for i in range(len(sextracted_coords)):
        bin_x = sextracted_x_bins[i]
        bin_y = sextracted_y_bins[i]
        #        print("Index: " + str(i))
        index_list.append(recursive_match(i, -1, -1, 0))
    return np.array(index_list)


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
