import sys
import os
import numpy as np
from scipy.spatial import distance as dt
import pandas as pd
import SciServer
from SciServer import CasJobs     # Communicate between SciServer Compute and CasJobs
from astroquery.gaia import Gaia

class Plate:
    def __init__(self, path)
        self.path = path
        self.hdul = fits.open(path)
        self.data = np.copy(hdul[0].data)
        if self.hdul[0].header['NAXIS'] == 3:
            self.data = self.data[0]
        self.bits = self.hdul[0].header['BITPIX']


# quick_match: implementation of TOPCAT matching (bin-based approach)
#     NOTE: ensure that all inputs use the same units (pixels or degrees)
#
#     sextracted_coords: array/list of tuples ('x' and 'y' columns from sep.extract in [x, y] form)
#     catalog_coords: array/list of tuples ('x' and 'y' columns from a star catalog in [x, y] form)
#     bin_size: float (width of each bin)
#     tolerance: float (search radius around an object)
#
#     return: array of ints (elements are indexes of closest objects from catalog_coords, array indexed by sextracted_coords)

def quick_match(sextracted_coords, catalog_coords, bin_size, tolerance):
    minima = np.amin(sextracted_coords, axis = 0)
    maxima = np.amax(sextracted_coords, axis = 0)
#    print(minima, maxima)
    width_x = maxima[0] - minima[0]
    width_y = maxima[1] - minima[1]
#    print(width_x, width_y)
    bin_num_x = np.ceil(width_x / bin_size).astype(int)
    bin_num_y = np.ceil(width_y / bin_size).astype(int)
#    print(bin_num_x, bin_num_y)
    bins = [[[] for j in range(bin_num_y + 2)] for i in range(bin_num_x + 2)]
    bins_indices = [[[] for j in range(bin_num_y + 2)] for i in range(bin_num_x + 2)]
    sextracted_x, sextracted_y = [loc[0] for loc in sextracted_coords], [loc[1] for loc in sextracted_coords]
    catalog_x, catalog_y = [loc[0] for loc in catalog_coords], [loc[1] for loc in catalog_coords]
    bin_range_x = np.linspace(minima[0], maxima[0], bin_num_x + 1)
    bin_range_y = np.linspace(minima[1], maxima[1], bin_num_y + 1)
#    print(bin_range_x, bin_range_y)
    sextracted_x_bins = np.digitize(sextracted_x, bin_range_x)
    sextracted_y_bins = np.digitize(sextracted_y, bin_range_y)
    catalog_x_bins = np.digitize(catalog_x, bin_range_x)
    catalog_y_bins = np.digitize(catalog_y, bin_range_y)
    for i in range(len(catalog_coords)): # binning catalog stars
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
        for j in range(np.max([0, bin_x - center_distance]), np.min([bin_num_x + 2, bin_x + center_distance + 1])):
            for k in range(np.max([0, bin_y - center_distance]), np.min([bin_num_y + 2, bin_y + center_distance + 1])):
                if (((np.abs(bin_x - j) == center_distance) or (np.abs(bin_y - k) == center_distance)) and (len(bins[j][k]) > 0)):
#                    print(current_index, j, k)
#                    print(bins[j][k])
                    match_distances = dt.cdist(bins[j][k], [sextracted_coords[current_index]])
#                    print(match_distances)
                    match_index_bin = np.argmin(match_distances)
                    if (match_distances[match_index_bin] < closest_match_distance) or (closest_match_distance == -1):
                        closest_match_distance = match_distances[match_index_bin]
                        closest_match_index = bins_indices[j][k][match_index_bin]
#                        print("UPDATE: " + str(closest_match_distance))
#                        print(catalog_coords[closest_match_index])
        if closest_match_distance == -1: # if no match found
            if ((center_distance + 1) * bin_size) > tolerance: # if outside tolerance
                return -1
            else: # if inside tolerance
                return recursive_match(current_index, closest_match_distance, closest_match_index, center_distance + 1)
        else: # if match found
            if ((first_match_found == 0) or (closest_match_distance > (bin_size * (center_distance)))): # check adjacent bins after first match found
                return recursive_match(current_index, closest_match_distance, closest_match_index, center_distance + 1)
            elif closest_match_distance > tolerance: # check if closest match is within tolerance
                return -1
        return closest_match_index
    
    index_list = []
    for i in range(len(sextracted_coords)):
        bin_x = sextracted_x_bins[i]
        bin_y = sextracted_y_bins[i]
#        print("Index: " + str(i))
        index_list.append(recursive_match(i, -1, -1, 0))
    return np.array(index_list)


def index_print(indexed_list):
    for i in range(len(indexed_list)):
        print(i, indexed_list[i])


def dict_print(dictionary):
    for key in dictionary:
        print(key, dictionary[key])

def crop(image, xa, ya, xb, yb): 
    cropped_image = np.copy(image)
    return cropped_image[ya:yb, xa:xb]

def crop_box(image, coords, xrad=20, yrad=20):
    return crop(image, coords[0] - xrad, coords[1] - yrad, coords[0] + xrad, coords[1] + yrad)

def gaussian_box(image, coords, xrad, yrad):
    return fitgaussian(crop_box(image, coords, xrad, yrad), 0, False)

def gaia_query_radius(loc, radius=0.2, threshold=20):
    jobquery = 'SELECT gaia_source.source_id,gaia_source.ra,gaia_source.dec,gaia_source.parallax,gaia_source.parallax_error,gaia_source.pmra,gaia_source.pmdec,gaia_source.phot_g_mean_mag,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.lum_val FROM gaiadr2.gaia_source  WHERE CONTAINS(POINT(\'ICRS\',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE(\'ICRS\',' + str(loc[0]) + ',' + str(loc[1]) + ',' + str((radius / 60)) + '))=1    AND  (phot_bp_mean_mag<= ' + str(threshold) + ')'
    job = Gaia.launch_job_async(jobquery, dump_to_file=True)
    return job.get_results()

def gaia_query_box(ra_0, ra_1, dec_0, dec_1, threshold=20):
    jobquery = 'SELECT gaia_source.source_id,gaia_source.ra,gaia_source.dec,gaia_source.parallax,gaia_source.parallax_error,gaia_source.pmra,gaia_source.pmdec,gaia_source.phot_g_mean_mag,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.lum_val FROM gaiadr2.gaia_source  WHERE CONTAINS(POINT(\'ICRS\',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX(\'ICRS\',' + str((ra_0 + ra_1) / 2) + ',' + str((dec_0 + dec_1) / 2) + ',' + str(np.absolute(ra_1 - ra_0)) + ',' + str(np.absolute(dec_1 - dec_0)) + '))=1    AND  (phot_bp_mean_mag<= ' + str(threshold) + ')'
    job = Gaia.launch_job_async(jobquery, dump_to_file=True)
    return job.get_results()

def SDSS_query(ra_0, ra_1, dec_0, dec_1, threshold=20, num_stars=20000):
    jobquery = 'SELECT TOP ' + str(num_stars) + ' p.objid,p.ra,p.dec,p.u,p.g,p.r,p.i,p.z,p.type,p.clean,pm.pmra,pm.pmdec FROM PhotoObj AS p JOIN propermotions pm ON p.objid = pm.objid WHERE p.ra BETWEEN ' + str(ra_0) +  ' AND ' + str(ra_1) + ' AND p.dec BETWEEN ' + str(dec_0) + ' AND ' + str(dec_1) + ' AND p.g < ' + str(threshold)
    return CasJobs.executeQuery(sql=jobquery, context="DR15")