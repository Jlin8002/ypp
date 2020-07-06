import sys
import os
import numpy as np
from scipy.spatial import distance as dt
import pandas as pd
import matplotlib.pyplot as plt
import SciServer
from SciServer import CasJobs     # Communicate between SciServer Compute and CasJobs
from astropy.io import fits, ascii
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astroquery.gaia import Gaia
from astropy.table import Table, Column
import sep

class Pipeline:
    steps = ['INIT', 'CROP', 'QUERY', 'SEXTRACT', 'MATCH', 'CALIBRATE']
    def __init__(self, image_list):
        self.images = glob.glob(image_list)
        self.current_step = None
        self.current_image = None

    def run_image():
        self.current_step = 'INIT'
        current_plate = Plate()
        while not(self.current_step == None):
            self.run_step()



class Plate:
    def __init__(self, path):
        self.path = path
        self.hdul = fits.open(path)
        if self.hdul[0].header['NAXIS'] == 3: # todo: find a way to deal with RGB images
            self.hdul[0].data = self.hdul[0].data[0]
            self.hdul[0].header['NAXIS'] = 2
        self.data_original = np.copy(self.hdul[0].data)
        self.wcs = WCS(self.hdul[0].header)
        self.bits = self.hdul[0].header['BITPIX']
        self.data = np.copy(self.data_original)

        self.crop_center = None
        self.crop_radius = None
        self.inverted = False
        self.center = (np.shape(self.data)[1] / 2, np.shape(self.data)[0] / 2) 
        # todo: look at SIP distortion to figure out center
        self.wcs_bounds = self.wcs.all_pix2world(np.asarray([(0, 0), (len(self.data[0]), 0), (0, len(self.data)), (len(self.data[0]), len(self.data))]), 1)
        self.pix_per_deg = len(self.data) / np.hypot(self.wcs_bounds[0][0] - self.wcs_bounds[1][0], self.wcs_bounds[0][1] - self.wcs_bounds[1][1])

        self.bkg = None
        self.objects = None
        self.catalog = None

        self.fluxes = None
        self.nan_mask = None
        self.mag_plate = None
        self.table = None

    # todo: add steps for astrometry.net and neural network

    def invert(self, display=True):
        self.data = ((1 << self.bits) - 1) - self.data

    def crop(self, center, radius, display=True):
        x = center[0]
        y = center[1]
        self.crop_center = center
        self.crop_radius = radius
        cutout = Cutout2D(self.data_original, center, 2*radius, WCS(self.hdul[0].header), copy=True)
        self.data = cutout.data
        self.wcs = cutout.wcs
        self.center = (x - self.center[0], y - self.center[1])
        self.wcs_bounds = self.wcs.all_pix2world(np.asarray([(0, 0), (len(self.data[0]), 0), (0, len(self.data)), (len(self.data[0]), len(self.data))]), 1)

    def query(self, catalog='SDSS', threshold=20, display=True):
        ra_bounds = []
        dec_bounds = []
        for bound in self.wcs_bounds:
            ra_bounds.append(bound[0])
            dec_bounds.append(bound[1])
        if catalog == 'SDSS':
            pandas_table = SDSS_query(np.amin(ra_bounds), np.amax(ra_bounds), np.amin(dec_bounds), np.amax(dec_bounds), threshold=threshold, num_stars=40000)
            self.catalog = Table.from_pandas(pandas_table)
            self.catalog.sort('g')
            print("SDSS stars found: " + str(len(self.catalog)))
        elif catalog == 'gaia':
            self.catalog = gaia_query_box(np.amin(ra_bounds), np.amax(ra_bounds), np.amin(dec_bounds), np.amax(dec_bounds), threshold=threshold)
            self.catalog.sort('phot_bp_mean_mag')
            print("Gaia stars found: " + str(len(self.catalog)))

    def sextract(self, sigma=3, apt_radius=0.004, buf_radius=0.006, bkg_radius=0.008, zero_point=22.5, display=True):
        self.bkg = sep.Background(self.data.astype('float'))
        self.data = self.data - self.bkg
        self.objects = sep.extract(self.data.astype('float'), sigma, err=self.bkg.globalrms)

        aperture_radius = apt_radius * self.pix_per_deg
        buffer_radius = buf_radius * self.pix_per_deg
        background_radius = bkg_radius * self.pix_per_deg
        aperture_area = (np.pi * (aperture_radius ** 2))
        buffer_area = (np.pi * (buffer_radius ** 2))
        background_area = (np.pi * (background_radius ** 2)) - buffer_area
        constant = zero_point
        print("Aperture radius: ", aperture_radius)
        print("Buffer radius: ", buffer_radius)
        print("Background radius: ", background_radius)

        aperture_list, _, _ = sep.sum_circle(self.data, self.objects['x'], self.objects['y'], aperture_radius, err=self.bkg.globalrms, gain=1.0)
        background_list, _, _ = sep.sum_circann(self.data, self.objects['x'], self.objects['y'], buffer_radius, background_radius, err=self.bkg.globalrms, gain=1.0)
        background_densities = [annulus_sum / background_area for annulus_sum in background_list]
        self.fluxes = np.array([aperture_list[i] - (aperture_area * background_densities[i]) for i in range(len(aperture_list))])
        self.mag_plate = constant - 2.5*np.log10(self.fluxes)


    def match(self, bin_size=0.1, tolerance=0.5, display=True):
        object_coords = np.array(self.wcs.all_pix2world([[self.objects['x'][i], self.objects['y'][i]] for i in range(len(self.objects))], 1))
        catalog_coords = [(self.catalog['ra'][i], self.catalog['dec'][i]) for i in range(len(self.catalog))]
        match_indices = quick_match(object_coords, catalog_coords, bin_size=0.1, tolerance=0.5)
        self.table = self.catalog[match_indices]
        self.table['flux'] = self.fluxes
        self.table['mag_plate'] = self.mag_plate
        self.table['plate_ra'] = object_coords[:, 0]
        self.table['plate_dec'] = object_coords[:, 1]
        self.table['plate_x'] = self.objects['x']
        self.table['plate_y'] = self.objects['y']
        nan_mask = self.fluxes < 0
        self.table.remove_rows(np.arange(0, len(nan_mask))[nan_mask])
        clean_mask = self.table['clean'] == 0
        self.table.remove_rows(np.arange(0, len(clean_mask))[clean_mask])


# work on this to make it look prettier
def display(image):
    plt.figure()
    plt.imshow(image, cmap='gray')
    plt.colorbar()

'''
quick_match: implementation of TOPCAT matching (bin-based approach)
    NOTE: ensure that all inputs use the same units (pixels or degrees)

    sextracted_coords: array/list of tuples ('x' and 'y' columns from sep.extract in [x, y] form)
    catalog_coords: array/list of tuples ('x' and 'y' columns from a star catalog in [x, y] form)
    bin_size: float (width of each bin)
    tolerance: float (search radius around an object)

    return: array of ints (elements are indexes of closest objects from catalog_coords, array indexed by sextracted_coords)
'''

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

def gaussian_box(image, coords, xrad, yrad):
    return fitgaussian(crop_box(image, coords, xrad, yrad), 0, False)

def gaia_query(ra_0, ra_1, dec_0, dec_1, threshold=20):
    jobquery = 'SELECT gaia_source.source_id,gaia_source.ra,gaia_source.dec,gaia_source.parallax,gaia_source.parallax_error,gaia_source.pmra,gaia_source.pmdec,gaia_source.phot_g_mean_mag,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.lum_val FROM gaiadr2.gaia_source  WHERE CONTAINS(POINT(\'ICRS\',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX(\'ICRS\',' + str((ra_0 + ra_1) / 2) + ',' + str((dec_0 + dec_1) / 2) + ',' + str(np.absolute(ra_1 - ra_0)) + ',' + str(np.absolute(dec_1 - dec_0)) + '))=1    AND  (phot_bp_mean_mag<= ' + str(threshold) + ')'
    job = Gaia.launch_job_async(jobquery, dump_to_file=True)
    return job.get_results()

def SDSS_query(ra_0, ra_1, dec_0, dec_1, threshold=20, num_stars=20000):
    jobquery = 'SELECT TOP ' + str(num_stars) + ' p.objid,p.ra,p.dec,p.u,p.g,p.r,p.i,p.z,p.type,p.clean,pm.pmra,pm.pmdec FROM PhotoObj AS p JOIN propermotions pm ON p.objid = pm.objid WHERE p.ra BETWEEN ' + str(ra_0) +  ' AND ' + str(ra_1) + ' AND p.dec BETWEEN ' + str(dec_0) + ' AND ' + str(dec_1) + ' AND p.g < ' + str(threshold)
    return CasJobs.executeQuery(sql=jobquery, context="DR15")

def gaussian_x(x, height, center, width, bgoffset):
#Returns a gaussian function with the given parameters
    width = float(width)
    return (height*np.exp(-((((x - center)/width)**2)/ 2.))) + bgoffset

def gauss_compare(dr, sl, g_dex):
    plt.plot(dr, sl[g_dex])
    plt.plot(dr, [gaussian_x(f, gaussian_fits[g_dex][0], gaussian_fits[g_dex][1], gaussian_fits[g_dex][2], gaussian_fits[g_dex][3]) for f in dr])

def sex_plot(image):
    from matplotlib.patches import Ellipse

    data_sub = np.copy(image)

    # plot background-subtracted image
    fig, ax = plt.subplots(figsize=(12,12))
    m, s = np.mean(data_sub), np.std(data_sub)
    im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
                vmin=m-s, vmax=m+s, origin='lower')

    # plot an ellipse for each object
    for i in range(len(objects)):
        e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                    width=6*objects['a'][i],
                    height=6*objects['b'][i],
                    angle=objects['theta'][i] * 180. / np.pi)
        e.set_facecolor('none')
        e.set_edgecolor('red')
        ax.add_artist(e)