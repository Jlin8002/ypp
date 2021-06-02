import numpy as np
from astroquery.sdss import SDSS
from astroquery.gaia import Gaia


def get_plate_bounds(data, wcs):
    coord_bounds = [
        (0, 0),
        (len(data[0]), 0),
        (0, len(data)),
        (len(data[0]), len(data)),
    ]
    wcs_bounds = wcs.all_pix2world(np.asarray(coord_bounds), 1)
    ra_bounds = (np.amin(wcs_bounds[:, 0]), np.amax(wcs_bounds[:, 0]))
    dec_bounds = (np.amin(wcs_bounds[:, 1]), np.amax(wcs_bounds[:, 1]))
    return ra_bounds, dec_bounds


def query(catalog="SDSS", threshold=20, display=True):
    if catalog == "SDSS":
        pass
    elif catalog == "gaia":
        pass


def gaia_query(ra_0, ra_1, dec_0, dec_1, threshold=20):
    jobquery = (
        "SELECT gaia_source.source_id,gaia_source.ra,gaia_source.dec,gaia_source.parallax,gaia_source.parallax_error,gaia_source.pmra,gaia_source.pmdec,gaia_source.phot_g_mean_mag,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.lum_val FROM gaiadr2.gaia_source  WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS',"
        + str((ra_0 + ra_1) / 2)
        + ","
        + str((dec_0 + dec_1) / 2)
        + ","
        + str(abs(ra_1 - ra_0))
        + ","
        + str(abs(dec_1 - dec_0))
        + "))=1    AND  (phot_bp_mean_mag<= "
        + str(threshold)
        + ")"
    )
    job = Gaia.launch_job_async(jobquery, dump_to_file=True)
    return job.get_results()


def SDSS_query(ra_bounds, dec_bounds, threshold=20, num_stars=20000):
    jobquery = (
        "SELECT TOP "
        + str(num_stars)
        + " p.objid,p.ra,p.dec,p.u,p.g,p.r,p.i,p.z,p.type,p.clean,pm.pmra,pm.pmdec FROM PhotoObj AS p JOIN propermotions pm ON p.objid = pm.objid WHERE p.ra BETWEEN "
        + str(ra_bounds[0])
        + " AND "
        + str(ra_bounds[1])
        + " AND p.dec BETWEEN "
        + str(dec_bounds[0])
        + " AND "
        + str(dec_bounds[1])
        + " AND p.g < "
        + str(threshold)
    )
    res = SDSS.query_sql(jobquery, data_release=15)
    return res
