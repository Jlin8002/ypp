def query(self, catalog="SDSS", threshold=20, display=True):

    ra_bounds = []
    dec_bounds = []
    for bound in self.wcs_bounds:
        ra_bounds.append(bound[0])
        dec_bounds.append(bound[1])
    if catalog == "SDSS":
        self.catalog = Table.from_pandas(pandas_table)
        self.catalog.sort("g")
        print("SDSS stars found: " + str(len(self.catalog)))
    elif catalog == "gaia":
        self.catalog = gaia_query(
            np.amin(ra_bounds),
            np.amax(ra_bounds),
            np.amin(dec_bounds),
            np.amax(dec_bounds),
            threshold=threshold,
        )
        self.catalog.sort("phot_bp_mean_mag")
        print("Gaia stars found: " + str(len(self.catalog)))


def gaia_query(ra_0, ra_1, dec_0, dec_1, threshold=20):
    jobquery = (
        "SELECT gaia_source.source_id,gaia_source.ra,gaia_source.dec,gaia_source.parallax,gaia_source.parallax_error,gaia_source.pmra,gaia_source.pmdec,gaia_source.phot_g_mean_mag,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.lum_val FROM gaiadr2.gaia_source  WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS',"
        + str((ra_0 + ra_1) / 2)
        + ","
        + str((dec_0 + dec_1) / 2)
        + ","
        + str(np.absolute(ra_1 - ra_0))
        + ","
        + str(np.absolute(dec_1 - dec_0))
        + "))=1    AND  (phot_bp_mean_mag<= "
        + str(threshold)
        + ")"
    )
    job = Gaia.launch_job_async(jobquery, dump_to_file=True)
    return job.get_results()


def SDSS_query(ra_0, ra_1, dec_0, dec_1, threshold=20, num_stars=20000):
    jobquery = (
        "SELECT TOP "
        + str(num_stars)
        + " p.objid,p.ra,p.dec,p.u,p.g,p.r,p.i,p.z,p.type,p.clean,pm.pmra,pm.pmdec FROM PhotoObj AS p JOIN propermotions pm ON p.objid = pm.objid WHERE p.ra BETWEEN "
        + str(ra_0)
        + " AND "
        + str(ra_1)
        + " AND p.dec BETWEEN "
        + str(dec_0)
        + " AND "
        + str(dec_1)
        + " AND p.g < "
        + str(threshold)
    )
    return CasJobs.executeQuery(sql=jobquery, context="DR15")
