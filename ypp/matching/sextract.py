import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import sep
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u


def sextract_bkg(data):
    data = data.astype("float")
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


def plot_objects(image, objects):

    # plot background-subtracted image
    fig, ax = plt.subplots(figsize=(12, 12))
    m, s = np.mean(image), np.std(image)
    im = ax.imshow(
        image,
        interpolation="nearest",
        cmap="gray",
        vmin=m - s,
        vmax=m + s,
        origin="lower",
    )

    # plot an ellipse for each object
    for i in range(len(objects)):
        e = Ellipse(
            xy=(objects["x"][i], objects["y"][i]),
            width=6 * objects["a"][i],
            height=6 * objects["b"][i],
            angle=objects["theta"][i] * 180.0 / np.pi,
        )
        e.set_facecolor("none")
        e.set_edgecolor("red")
        ax.add_artist(e)

    plt.show()
