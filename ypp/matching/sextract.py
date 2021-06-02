import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import sep
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u


def sextract_pix_locs(data: np.ndarray, sigma: float):
    data = data.astype("float")
    bkg = sep.Background(data)
    data = data - bkg
    objects = sep.extract(data, sigma, err=bkg.globalrms)
    return np.array([objects["x"], objects["y"]]).T


def pix_to_wcs(pix_locs, wcs: WCS):
    object_coords = wcs.pixel_to_world_values(pix_locs)
    return SkyCoord(object_coords, unit=u.deg, frame="icrs")


def sextract_world_locs(data, sigma, wcs):
    object_locs = sextract_pix_locs(data, sigma)
    object_coords = pix_to_wcs(object_locs, wcs)
    return object_coords


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
