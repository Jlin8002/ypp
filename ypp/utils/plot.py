import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Ellipse


def plot_plate(data, figsize=None, logscale=True):
    if figsize:
        plt.figure(figsize=figsize)
    if logscale:
        plt.imshow(data, cmap="gray", norm=LogNorm())
    else:
        plt.imshow(data, cmap="gray")
    plt.colorbar()


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
