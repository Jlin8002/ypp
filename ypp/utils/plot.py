import numpy as np
import matplotlib as plt

# work on this to make it look prettier
def display(image):
    plt.figure()
    plt.imshow(image, cmap="gray")
    plt.colorbar()


def sex_plot(image):
    from matplotlib.patches import Ellipse

    data_sub = np.copy(image)

    # plot background-subtracted image
    fig, ax = plt.subplots(figsize=(12, 12))
    m, s = np.mean(data_sub), np.std(data_sub)
    im = ax.imshow(
        data_sub,
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
