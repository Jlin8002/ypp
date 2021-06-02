import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def plot_plate(data, figsize=None, logscale=True):
    if figsize:
        plt.figure(figsize=figsize)
    if logscale:
        plt.imshow(data, cmap="gray", norm=LogNorm())
    else:
        plt.imshow(data, cmap="gray")
    plt.colorbar()
