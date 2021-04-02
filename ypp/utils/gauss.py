def gaussian_box(image, coords, xrad, yrad):
    return fitgaussian(crop_box(image, coords, xrad, yrad), 0, False)


def gaussian_x(x, height, center, width, bgoffset):
    # Returns a gaussian function with the given parameters
    width = float(width)
    return (height * np.exp(-((((x - center) / width) ** 2) / 2.0))) + bgoffset


def gauss_compare(dr, sl, g_dex):
    plt.plot(dr, sl[g_dex])
    plt.plot(
        dr,
        [
            gaussian_x(
                f,
                gaussian_fits[g_dex][0],
                gaussian_fits[g_dex][1],
                gaussian_fits[g_dex][2],
                gaussian_fits[g_dex][3],
            )
            for f in dr
        ],
    )
