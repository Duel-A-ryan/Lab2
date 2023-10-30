from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import plotting_functions as pf
from scipy.signal import find_peaks
import scipy.optimize as sc

if __name__ == '__main__':

    # Opening fits files for Neon and Darks and Flats
    neon_30s = fits.open("data/part3/Raw Data/Neon30sfit.fit")
    darks_30s = fits.open("data/part3/Darks/Darks__30s_01.fit")

    # Opening fits files for raw stellar data
    mirkak_raw_50 = fits.open("data/part3/Raw Data/Mirfak_SpectroPicture_50s.fits")
    jupiter_raw_5 = fits.open("data/part3/Raw Data/Jupiter5s.fit")
    moon_raw_1 = fits.open("data/part3/Raw Data/Moon_SpectroPicture_1s.fits")
    dark_15 = fits.open("data/part3/Darks/Dark-00015s.fit")
    flats_5 = fits.open("data/part3/Flats/Flats_5s_01.fit")

    # Cleaning Data
    neon_30s[0].data = pf.cleaning(neon_30s[0].data, 1, darks_30s[0].data, 30, flats_5[0].data, 5)
    mirkak_raw_50[0].data = pf.cleaning(mirkak_raw_50[0].data, 50, darks_30s[0].data, 30, flats_5[0].data, 5)
    moon_raw_1[0].data = pf.cleaning(moon_raw_1[0].data, 1, dark_15[0].data, 15, flats_5[0].data, 5)
    jupiter_raw_5[0].data = pf.cleaning(jupiter_raw_5[0].data, 5, dark_15[0].data, 15, flats_5[0].data, 5)

    # New fits files with Cleaned Data
    mirkak_raw_50.writeto("data/part3/Cleaned Data/mirfak_test.fits", overwrite=True)
    neon_30s.writeto("data/part3/Cleaned Data/Neon_data.fits", overwrite=True)
    moon_raw_1.writeto("data/part3/Cleaned Data/moon_test.fits", overwrite=True)
    jupiter_raw_5.writeto("data/part3/Cleaned Data/jupiter_test.fits", overwrite=True)

    # Plotting Neon Data
    neon_x, neon_y = np.loadtxt("data/part3/Cleaned Data/Neon_data", unpack=True)
    pf.line(True, neon_x, neon_y)
    peak_pixels, peak_heights = find_peaks(neon_y, height=0.9)
    pf.scatter(False, peak_pixels, peak_heights["peak_heights"])
    pf.show("Plots/Intensity Plots/3_Neon")

    wavelengths = np.array([540.056, 640.225, 650.653, 653.288, 659.885, 656.5, 692.647, 703.241, 717.395])

    # TODO: Something wrong with the uncertainties. Check this out
    # sc.curve_fit(pf.line_eq)

    """
    moon_x, moon_y = np.loadtxt("data/moon_test", unpack=True)
    jupiter_x, jupiter_y = np.loadtxt("data/jupiter_test", unpack=True)

    plt.plot(moon_x, moon_y / 5, label="Moon")
    plt.plot(jupiter_x, jupiter_y, label="Jupiter")
    plt.show()
    """
