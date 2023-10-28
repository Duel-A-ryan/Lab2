from astropy.io import fits
import numpy as np


if __name__ == '__main__':

    # Reading fits files into variables for data cleaning
    moon_1s = fits.open("data/part3/Moon_SpectroPicture_1s.fits")
    dark_25s = fits.open("data/part3/Darks/Dark-00025s.fit")
    dark_1s = fits.open("data/part3/Darks/Dark-0002_1s.fit")
    flats_5s = fits.open("data/part3/Flats_5s_01.fit")

    #print(mirfak_50s[0].data)
    #rint(dark_25s[0].data)

    print(moon_1s[0].data * 25 - dark_25s[0].data)

"""
from astropy.io import fits
import numpy as np

mirfak_raw = fits.open("data/part3/Mirfak_SpectroPicture_50s.fits")
dark_25 = fits.open("data/part3/Darks/Dark-00025s.fit")
dark_1 = fits.open("data/part3/Darks/Dark-0002_1s.fit")

# Finding the average dark photon data over different durations
avg_dark_25s = np.mean(dark_25[0].data)
avg_dark_1s = np.mean(dark_1[0].data)
"""