import plotting_functions as pf
import matplotlib.pyplot as plt
import numpy as np


if __name__ == '__main__':
    """
    # Opening the two files to write their data into numpy arrays
    # The first file is the spectrum of a blackbody source (Group_X_BB.dat)
    # The second file is the Neon calibration (Ne_calib.dat) which contains the 
      spectrum of a Neon lamp using a spectrograph.
    """

    bb_data = pf.readfile("data/Group_X_BB.dat")
    ne_data = pf.readfile("data/Ne_calib.dat")

    bb_pixel = np.arange(0, len(bb_data))
    ne_pixel = np.arange(0, len(ne_data))

    pf.scatter(True, bb_pixel, bb_data, "test", "xtest", "ytest", "label", ".")
    pf.line(True, ne_pixel, ne_data, "test", "xtest", "ytest", "label", None)

    plt.legend()
    plt.show()
