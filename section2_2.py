from astropy.io import fits
import matplotlib.pyplot as plt
import plotting_functions as pf
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import numpy as np

if __name__ == '__main__':
    """
    ----------------------------------------------------------------------------------------------------------------    
    # Opening the TelluricOH.dat file and reading its contents into two arrays.
      pixel contains the pixels_of_peaks of the plots and intensity contains the intensity found at those pixels_of_peaks  
    # 
    --------------------------------------------------------------------------------------------------------------    
    """

    pixel, intensity = np.loadtxt("data/TelluricOH_Data/TelluricOH_2", unpack=True)

    peak_pixels, peak_heights = find_peaks(intensity, height=1700)

    pf.line(True, pixel, intensity, "TITLE", "X", "Y", "UGH", None, "black")
    pf.scatter(False, peak_pixels, peak_heights['peak_heights'], color="red", style='.')

    wavelengths = np.array([16079.753, 16128.608, 16194.615, 16235.376, 16317.161, 16350.650, 16388.492, 16442.155,
                            16475.648, 16512.365, 16513, 16525, 16553.814])

    popt, pcov = curve_fit(pf.line_eq, peak_pixels, wavelengths)
    m, c = popt[0], popt[1]
    u_m, u_c = np.sqrt(np.diag(pcov))

    print(f"The wavelength equation is y = {popt[0]:.5}x + {popt[1]:.5}")
    print(f"The uncertainties of the slope is {u_m:.5} and the uncertainty in the intercept is {u_c:.5}")

    """Calculating Residuals"""
    # Find the predicted values of the max_pixels and then gather the residuals
    predicted_wavelengths = [pf.line_eq(p, m, c) for p in peak_pixels]
    residual_wavelengths = predicted_wavelengths - wavelengths

    # Plot the residuals to see the variation in the actual values and the predicted values from the line equation
    pf.scatter(True, np.arange(0, len(residual_wavelengths)), residual_wavelengths, label="Residuals of Gathered Data")
    plt.hlines(0, xmin=-1, xmax=14, colors="green")
    plt.legend()
    # plt.show()

    """
    --------------------------------------------------------------------------------------------------------------------
    # Reading the data from the Fe II topic to determine the wavelength
    --------------------------------------------------------------------------------------------------------------------
    """

    fe_pixel, fe_intensity = np.loadtxt("data/TelluricOH_Data/Fe_Data", unpack=True)  # COULD TOTALLY TRY FOR BETTER

    pf.line(True, fe_pixel, fe_intensity, "TITLE", "X", "Y", "LABEL", color="black")
    plt.show()

    max_intensity = max(fe_intensity)
    max_elem = np.where(fe_intensity == max_intensity)[0][0]

    fe_wave = pf.line_eq(fe_pixel[max_elem], m, c)
    print(fe_wave)


