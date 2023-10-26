from astropy.io import fits
import matplotlib.pyplot as plt
import plotting_functions as pf
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import numpy as np

if __name__ == '__main__':
    """
    ----------------------------------------------------------------------------------------------------------------    
    # Opening the TelluricOH.dat file and reading its contents into two arrays. Variable pixel contains the 
      pixels_of_peaks of the plots and intensity contains the intensity found at those pixels_of_peaks  
    # 
    --------------------------------------------------------------------------------------------------------------    
    """

    # Read and write the .dat file into two arrays
    pixel, intensity = np.loadtxt("data/TelluricOH_Data/TelluricOH_2", unpack=True)

    # Find the peaks and the pixels associated
    peak_pixels, peak_heights = find_peaks(intensity, height=1700)

    # Plotting the data as well as the points where the peaks occur
    pf.line(True, pixel, intensity, "Cross Section of Pixel vs. Intensity of Supernova Cloud", "Pixel", "Intensity",
            "Intensity of Pixel", None, "black")
    pf.scatter(False, peak_pixels, peak_heights['peak_heights'], color="orange", style='o', label="Peaks")
    plt.legend()

    # Write the wavelengths associated with the pixel
    wavelengths = np.array([16079.753, 16128.608, 16194.615, 16235.376, 16317.161, 16350.650, 16388.492, 16442.155,
                            16475.648, 16512.365, 16513, 16525, 16553.814])

    # Perform linear regression on pixel s wavelength data
    popt, pcov = curve_fit(pf.line_eq, peak_pixels, wavelengths)
    m, c = popt[0], popt[1]
    u_m, u_c = np.sqrt(np.diag(pcov))

    print(f"The wavelength equation is y = {popt[0]:.5}x + {popt[1]:.5}")
    print(f"The uncertainties of the slope is {u_m:.5} and the uncertainty in the intercept is {u_c:.5}")

    # Calculating Residuals
    # Find the predicted values of the max_pixels and then gather the residuals
    predicted_wavelengths = [pf.line_eq(p, m, c) for p in peak_pixels]
    residual_wavelengths = predicted_wavelengths - wavelengths

    # Plot the residuals to see the variation in the actual values and the predicted values from the line equation
    pf.scatter(True, np.arange(0, len(residual_wavelengths)), residual_wavelengths, "Residuals of Wavelength Data",
               "Index", "Difference from Actual Value", "Residuals of Gathered Data")
    plt.hlines(0, xmin=-1, xmax=14, colors="red")
    plt.legend()
    plt.show()

    """
    --------------------------------------------------------------------------------------------------------------------
    # Reading the data from the Fe II topic to determine the wavelength
    # Plotting the data from the Fe line
    # Finding the peak in the plot which indicates the gas light emission and finding the correlated pixel
    # Calculating the velocity of the gas by using the doppler shift equation for light waves
    --------------------------------------------------------------------------------------------------------------------
    """

    fe_pixel, fe_intensity = np.loadtxt("data/TelluricOH_Data/Fe_Data", unpack=True)  # COULD TOTALLY TRY FOR BETTER

    # Plotting the data
    pf.line(True, fe_pixel, fe_intensity, "Pixel vs. Intensity For Cross Section of Supernova Gas Spectra", "Pixel",
            "Intensity", "Data", color="black")
    plt.legend()
    plt.show()

    # Finding the pixel connected to the max intensity
    max_intensity = max(fe_intensity)
    max_elem = np.where(fe_intensity == max_intensity)[0][0]

    # Gives the wavelength in μm
    fe_wave = pf.line_eq(fe_pixel[max_elem], m, c) * 1e-4
    print(fe_wave)

    # Constants needed for velocity calculation
    acc_wavelength = 1.6439981      # Units of μm
    light_speed = 3.0e8             # Units on m/s

    # Calculation of velocity of gas
    speed = (fe_wave / acc_wavelength - 1) * light_speed
    print(f"The speed of the Fe (II) gas cloud was found to be {speed:.4} m/s")
