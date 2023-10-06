import plotting_functions as pf
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as sc


if __name__ == '__main__':

    """
    ----------------------------------------------------------------------------------------------------------------    
    # Opening the two files to write their data into numpy arrays:
    
    # The first file is the spectrum of a blackbody source (Group_X_BB.dat)
    # The second file is the Neon calibration (Ne_calib.dat) which contains the 
      spectrum of a Neon lamp using a spectrograph
      --------------------------------------------------------------------------------------------------------------    
    """

    bb_data = pf.readfile("data/Group_X_BB.dat") * 1e-9
    ne_data = pf.readfile("data/Ne_calib.dat") * 1e-9
    ne_wave_data = pf.readfile("data/Ne_wavelengths") * 1e-9

    """ 
    ----------------------------------------------------------------------------------------------------------------
    # Plotting the two data sets of Pixel vs. Intensity plots:
    
    # Set up pixels to plot across the x-axis to then relate to the intensity gathered 
    # Take the length of the data set to find correlating pixels to each value
    ----------------------------------------------------------------------------------------------------------------
    """
    bb_pixel = np.arange(0, len(bb_data))
    ne_pixel = np.arange(0, len(ne_data))

    pf.scatter(True, bb_pixel, bb_data, "Black Body Spectrum", "Pixels", "Intensity", "Black Body", ".")
    pf.show()

    pf.line(True, ne_pixel, ne_data, "Neon Lamp Calibration Spectrum", "Pixels", "Relative Intensity", "TODO", None)

    """
    ----------------------------------------------------------------------------------------------------------------
    # Finding the Centroids:
    
    # Set a limit for minimum required intensity (min_intensity) to be detected by centroid calculation
    # Create a list to store values of pixels and their values reaching above the min_intensity 
    # Run a loop to compare data values to find where peaks fall
    # Finding values where the values to the left and right are smaller and higher than the min_intensity
    # Plot vertical lines on the graph to pin point centroids
    ----------------------------------------------------------------------------------------------------------------
    """

    min_intensity = 2e-7

    peak_pixels, peak_values = pf.centroids(min_intensity, ne_pixel, ne_data)

    plt.vlines(peak_pixels, ymin=min_intensity, ymax=1500e-9, color='red', linewidth=0.8, ls='-.')
    pf.show()

    """
    ----------------------------------------------------------------------------------------------------------------
    # Linear Square Fitting of Data:
    
    # Using curvefit() to take the centroids and their correlated pixel and find a linear fit for it
    # This correlates the pixel value and the approximate wavelength that follows it
    ----------------------------------------------------------------------------------------------------------------
    """

    popt, pcov = sc.curve_fit(pf.line_eq, peak_pixels, peak_values)
    m, b = popt[0], popt[1]

    m_err, b_err = np.sqrt(np.diag(pcov))

    predicted = [pf.line_eq(p, m, b) for p in peak_pixels]

    pf.line(True, peak_pixels, predicted, "Wavelength vs. Pixels", 'Pixels', r'Wavelengths ($\lambda$)'
            , 'Wavelength Solution', None)
    plt.text(40, 450e-9, '$y=mx+b$')
    plt.text(40, 435e-9, f"$m = {m:.4} \u00B1 {m_err:.4}$")
    plt.text(40, 420e-9, f"$b = {b:.4} \u00B1 {b_err:.4}$")
    pf.show()

    """
    ----------------------------------------------------------------------------------------------------------------
    # Applying the Wavelength Solution:
    
    # Using the line_eq() the wavelength at the peak can be found 
    # Using the ___WHICH EQUAITON?___ find the temperature
    ----------------------------------------------------------------------------------------------------------------
    """

    prop_const = 2.987e-3

    wave_pixels, wave_max = pf.centroids(5.1e-6, bb_pixel, bb_data)
    print(wave_max)

    T_list = []

    for wave in wave_max:
        T_list.append(prop_const / wave)
    print(T_list)
