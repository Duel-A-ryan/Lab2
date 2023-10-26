import plotting_functions as pf
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as sc
from scipy.signal import find_peaks


if __name__ == '__main__':

    """
    ----------------------------------------------------------------------------------------------------------------    
    # Opening the two files to write their data into numpy arrays:
    
    # The first file is the spectrum of a blackbody source (Group_X_BB.dat)
    # The second file is the Neon calibration (Ne_calib.dat) which contains the 
      spectrum of a Neon lamp using a spectrograph
      --------------------------------------------------------------------------------------------------------------    
    """

    # Opening the .dat files
    bb_data = pf.readfile("data/Group_X_BB.dat")
    ne_data = pf.readfile("data/Ne_calib.dat")
    ne_wave_data = pf.readfile("data/Ne_wavelengths")

    """ 
    ----------------------------------------------------------------------------------------------------------------
    # Plotting the two data sets of Pixel vs. Intensity plots:
    
    # Set up pixels_of_peaks to plot across the x-axis to then relate to the intensity gathered 
    # Take the length of the data set to find correlating pixels_of_peaks to each value
    ----------------------------------------------------------------------------------------------------------------
    """
    # Making arrays to associate to the pixels_of_peaks of the spectrometer
    pixels = np.arange(0, len(bb_data))

    # Plotting the data of the Black Body
    pf.scatter(True, pixels, bb_data, "Black Body Spectrum", "Pixels", "Intensity", "Black Body", ".", color='black')

    # Plotting the data of the Neon Lamp
    pf.line(True, pixels, ne_data, "Neon Lamp Calibration Spectrum", "Pixels", "Relative Intensity", "TODO", None,
            "black")

    """
    ----------------------------------------------------------------------------------------------------------------
    # Finding the Centroids:
    
    # Use the find_peaks function from scipy.signal to find the centroids
    # Use a height of 22 setting to find all 33 points
    # Manually correlate each peak with its wavelengths
    # Plot the dots on the Neon lamp plot 
    ----------------------------------------------------------------------------------------------------------------
    """
    # Finding peak pixels_of_peaks and the heights
    peak_pixels, peak_heights = find_peaks(ne_data, height=22)

    # Plot the Neon Lamp data again with the peak intensities labelled
    pf.line(True, pixels, ne_data, "Neon Lamp Calibration Spectrum", "Pixels", "Relative Intensity", "TODO", None,
            "black")
    pf.scatter(False, peak_pixels, peak_heights['peak_heights'], label="Peak Intensities", style='.', color='red')
    plt.show()

    """
    ----------------------------------------------------------------------------------------------------------------
    # Linear Square Fitting of Data:
    
    # Create two new arrays filled with the pixels_of_peaks and wavelength values observed in the plots
    # Using curvefit() to take the centroids and their correlated pixel and find a linear fit for it and the expected 
      wavelengths
    # This correlates the pixel value and the approximate wavelength that follows it
    ----------------------------------------------------------------------------------------------------------------
    """

    pixels_of_peaks = np.array([53, 205, 227, 241, 258, 284, 292, 310, 339, 374, 394, 418, 447, 464, 475, 522, 533, 560, 595,
                                611, 698, 742, 804, 830, 916, 928, 938, 954, 963])
    waves = np.array([540.056, 576.441, 582.015, 585.249, 588.189, 594.483, 597.553, 602.000, 607.433, 614.306, 616.359,
                      621.728, 626.649, 638.299, 640.255, 650.653, 653.288, 659.895, 667.828, 671.704, 692.947, 703.241,
                      717.394, 724.512, 743.890, 747.224, 748.887, 753.377, 754.404]) * 1e-9

    popt, pcov = sc.curve_fit(pf.line_eq, pixels_of_peaks, waves)
    m, b = popt[0], popt[1]
    m_err, b_err = np.sqrt(np.diag(pcov))

    predicted = [pf.line_eq(p, m, b) for p in peak_pixels]
    unc_predicted = [pf.m_d_unc([m, p], [m_err, (np.sqrt(2))/(2)], m*p) for p in pixels_of_peaks] + b_err

    pf.line(True, peak_pixels, predicted, "Wavelength vs. Pixels", 'Pixels', r'Wavelengths ($m$)',
            'Wavelength Solution', None, "black")
    plt.text(40, 720e-9, '$y=mx+b$')
    plt.text(40, 705e-9, f"$m = {m:.4} \u00B1 {m_err:.4}$")
    plt.text(40, 690e-9, f"$b = {b:.4} \u00B1 {b_err:.4}$")
    pf.show()

    """
    #----------------------------------------------------------------------------------------------------------------
    # Applying the Wavelength Solution:
    
    # Using the line_eq() the wavelength at the peak can be found 
    # Using Wien's Law to find the temperature in Kelvin
    #----------------------------------------------------------------------------------------------------------------
    """

    prop_const = 2.987e-3

    max_bb = max(bb_data)

    max_elem = np.where(bb_data == max_bb)[0][0]
    max_wave = pf.line_eq(pixels[max_elem], m, b)

    unc_max_wave = pf.m_d_unc([m, max_elem], [m_err, (np.sqrt(2))/(2)], max_elem*m) + b_err

    temperature = prop_const / max_wave
    unc_temperature = pf.m_d_unc([max_wave], [unc_max_wave], temperature)
    print(f"The temperature found for the blackbody data is {temperature:.5} \u00B1 {unc_temperature:.3} Kelvin")

