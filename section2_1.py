import plotting_functions as pf
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as sc
from scipy.signal import find_peaks
from scipy.stats import linregress

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
    bb_data = pf.readfile("data/part1and2/Group_X_BB.dat")
    ne_data = pf.readfile("data/part1and2/Ne_calib.dat")
    ne_wave_data = pf.readfile("data/part1and2/Ne_wavelengths")

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
    pf.show("Plots/Intensity Plots/Neon Lamp Data")

    """
    ----------------------------------------------------------------------------------------------------------------
    # Linear Square Fitting of Data:
    
    # Create two new arrays filled with the pixels_of_peaks and wavelength values observed in the plots
    # Using curvefit() to take the centroids and their correlated pixel and find a linear fit for it and the expected 
      wavelengths
    # Repeating previous step but with quadratic and cubic fits
    # Shows the wavelength correlated with each pixel through 3 different functions
    ----------------------------------------------------------------------------------------------------------------
    """

    # Arrays containing the wavelengths and the associated pixels
    pixels_of_peaks = np.array([53, 205, 227, 241, 258, 284, 292, 310, 339, 374, 394, 418, 447, 464, 475, 522, 533, 560,
                                595, 611, 698, 742, 804, 830, 916, 928, 938, 954, 963])
    waves = np.array([540.056, 576.441, 582.015, 585.249, 588.189, 594.483, 597.553, 602.000, 607.433, 614.306, 616.359,
                      621.728, 626.649, 638.299, 640.255, 650.653, 653.288, 659.895, 667.828, 671.704, 692.947, 703.241,
                      717.394, 724.512, 743.890, 747.224, 748.887, 753.377, 754.404]) * 1e-9

    # Linear Fit of Data
    popt_linear, pcov_linear = sc.curve_fit(pf.line_eq, pixels_of_peaks, waves)
    slope, intercept = popt_linear[0], popt_linear[1]
    u_slope, u_intercept = np.sqrt(np.diag(pcov_linear))

    # Linear Data Plotting
    predicted_line = [pf.line_eq(p, slope, intercept) for p in pixels_of_peaks]

    # R^2 value of Data
    r2_data = linregress(waves, predicted_line)[2]

    # Plotting the previous plot with peaks added
    pf.scatter(True, pixels_of_peaks, waves)
    pf.line(False, pixels_of_peaks, predicted_line, "Wavelength vs. Pixels", 'Pixels', r'Wavelengths ($m$)',
            'Wavelength Solution', None, color='red')

    pf.show("Plots/Fits/1-Wavelength vs Pixels")

    """-----Residuals-----"""
    u_wavelengths_line = [pf.a_s_unc([pf.m_d_unc(slope, u_slope, p, 0.5, 1), u_intercept]) for p in pixels_of_peaks]

    residual_line = predicted_line - waves

    pf.scatter(True, pixels_of_peaks, residual_line, "Residuals of Wavelength Data from Linear Fit",
               "Index", "Difference from Actual Value", "Residuals of Gathered Data")
    plt.errorbar(pixels_of_peaks, residual_line, xerr=np.ones(len(residual_line)) * 0.5, yerr=u_wavelengths_line, label="Uncertainties")
    plt.hlines(0, xmin=-1, xmax=1000, colors="red", label="Zero Line")

    pf.show("Plots/Residuals/1-Residuals Linear")
    """
    #----------------------------------------------------------------------------------------------------------------
    # Applying the Wavelength Solution:
    
    # Using the line_eq() the wavelength at the peak can be found 
    # Using Wien's Law to find the temperature in Kelvin
    #----------------------------------------------------------------------------------------------------------------
    """

    prop_const = 2.987e-3           # Units of meter Kelvins

    # Finding the Max Intensity
    max_bb = max(bb_data)
    max_elem = np.where(bb_data == max_bb)[0][0]
    max_wave = pf.line_eq(pixels[max_elem], slope, intercept)

    # Uncertainty in the observed wavelength
    unc_max_wave = pf.a_s_unc((pf.m_d_unc(slope, u_slope, max_elem, 0.5, 1), u_intercept))

    # Temperature Calculation and Uncertainty
    temperature = prop_const / max_wave
    unc_temperature = (unc_max_wave / max_wave) * temperature
    print(f"The temperature found for the blackbody data is {temperature:.5} \u00B1 {unc_temperature:.3} Kelvin")
