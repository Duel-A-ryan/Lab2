from scipy.stats import linregress
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
    pixel, intensity = np.loadtxt("data/New folder/TelluricOH_2", unpack=True)

    # Find the peaks and the pixels associated
    peak_pixels, peak_heights = find_peaks(intensity, height=1700)

    # Plotting the data as well as the points where the peaks occur
    pf.line(True, pixel, intensity, "Cross Section of Pixel vs. Intensity of Supernova Cloud", "Pixel", "Intensity",
            "Intensity of Pixel", None, "black")
    pf.scatter(False, peak_pixels, peak_heights['peak_heights'], color="orange", style='o', label="Peaks")
    pf.show("Plots/Intensity Plots/OH Lines Data")

    # Write the wavelengths associated with the pixel
    wavelengths = np.array([16079.753, 16128.608, 16194.615, 16235.376, 16317.161, 16350.650, 16388.492, 16442.155,
                            16475.648, 16512.365, 16513, 16525, 16553.814]) / 10000

    # Perform linear regression on pixel s wavelength data
    popt_linear, pcov_linear = curve_fit(pf.line_eq, peak_pixels, wavelengths)
    slope, intercept = popt_linear
    u_slope, u_intercept = np.sqrt(np.diag(pcov_linear))

    # Quadratic Fit of Data
    popt_quad, pcov_quad = curve_fit(pf.quad_eq, peak_pixels, wavelengths)
    a_quad, b_quad, c_quad = popt_quad
    u_a_quad, u_b_quad, u_c_quad = np.sqrt(np.diag(pcov_quad))

    # Cubic Fit of Data
    popt_cube, pcov_cube = curve_fit(pf.cubic_eq, peak_pixels, wavelengths)
    a_cube, b_cube, c_cube, d_cube = popt_cube
    u_a_cube, u_b_cube, u_c_cube, u_d_cube = np.sqrt(np.diag(pcov_cube))

    # Plotting Data
    pf.scatter(True, peak_pixels, wavelengths, "Wavelength vs. Pixels with 3 Polynomial Fits", 'Pixels',
               r'Wavelengths ($m$)', 'Linear')

    # Linear Data Plotting and R Squared Calculation
    predicted_line = [pf.line_eq(p, slope, intercept) for p in peak_pixels]
    pf.line(False, peak_pixels, predicted_line, label="Linear Fit", color='red')
    r2_linear = linregress(wavelengths, predicted_line)[2]

    # Quadratic Data Plotting and R Squared Calculation
    predicted_quad = [pf.quad_eq(p, a_quad, b_quad, c_quad) for p in peak_pixels]
    pf.line(False, peak_pixels, predicted_quad, label="Quadratic Fit", color="blue")
    r2_quad = linregress(wavelengths, predicted_quad)[2]

    # Cubic Data Plotting and R Squared Calculation
    predicted_cube = [pf.cubic_eq(p, a_cube, b_cube, c_cube, d_cube) for p in peak_pixels]
    pf.line(False, peak_pixels, predicted_cube, label="Cubic Fit", color="green")
    r2_cube = linregress(wavelengths, predicted_cube)[2]

    pf.show("Plots/Fits/Fitted Data")

    # Find the predicted values of the max_pixels and then gather the residuals
    p_wavelengths_line = [pf.line_eq(p, slope, intercept) for p in peak_pixels]
    p_wavelengths_quad = [pf.quad_eq(p, a_quad, b_quad, c_quad) for p in peak_pixels]
    p_wavelengths_cube = [pf.cubic_eq(p, a_cube, b_cube, c_cube, d_cube) for p in peak_pixels]

    # Finding the uncertainties in the values for all three fits
    u_wavelengths_line = [pf.a_s_unc([pf.m_d_unc(slope, u_slope, p, 0.5, 1), u_intercept]) for p in peak_pixels]

    u_wavelengths_quad = [pf.a_s_unc([pf.m_d_unc(a_quad, u_a_quad, p**2, 0.5, 2),
                                      pf.m_d_unc(b_quad, u_b_quad, p, 0.5, 1),
                                      u_c_quad]) for p in peak_pixels]

    u_wavelengths_cube = [pf.a_s_unc([pf.m_d_unc(a_cube, u_a_cube, p**3, 0.5, 3),
                                      pf.m_d_unc(b_cube, u_b_cube, p**2, 0.5, 2),
                                      pf.m_d_unc(c_cube, u_c_cube, p, 0.5, 1),
                                      d_cube]) for p in peak_pixels]

    # Calculating Residuals for all 3 Fits
    residual_line = predicted_line - wavelengths
    residual_quad = predicted_quad - wavelengths
    residual_cube = predicted_cube - wavelengths

    """-----Plot the Residuals-----"""

    # Residuals for Linear
    pf.scatter(True, peak_pixels, residual_line, "Residuals of Wavelength Data from Linear Fit",
               "Index", "Difference from Actual Value", "Residuals of Gathered Data")
    # plt.errorbar(peak_pixels, residual_line, xerr=np.ones(len(peak_pixels)) * 0.5, yerr=u_wavelengths_line)
    plt.hlines(0, xmin=-1, xmax=180, colors="red", label="Zero Line")

    pf.show("Plots/Residuals/2-2 Linear_Residuals")

    # Residuals for Quadratic
    pf.scatter(True, peak_pixels, residual_quad, "Residuals of Wavelength Data from Quadratic Fit",
               "Index", "Difference from Actual Value", "Residuals of Gathered Data")
    # plt.errorbar(peak_pixels, residual_quad, xerr=np.ones(len(peak_pixels)) * 0.5, yerr=u_wavelengths_quad)
    plt.hlines(0, xmin=-1, xmax=180, colors="red", label="Zero Line")
    #pf.show("Plots/Residuals/2-2 Quadratic_Residuals")

    # Residuals for Cubic
    pf.scatter(True, peak_pixels, residual_cube, "Residuals of Wavelength Data from Cubic Fit",
               "Index", "Difference from Actual Value", "Residuals of Gathered Data")
    # plt.errorbar(peak_pixels, residual_cube, xerr=np.ones(len(peak_pixels)) * 0.5, yerr=u_wavelengths_cube)
    plt.hlines(0, xmin=-1, xmax=180, colors="red", label="Zero Line")
    #pf.show("Plots/Residuals/2-2 Cubic_Residuals")
    plt.show()
    """
    --------------------------------------------------------------------------------------------------------------------
    # Reading the data from the Fe II topic to determine the wavelength
    # Plotting the data from the Fe line
    # Finding the peak in the plot which indicates the gas light emission and finding the correlated pixel
    # Calculating the velocity of the gas by using the doppler shift equation for light waves
    --------------------------------------------------------------------------------------------------------------------
    """
    # Reading Fe II Data
    fe_pixel, fe_intensity = np.loadtxt("data/New folder/Fe_Data", unpack=True)  # TODO: Do a vertical line

    # Plotting the data
    pf.line(True, fe_pixel, fe_intensity, "Pixel vs. Intensity For Cross Section of Supernova Gas Spectra", "Pixel",
            "Intensity", "Data", color="black")
    pf.show("Plots/Intensity Plots/FeII Data")

    # Finding the pixel connected to the max intensity
    max_intensity = max(fe_intensity)
    max_elem = np.where(fe_intensity == max_intensity)[0][0]
    max_pixel = fe_pixel[max_elem]

    # Finding wavelength and uncertainty from linear fit
    wave_line = pf.line_eq(max_pixel, slope, intercept)
    unc_wave_line = pf.a_s_unc((pf.m_d_unc(slope, u_slope, max_pixel, 0.5, 1), u_intercept))

    # Finding wavelength and uncertainty from quadratic fit
    wave_quad = pf.quad_eq(max_pixel, a_quad, b_quad, c_quad)
    unc_wave_quad = pf.a_s_unc((pf.m_d_unc(a_quad, u_a_quad, max_pixel, 0.5, 2),
                                pf.m_d_unc(b_quad, u_b_quad, max_pixel, 0.5, 1), u_c_quad))

    # Finding wavelength and uncertainty from cubic fit
    wave_cube = pf.cubic_eq(max_pixel, a_cube, b_cube, c_cube, d_cube)
    unc_wave_cube = pf.a_s_unc((pf.m_d_unc(a_cube, u_a_cube, max_pixel, 0.5, 3),
                                pf.m_d_unc(b_cube, u_b_cube, max_pixel, 0.5, 2),
                                pf.m_d_unc(c_cube, u_c_cube, max_pixel, 0.5, 1), u_d_cube))

    print(f"Linear: {wave_line} \u00B1 {unc_wave_line}")
    print(f"Quadratic: {wave_quad} \u00B1 {unc_wave_quad}")
    print(f"Cubic: {wave_cube} \u00B1 {unc_wave_cube}")

    # Constants needed for velocity calculation
    acc_wavelength = 1.6439981  # Units of μm
    light_speed = 3.0e8  # Units on m/s

    # Calculation of velocity of gas
    speed_line = (wave_line / acc_wavelength - 1) * light_speed
    speed_quad = (wave_quad / acc_wavelength - 1) * light_speed
    speed_cube = (wave_cube / acc_wavelength - 1) * light_speed

    # unc_speed_line = unc_wave_line * acc_wavelength * light_speed
    unc_speed_line = unc_wave_line / acc_wavelength * light_speed
    unc_speed_quad = unc_wave_quad / acc_wavelength * light_speed
    unc_speed_cube = unc_wave_cube / acc_wavelength * light_speed

    print(f"The Linear speed of the Fe (II) gas cloud was found to be {speed_line:.4} m/s \u00B1 {unc_speed_line:.4}")
    print(f"The Quadratic speed of the Fe (II) gas cloud was found to be {speed_quad:.4} m/s \u00B1 {unc_speed_quad:.4}")
    print(f"The Cubic speed of the Fe (II) gas cloud was found to be {speed_cube:.4} m/s \u00B1 {unc_speed_cube:.4}")
