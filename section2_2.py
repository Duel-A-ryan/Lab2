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
    pixel, intensity = np.loadtxt("data/New folder/TelluricOH_2", unpack=True)

    # Find the peaks and the pixels associated
    peak_pixels, peak_heights = find_peaks(intensity, height=1700)

    # Plotting the data as well as the points where the peaks occur
    pf.line(True, pixel, intensity, "Cross Section of Pixel vs. Intensity of Supernova Cloud", "Pixel", "Intensity",
            "Intensity of Pixel", None, "black")
    pf.scatter(False, peak_pixels, peak_heights['peak_heights'], color="orange", style='o', label="Peaks")
    plt.legend()

    # Write the wavelengths associated with the pixel
    wavelengths = np.array([16079.753, 16128.608, 16194.615, 16235.376, 16317.161, 16350.650, 16388.492, 16442.155,
                            16475.648, 16512.365, 16513, 16525, 16553.814]) / 10

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

    """
    # Plotting Data
    pf.scatter(True, peak_pixels, wavelengths)

    # Linear Data Plotting
    predicted_line = [pf.line_eq(p, slope, intercept) for p in peak_pixels]
    # unc_predicted_line = [pf.m_d_unc([slope, p], [u_slope, 0.5], slope * p) for p in peak_pixels] + u_intercept

    pf.line(False, peak_pixels, predicted_line, "Wavelength vs. Pixels", 'Pixels', r'Wavelengths ($m$)',
            'Linear', None, color='red')

    # Quadratic Data Plotting
    predicted_quad = [pf.quad_eq(p, a_quad, b_quad, c_quad) for p in peak_pixels]

    u_ap_quad = [pf.m_d_unc(a_quad, u_a_quad, p**2, 0.5, 2) for p in peak_pixels]
    u_bp_quad = [pf.m_d_unc(a_quad, u_a_quad, p, 0.5, 1) for p in peak_pixels]
    unc_predicted_quad = [pf.a_s_unc([u_ap_quad[i], u_bp_quad[i], u_c_quad]) for i in range(0, len(u_ap_quad))]


    pf.line(False, peak_pixels, predicted_quad, "Wavelength vs. Pixels", 'Pixels', r'Wavelengths ($m$)',
            label='Quadratic', style=None, color="blue")

    # Cubic Data Plotting
    predicted_cube = [pf.cubic_eq(p, a_cube, b_cube, c_cube, d_cube) for p in peak_pixels]

    u_ap_cube = [pf.m_d_unc(a_cube, u_a_cube, p**3, 0.5, 3) for p in pixel]
    u_bp_cube = [pf.m_d_unc(a_quad, u_a_quad, p ** 2, 0.5, 2) for p in pixel]
    unc_predicted_cube = [pf.m_d_unc([a_cube, p ** 3], [u_a_cube, 0.5], a_cube * p ** 3) + pf.m_d_unc([b_cube, p ** 2],
                            [u_b_cube, 0.5],
                                                                                                      b_cube * p ** 2) + pf.m_d_unc(
        [c_cube, p], [u_c_cube, 0.5], c_cube * p) for p
                          in peak_pixels] + u_c_quad

    pf.line(False, peak_pixels, predicted_cube, "Wavelength vs. Pixels", 'Pixels', r'Wavelengths ($m$)',
            label='Cubic', style=None, color="green")

    pf.show()
    """

    # Calculating Residuals

    # Find the predicted values of the max_pixels and then gather the residuals
    predicted_wavelengths = [pf.line_eq(p, slope, intercept) for p in peak_pixels]
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

    fe_pixel, fe_intensity = np.loadtxt("data/New folder/Fe_Data", unpack=True)  # COULD TOTALLY TRY FOR BETTER

    # Plotting the data
    pf.line(True, fe_pixel, fe_intensity, "Pixel vs. Intensity For Cross Section of Supernova Gas Spectra", "Pixel",
            "Intensity", "Data", color="black")
    plt.legend()
    plt.show()

    # Finding the pixel connected to the max intensity
    max_intensity = max(fe_intensity)
    max_elem = np.where(fe_intensity == max_intensity)[0][0]

    max_pixel = fe_pixel[max_elem]

    # Finding wavelength and uncertainty from linear fit
    wave_line = pf.line_eq(max_pixel, slope, intercept)
    unc_wave_line = pf.a_s_unc((pf.m_d_unc(slope, u_slope, max_pixel, 0.5, 1), u_intercept))

    wave_quad = pf.quad_eq(max_pixel, a_quad, b_quad, c_quad)
    unc_wave_quad = pf.a_s_unc((pf.m_d_unc(a_quad, u_a_quad, max_pixel, 0.5, 2),
                                pf.m_d_unc(b_quad, u_b_quad, max_pixel, 0.5, 1), u_intercept))

    wave_cube = pf.cubic_eq(max_pixel, a_cube, b_cube, c_cube, d_cube)
    unc_wave_cube = pf.a_s_unc((pf.m_d_unc(a_cube, u_a_cube, max_pixel, 0.5, 3),
                                pf.m_d_unc(b_cube, u_b_cube, max_pixel, 0.5, 2),
                                pf.m_d_unc(c_cube, u_c_cube, max_pixel, 0.5, 1), u_intercept))

    print(f"Linear: {wave_line} \u00B1 {unc_wave_line}")
    print(f"Quadratic: {wave_quad} \u00B1 {unc_wave_quad}")
    print(f"Cubic: {wave_cube} \u00B1 {unc_wave_cube}")

    # Constants needed for velocity calculation
    acc_wavelength = 1.6439981  # Units of μm
    light_speed = 3.0e8  # Units on m/s

    # Calculation of velocity of gas
    speed = (wave_line / acc_wavelength - 1) * light_speed
    print(f"The speed of the Fe (II) gas cloud was found to be {speed:.4} m/s")
