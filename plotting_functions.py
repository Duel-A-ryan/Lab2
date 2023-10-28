import numpy as np
import matplotlib.pyplot as plt


def readfile(filename: str):
    """
    Reads file information using numpy.genfromtxt()
    :param filename: Name of file being read
    :return: Array of data from file
    """

    return np.genfromtxt(filename)


def show():
    """
    Makes legend and shows plots
    """
    plt.legend()
    plt.show()


def plot_info(title, x_name, y_name):
    """
    Adds axis titles to plot

    :param title: Label of plot
    :param x_name: Label for x-axis
    :param y_name: Label for y-axis
    """

    plt.title(title)
    plt.xlabel(x_name)
    plt.ylabel(y_name)


def scatter(new, data1, data2, title="TITLE", x_name="X NAME", y_name="Y NAME", label=None, style=None,
            color="black"):
    """
    Plots a scatter plot from given information

    :param new: Whether to make a new figure to plot on
    :param data1: First array of data
    :param data2: Second array of data
    :param title: Title for plot
    :param x_name: Label for x-axis
    :param y_name: Label for y-axis
    :param label: Name of plotted data
    :param style: Style choice for scatter plots
    :param color: Colour of the points
    """
    if new is True:
        plt.figure()
        plot_info(title, x_name, y_name)

    plt.scatter(data1, data2, label=label, marker=style, color=color)


def line(new, data1, data2, title, x_name, y_name, label, style=None, color='black'):
    """
    Plots a line plot from given information

    :param color:
    :param new: Whether to make a new figure to plot on
    :param data1: First array of data
    :param data2: Second array of data
    :param title: Title for plot
    :param x_name: Label for x-axis
    :param y_name: Label for y-axis
    :param label: Name of plotted data
    :param style: Style choice for scatter plots
    """
    if new is True:
        plt.figure()
        plot_info(title, x_name, y_name)

    plt.plot(data1, data2, label=label, marker=style, color=color)


def line_eq(x, m, b):
    """

    :param x: Input values
    :param m: Slope value
    :param b: Y-intercept value
    :return: y = mx + b
    """
    return m * x + b


def quad_eq(x, a, b, c):
    """

    :param x: Input values
    :param a: Coefficient on x squared
    :param b: Coefficient on x
    :param c: Y-intercept
    :return: y = ax^2 + bx + c
    """
    return a*x**2 + b*x + c


def cubic_eq(x, a, b, c, d):
    """

    :param x: Input values
    :param a: Coefficient on x cubed
    :param b: Coefficient on x squared
    :param c: Coefficient on x
    :param d: Y-intercept
    :return: y = ax^3 + bx^2 + cx + d
    """
    return a*x**3 + b*x**2 + c*x + d


def centroids(min_intensity, ne_pixel, intensity_data):
    """
    Looks for centroids in a plot above the minimum value

    :param min_intensity: minimum intensity needed to consider looking at pixel for centroid
    :param ne_pixel: data of pixels_of_peaks
    :param intensity_data: data of intensity
    """
    peak_pixels = []
    peak_values = []

    """
    for i in range(3, len(pixels_of_peaks) - 4):
        if (intensity_data[i] > min_intensity) and (intensity_data[i] > (intensity_data[i + 1] and intensity_data[i + 2] and intensity_data[i + 3])) and \
                (intensity_data[i] < (intensity_data[i - 1] and intensity_data[i - 2] and intensity_data[i - 3])):
            peak_pixels.append(pixels_of_peaks[i])
            peak_values.append(intensity_data[i])
    """

    for i in range(3, len(ne_pixel) - 4):
        if (intensity_data[i] > min_intensity) and \
                (intensity_data[i] > (intensity_data[i + 1] and intensity_data[i + 2]))\
                and (intensity_data[i] < (intensity_data[i - 1] and intensity_data[i - 2])):

            peak_pixels.append(ne_pixel[i])
            peak_values.append(intensity_data[i])

    return peak_pixels, peak_values


def m_d_unc(a, u_a, b, u_b, p_b):

    """
    Returns uncertainty in multiplication of values with individual uncertainties

    :param values:
    :param unc_values:
    :return: Uncertainty in final value
    """

    return np.sqrt((u_a / a) ** 2 + (p_b*u_b / b) ** 2) * (a*b**p_b)


def a_s_unc(unc_values):
    """
    Returns uncertainty in multiplication of values with individual uncertainties

    :param unc_values:
    :return: Uncertainty in final value
    """
    total = 0.0

    for i in range(len(unc_values)):
        total += (unc_values[i]) ** 2

    return np.sqrt(total)


def model_fit(measured: list, predicted: list, uncertainty: list):
    """
    Precondition:
    len(measured) == len(predicted) == len(uncertainty)

    measured: list of measured values
    predicted: list of predicted values
    param uncertainty: list of uncertainties in measured values
    returns: X_r squared metric
    """
    sums = 0

    for i in range(0, len(measured)):
        sums += ((measured[i] - predicted[i]) / (uncertainty[i])) ** 2

    return sums / (len(measured) - 2)
