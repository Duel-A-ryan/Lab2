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
    :param c: Y-intercept value
    :return: y = mx + b
    """
    return m * x + b


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


def m_d_unc(values, unc_values, actual):
    """
    Returns uncertainty in multiplication of values with individual uncertainties

    :param values:
    :param unc_values:
    :return: Uncertainty in final value
    """
    total = 0.0

    for i in range(len(values)):
        total += (unc_values[i] / values[i]) ** 2

    return np.sqrt(total) * actual


def a_s_unc(unc_values):
    """

    :param unc_values: list of uncertainties
    :return: Sum of uncertainties
    """
    return np.sum(unc_values)
