import numpy as np
import matplotlib.pyplot as plt


def readfile(filename: str):
    """
    Reads file information using numpy.genfromtxt()
    :param filename: Name of file being read
    :return: Array of data from file
    """

    return np.genfromtxt(filename)


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


def scatter(new, data1, data2, title, x_name, y_name, label, style):
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
    """
    if new is True:
        plt.figure()
        plot_info(title, x_name, y_name)

    plt.scatter(data1, data2, label=label, marker=style)


def line(new, data1, data2, title, x_name, y_name, label, style):
    """
    Plots a line plot from given information

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

    plt.plot(data1, data2, label=label, marker=style)






