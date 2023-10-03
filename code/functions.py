import numpy as np


def readfile(filename: str):
    """
    Reads file information using numpy.genfromtxt()
    :param filename: Name of file being read
    :return: Array of data from file
    """

    return np.genfromtxt(f"data/{filename}")
