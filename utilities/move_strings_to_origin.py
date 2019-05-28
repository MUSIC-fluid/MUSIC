#!/usr/bin/env python
"""This script move all the strings to origin"""

import sys
import numpy as np

def print_usage():
    """This function prints out help message"""
    print("\U0000269B  Usage: {} strings_file".format(sys.argv[0]))

def main(filepath):
    """This is the main funciton"""
    data = np.loadtxt(filepath)
    data[:, 3] = 0.0
    data[:, 4] = 0.0
    np.savetxt(filepath, data, fmt="%.8e", delimiter="  ")


if __name__ == "__main__":
    try:
        FILENAME = str(sys.argv[1])
    except IndexError:
        print_usage()
        exit(0)

    main(FILENAME)
