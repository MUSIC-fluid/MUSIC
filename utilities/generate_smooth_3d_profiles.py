#! /usr/bin/env python

import sys, shutil
from numpy import *
from os import path, makedirs
from glob import glob
import subprocess

class color:
    """
    define colors in the terminal
    """
    purple = '\033[95m'
    cyan = '\033[96m'
    darkcyan = '\033[36m'
    blue = '\033[94m'
    green = '\033[92m'
    yellow = '\033[93m'
    red = '\033[91m'
    bold = '\033[1m'
    underline = '\033[4m'
    end = '\033[0m'

grid_nrap = 11

rand_flag = 1

# the width of the Gaussian in the transverse plane
sigma_perp = 0.5

# peak position in the longitudinal direction
eta_0 = 2.0
# the width of the Gaussian in the longitudinal direction
sigma_eta = 0.5

centrality_list = ['0-5', '5-10', '10-20', '20-30',
                   '30-40', '40-50', '50-60', '60-70', '70-80']

def get_eta_factor_left(eta_local, beam_rapidity):
    exp_factor = 1.0
    if abs(eta_local) > eta_0:
        exp_factor = exp(-(abs(eta_local) - eta_0)**2./(2.*sigma_eta**2.))
    eta_left = 0.5*(1. - eta_local/beam_rapidity)*exp_factor
    return(eta_left)

def get_eta_factor_right(eta_local, beam_rapidity):
    exp_factor = 1.0
    if abs(eta_local) > eta_0:
        exp_factor = exp(-(abs(eta_local) - eta_0)**2./(2.*sigma_eta**2.))
    eta_right = 0.5*(1. + eta_local/beam_rapidity)*exp_factor
    return(eta_right)

def generate_3d_profile(data_path, ecm):
    beam_rapidity = arctanh(sqrt(1. - 1./((ecm/2.)**2.)))
    grid_eta = linspace(-beam_rapidity, beam_rapidity, grid_nrap)
    for icen, cen_string in enumerate(centrality_list):
        TA = loadtxt(
            path.join(data_path, 'nuclear_thickness_TA_fromSd_order_2_C%s.dat'
                                  % cen_string))
        TB = loadtxt(
            path.join(data_path, 'nuclear_thickness_TB_fromSd_order_2_C%s.dat'
                                  % cen_string))
        entropy_density = []
        for ieta in range(len(grid_eta)):
            eta_local = grid_eta[ieta]
            temp_density = (
                get_eta_factor_left(eta_local, beam_rapidity)*TA
                + get_eta_factor_right(eta_local, beam_rapidity)*TB)
            entropy_density.append(temp_density)
         
        with file('sdAvg_order_2_C%s_block_3d.dat' % cen_string, 'w') as outfile:
            for slice_2d in entropy_density:
                savetxt(outfile, slice_2d)

def print_help_message():
    print "Usage : "
    print(color.bold
          + "./generate_smooth_3d_profiles.py folder ecm "
          + color.end)
    print "Usage of generate_smooth_3d_profiles.py command line arguments: "
    print(color.bold + "-folder" + color.end + "  folder path")
    print(color.bold + "-ecm" + color.end
          + "     collision energy (GeV): "
          + color.purple + "7.7, 11.5, 19.6, 27, 39, 62.4, 200, 2760, 5500"
          + color.end)

if __name__ == "__main__":
    try:
        data_path = path.abspath(str(sys.argv[1]))
        ecm = float(sys.argv[2])
    except(IndexError):
        print_help_message()
        exit(0)
    generate_3d_profile(data_path, ecm)
