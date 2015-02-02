#!/usr/bin/env python

import sys
from music_parameters_dict import *

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

initial_condition_dict.update({
    'Initial_profile': 11,          # type of initial condition
    'initialize_with_entropy': 0,   # 0: with energy density
    'Initial_Distribution_Filename': 'initial/edAvg_order_2_C0-5.dat',
    'Initial_Rhob_Distribution_Filename':
        'initial/rhob_fromEd_order_2_C0-5.dat',
    's_factor': 10.515679450,   # normalization factor read in initial data file

    #  envelope function in eta_s direction
    'Eta_plateau_size': 20.,          # size of the plateau in eta_s direction
    'Eta_fall_off': 0.7,              # the scale of the fall off of the plateau in eta_s direction
    'initial_eta_rhob_profile': 1,    # type of the envelope profile for rho_B's eta_s distribution
    'eta_rhob_0': 3.0,                # peak position of rho_B in eta_s direction
    'eta_rhob_width': 1.0,            # the width of the Gaussian (for initial_eta_rhob_profile == 1)
    'eta_rhob_plateau_height': 1.0,   # the relative height of the central plateau (for initial_eta_rhob_profile == 2)
    'eta_rhob_width_1': 1.0,          # the width of the Gaussian for the outside tail (for initial_eta_rhob_profile == 2)
    'eta_rhob_width_2': 0.5,          # the width of the Gaussian for the inside (for initial_eta_rhob_profile == 2)
})

hydro_dict.update({
    # grid information
    'Initial_time_tau_0': 0.6,   # starting time of the hydrodynamic evolution (fm/c)
    'Delta_Tau': 0.04,           # time step to use in the evolution [fm/c]

    'Eta_grid_size': 14.0,       # spatial rapidity range
    'Grid_size_in_eta': 4,       # number of the grid points in spatial rapidity direction
    'X_grid_size_in_fm': 26.0,   # spatial range along x direction in the transverse plane
    'Y_grid_size_in_fm': 26.0,   # spatial range along y direction in the transverse plane
    'Grid_size_in_y': 260,       # number of the grid points in y direction
    'Grid_size_in_x': 260,       # number of the grid points in x direction

    'EOS_to_use': 10,        # type of the equation of state
    'reconst_type': 1,       # the type of quantity that will be first reconstruct from T^0\mu and J^0
    'boost_invariant': 0,    # initial condition is boost invariant

    #viscosity and diffusion options
    'Viscosity_Flag_Yes_1_No_0': 1,         # turn on viscosity in the evolution
    'Include_Shear_Visc_Yes_1_No_0': 1,     # include shear viscous effect
    'Shear_to_S_ratio': 0.08,               # value of \eta/s
    'Include_Bulk_Visc_Yes_1_No_0': 0,      # include bulk viscous effect
    'Bulk_to_S_ratio': 0.1,                 # value of \zeta/s
    'Include_Rhob_Yes_1_No_0': 0,           # turn on propagation of baryon current
    'turn_on_baryon_diffusion': 0,          # turn on baryon current diffusion
    'Bulk_relaxation_time_tau_b_pi': 0.6,   # bulk relation time
    'Shear_relaxation_time_tau_pi': 0.01,   # shear relation time

    'output_evolution_data': 0,             # flag to output evolution history to file
})

freeze_out_dict.update({
    'Do_FreezeOut_Yes_1_No_0': 1,   # flag to find freeze-out surface
    'freeze_out_method': 4,         # method for hyper-surface finder

    'average_surface_over_this_many_time_steps': 5,   # the step skipped in the tau direction
    'Ncell_skip_x': 5,              # the step skipped in x direction
    'Ncell_skip_y': 5,              # the step skipped in y direction
    'N_freeze_out': 5,              # number of freeze-out surfaces
    'eps_freeze_max': 0.508,        # the maximum freeze-out energy density (GeV/fm^3)
    'eps_freeze_min': 0.100,        # the minimum freeze-out energy density (GeV/fm^3)

    'Include_deltaf': 1,        # flag to include delta f correction in Cooper-Frye formula
    'Include_deltaf_qmu': 0,    # flag to include delta f for qmu
    'Inlucde_deltaf_bulk': 0,   # flag to include delta f for bulk viscosity

    'number_of_particles_to_include': 320,  # number of thermal particles to compute for particle spectra and vn
    'particle_spectrum_to_compute': 0,      # 0: Do all up to number_of_particles_to_include

    'pseudofreeze': 1,          # calculated particle spectra in equally-spaced pseudorapidity
    'max_pseudorapidity': 2.5,  # particle spectra calculated from (0, max_pseudorapidity)
    'pseudo_steps': 11,         # number of lattice points along pseudo-rapidity
    'phi_steps': 40,            # number of points calculated in phi for Cooper-Frye
    'min_pt': 0.01,             # the minimum value of pT calculated in the Cooper-Frye
    'max_pt': 3.0,              # the maximum value of pT calculated in the Cooper-Frye
    'pt_steps': 15,             # number of the pT points calculated in the Cooper-Frye
})

collect_dict.update({
    'dNdy_y_min': -0.01,
    'dNdy_y_max': 0.01,
    'dNdy_eta_min': -2.5,
    'dNdy_eta_max': 2.5,
    'dNdy_nrap': 30,
    'dNdyptdpt_y_min': -0.01,
    'dNdyptdpt_y_max': 0.01,
    'dNdyptdpt_eta_min': -0.5,
    'dNdyptdpt_eta_max': 0.5,
})


def generate_music_input_file():
    # mode 2:
    control_dict.update({'mode': 2})
    f = open('music_input_2', 'w')
    dict_list = [control_dict, initial_condition_dict, hydro_dict,
                 freeze_out_dict, collect_dict]
    for idict in range(len(dict_list)):
        temp_list = dict_list[idict].items()
        for i in range(len(temp_list)):
            f.write('%s  %s \n' % (temp_list[i][0], str(temp_list[i][1])))
    f.close()

    # mode 3:
    control_dict.update({'mode': 3})
    f = open('music_input_3', 'w')
    dict_list = [control_dict, initial_condition_dict, hydro_dict,
                 freeze_out_dict, collect_dict]
    for idict in range(len(dict_list)):
        temp_list = dict_list[idict].items()
        for i in range(len(temp_list)):
            f.write('%s  %s \n' % (temp_list[i][0], str(temp_list[i][1])))
    f.close()

    # mode 4:
    control_dict.update({'mode': 4})
    f = open('music_input_4', 'w')
    dict_list = [control_dict, initial_condition_dict, hydro_dict,
                 freeze_out_dict, collect_dict]
    for idict in range(len(dict_list)):
        temp_list = dict_list[idict].items()
        for i in range(len(temp_list)):
            f.write('%s  %s \n' % (temp_list[i][0], str(temp_list[i][1])))
    f.close()

    # mode 13:
    control_dict.update({'mode': 13})
    f = open('music_input_13', 'w')
    dict_list = [control_dict, initial_condition_dict, hydro_dict,
                 freeze_out_dict, collect_dict]
    for idict in range(len(dict_list)):
        temp_list = dict_list[idict].items()
        for i in range(len(temp_list)):
            f.write('%s  %s \n' % (temp_list[i][0], str(temp_list[i][1])))
    f.close()

    # mode 14:
    control_dict.update({'mode': 14})
    f = open('music_input_14', 'w')
    dict_list = [control_dict, initial_condition_dict, hydro_dict,
                 freeze_out_dict, collect_dict]
    for idict in range(len(dict_list)):
        temp_list = dict_list[idict].items()
        for i in range(len(temp_list)):
            f.write('%s  %s \n' % (temp_list[i][0], str(temp_list[i][1])))
    f.close()



def print_help_message():
    print "Usage : "
    print(color.bold
          + "./generate_music_inputfile.py"
          + "[-cen centrality -shear_vis eta_over_s -EOS eos_type]"
          + color.end)
    print "Usage of runHydro.py command line arguments: "
    print(color.bold + "-shear_vis" + color.end
          + "   the specific shear viscosity used in the hydro simulation\n"
          + color.bold + "       eta/s = 0.08 [default]" + color.end)
    print(color.bold + "-EOS" + color.end
          + "   the equation of state for hydrodynamic simulation \n"
          + color.purple + color.bold + "       10 [default]"
          + color.end
          + color.purple + ", 2-5" + color.end)
    print(color.bold + "-cen" + color.end
          + "  specify the centrality bin: "
          + color.bold + "0-5 [default]" + color.end
          + color.purple + ', e.g. 20-30' + color.end)
    print(color.bold + "-h | -help" + color.end + "    This message")


if __name__ == "__main__":
    while len(sys.argv) > 1:
        option = sys.argv[1]
        del sys.argv[1]
        if option == '-cen':
            centrality = str(sys.argv[1])
            del sys.argv[1]
            initial_condition_dict.update({
                'Initial_Distribution_Filename':
                    'initial/edAvg_order_2_C%s.dat' % centrality,
                'Initial_Rhob_Distribution_Filename':
                    'initial/rhob_fromEd_order_2_C%s.dat' % centrality,})
        elif option == '-shear_vis':
            vis = float(sys.argv[1])
            del sys.argv[1]
            hydro_dict.update({'Shear_to_S_ratio': vis,})
        elif option == '-EOS':
            eos_type = int(sys.argv[1])
            del sys.argv[1]
            hydro_dict.update({'EOS_to_use': 10,})
        elif option == '-h':
            print_help_message()
            sys.exit(0)
        else:
            print sys.argv[0], ': invalid option', option
            print_help_message()
            sys.exit(1)

    generate_music_input_file()