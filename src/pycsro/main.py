from pycsro.io import *
from pycsro.helper import *
from pycsro.neighbor import *
from pycsro.sro import *
import numpy as np


################################################################################################
# main

def run_pycsro_pmsro(ion1, cutoff1, file_name, cutoff2=None, save_name=None, skip_distance=0.1,
                     plot_save=None, cal_same_pair=None, safe_mode=None, partial_neighbors=None, xyz=None):
    """
    Calculate the PM-SRO of element pairs.

    Args:
        ion1 (str): The selected elements for the PM-SRO calculation.
        cutoff1 (float): The cutoff of thr 1st shell.
        cutoff2 (float): The cutoff of the 2nd shell.
        file_name (str): The absolute path of input file.
        save_name (str): The absolute path of saved file.
        skip_distance (float): Skip the neighbor distance under 0.1 (Default).
        cal_same_pair (str): Whether calculate the wcp of same elements but different center atoms.
        safe_mode (str): Whether use the supercell selection function, which can reduce the calculation time.
        plot_save (str): Whether save the neighbor plot.
        partial_neighbors (str): Whether plot the partial neighbor distribution of atoms in the cell.
        xyz (str): Whether the format of input file is XYZ.
    """
    file_name, initial_file_name = check_readfile(file_name)
    ion1, cutoff1, cutoff2, skip_distance, plot_save, cal_same_pair, safe_mode, dual_cutoff, single_ele, \
    partial_neighbors, xyz = \
        default_settings(ion1, cutoff1, cutoff2, skip_distance, plot_save, cal_same_pair, safe_mode,
                         partial_neighbors, xyz, file_name)
    print_settings(ion1, cutoff1, cutoff2, file_name, save_name,
                   skip_distance, dual_cutoff, cal_same_pair, safe_mode, xyz, initial_file_name)
    # read file part
    if xyz:
        cell, pos_new, ele, pos_su_new, ele_su_new = readfile_xyz(file_name)  # read xyz file
    else:
        cell, pos_new, ele, cell_su, pos_su, ele_su, plane_n_su = readfile_pmsro(file_name, cutoff2)
        pos_su_new, ele_su_new = safemode(safe_mode, pos_su, ele_su, cutoff2, cell, cell_su, plane_n_su)
    # calculate neighbors
    neighbors_1, neighbors_ele, neighbors_2, neighbors_ele_2 = \
        cal_neighbors(ele, pos_new, ele_su_new, cutoff1, cutoff2, pos_su_new, skip_distance, ion1)
    if xyz:
        save_plot_data = plot_for_nieighbors(partial_neighbors, cutoff2, neighbors_1, neighbors_2, ele,
                                             neighbors_ele, neighbors_ele_2, ion1, plot_save)
    else:
        save_plot_data = plot_for_rdf(partial_neighbors, cutoff2, neighbors_1, neighbors_2,
                                      ele, neighbors_ele, neighbors_ele_2, ion1, plot_save, cell)
    wc_list, wc_list_2, ele_list, ele_list_temp, ele_list_temp_2 = cal_sro_default_settings(ion1)
    if single_ele:
        wc_list, wc_list_2 = single_sro_cal_fun(ion1, ele, neighbors_ele, neighbors_ele_2, wc_list,
                                                wc_list_2,dual_cutoff)
    else:
        wc_list, wc_list_2 = pmsro_cal_fun(wc_list, wc_list_2, ele_list, ele_list_temp, ele_list_temp_2, ele, ion1,
                                           neighbors_ele, neighbors_ele_2, cal_same_pair, dual_cutoff)
    save_file(save_name, dual_cutoff, wc_list, wc_list_2, save_plot_data, xyz)
    return ion1, cutoff1, file_name