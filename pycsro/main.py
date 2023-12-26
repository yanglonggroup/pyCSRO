import os
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import read, write
import copy
from operator import *
from ase.build import sort
import time
from pycsro.io import *
from pycsro.helper import *
from pycsro.neighbor import *
from pycsro.sro import *


################################################################################################
# main

def run_pycsro_pmsro(ion1, cutoff1, filename, cutoff2=None, savename=None, skip_distance=0.1,
                     plotsave=None, cal_same_pair=None, safe_mode=None, partial_neighbors=None):
    """
    Calculate the PM-SRO of element pairs.

    Args:
        ion1: The selected elements for the PM-SRO calculation.
        cutoff1: The cutoff of thr 1st shell.
        cutoff2: The cutoff of the 2nd shell.
        filename: The absolute path of input file.
        savename: The absolute path of saved file.
        skip_distance: Skip the neighbor distance under 0.1 (Default).
        cal_same_pair: Whether calculate the wcp of same elements but different center atoms.
        safe_mode: Whether use the supercell selection function, which can reduce the calculation time.
        plotsave: Whether save the neighbor plot.
        partial_neighbors: Whether plot the partial neighbor distribution of atoms in the cell.
    """
    ion1, cutoff1, cutoff2, skip_distance, plotsave, cal_same_pair, safe_mode, dual_cutoff, single_ele = \
        pycsro.pre.default_settings(ion1, cutoff1, cutoff2, skip_distance, plotsave, cal_same_pair, safe_mode)
    pycsro.pycsroio.print_settings(ion1, cutoff1, cutoff2, filename, savename,
                                   skip_distance, dual_cutoff, cal_same_pair, safe_mode)
    # read file part
    atoms = pycsro.pycsroio.readfile(filename)  # read file
    cell, pos, ele = pycsro.pycsroio.get_config(atoms)
    axis = pycsro.pre.get_axis(cell)
    supercell = pycsro.pre.periodicity(cutoff2, axis, atoms)  # make supercell
    cell_su, pos_su, ele_su = pycsroio.get_config(supercell)
    pos_new = pycsro.pre.move_original_positions(pos, cell, cell_su)  # move original cell
    plane_n_su = pycsro.pre.cal_plane_normal_vector(cell_su)
    if safe_mode:
        pos_su_new = pos_su
        ele_su_new = ele_su
    else:
        pos_su_new, ele_su_new = \
            pycsro.pre.select_supercell(pos_su, ele_su, cutoff2, cell, cell_su, plane_n_su)  # select supercell
    # calculate neighbors
    neighbors_1, neighbors_ele, neighbors_2, neighbors_ele_2 = \
        pycsro.neighbor.cal_neighbors(ele, pos_new, ele_su_new, cutoff1, cutoff2, pos_su_new, skip_distance, ion1)
    if partial_neighbors == None or partial_neighbors.lower() == 'no' or partial_neighbors.lower() == 'n':
        save_plot_data = pycsro.neighbor.plot_neighbors(cutoff2, neighbors_1, neighbors_2)
    elif partial_neighbors.lower() == 'yes' or partial_neighbors.lower() == 'y':
        save_plot_data = pycsro.neighbor.plot_partial_neighbors(cutoff2, neighbors_1, neighbors_2, ele,
                                                                neighbors_ele, neighbors_ele_2, ion1)
    if plotsave:
        plt.savefig('neighbors', dpi=300)
        plt.show()
    else:
        plt.show()
    # calculate PM-SRO parameter
    wc_list = [f'The PM-SRO parameter for {ion1} element group in the 1st shell']
    wc_list_2 = ['+-----------------------------------------------------------------------------',
                 f'The PM-SRO parameter for {ion1} element group in the 2nd shell']
    ele_list = copy.deepcopy(ion1)
    ele_list_temp = copy.deepcopy(ele_list)
    ele_list_temp_2 = copy.deepcopy(ele_list)
    print(f'| The PM-SRO parameter for {ion1} element group in the 1st shell')
    if single_ele:
        center_ele = ion1[0]
        second_ele = ion1[0]
        average_density = 1
        neighbor_density = pycsro.sro.cal_neighbor_density(ele, neighbors_ele, center_ele, second_ele)
        wc_list_temp = pycsro.sro.cal_wc_sro(neighbor_density, average_density, center_ele, second_ele)
        wc_list.append(f'{center_ele}-{second_ele} {wc_list_temp}')
        print(f'+-----------------------------------------------------------------------------')
        if dual_cutoff:
            print(f'| The PM-SRO parameter for {ion1} element group in the 2nd shell')
            neighbor_density = pycsro.sro.cal_neighbor_density(ele, neighbors_ele_2, center_ele, second_ele)
            wc_list_temp = pycsro.sro.cal_wc_sro(neighbor_density, average_density, center_ele, second_ele)
            wc_list_2.append(f'{center_ele}-{second_ele} {wc_list_temp}')
            print(f'+-----------------------------------------------------------------------------')
    else:
        for i in ele_list:  # calculate pm-sro parameter for 1st shell
            center_ele = i
            for j in ele_list_temp:
                second_ele = j
                average_density = pycsro.sro.cal_average_density(ele, second_ele, ion1)
                neighbor_density = pycsro.sro.cal_neighbor_density(ele, neighbors_ele, center_ele, second_ele)
                wc_list_temp = pycsro.sro.cal_pm_sro(neighbor_density, average_density, center_ele, second_ele)
                wc_list.append(f'{center_ele}-{second_ele} {wc_list_temp}')
            if cal_same_pair:
                del ele_list_temp[0]
        print(f'+-----------------------------------------------------------------------------')
        if dual_cutoff:  # check dual_cutoff
            print(f'| The PM-SRO parameter for {ion1} element group in the 2nd shell')
            for i in ele_list:  # calculate pm-sro parameter for 2nd shell
                center_ele = i
                for j in ele_list_temp_2:
                    second_ele = j
                    average_density = pycsro.sro.cal_average_density(ele, second_ele, ion1)
                    neighbor_density = pycsro.sro.cal_neighbor_density(ele, neighbors_ele_2, center_ele, second_ele)
                    wc_list_temp = pycsro.sro.cal_pm_sro(neighbor_density, average_density, center_ele, second_ele)
                    wc_list_2.append(f'{center_ele}-{second_ele} {wc_list_temp}')
                if cal_same_pair:
                    del ele_list_temp_2[0]
            print(f'+-----------------------------------------------------------------------------')
    if savename:  # savefile
        pycsro.pycsroio.save_wcp_2(os.path.abspath(savename), wc_list)
        if dual_cutoff:
            pycsro.pycsroio.save_wcp_2(os.path.abspath(savename), wc_list_2)
        end = ['+-----------------------------------------------------------------------------']
        pycsro.pycsroio.save_wcp_2(os.path.abspath(savename), end)
        pycsro.pycsroio.save_plot(savename, save_plot_data)
    return ion1, cutoff1, filename