import numpy as np
import matplotlib.pyplot as plt
import copy
import time
from scipy.signal import savgol_filter


################################################################################################
# calculate neighbors

def cal_distance(center_atom, second_atom, pos_new, pos_su):
    """
    Calculates the distance between atoms.

    Args:
        center_atom (str): The position of center atom.
        second_atom (str): The position of second atom.
        pos_new (numpy.ndarray): The position of atoms in the cell.
        pos_su (numpy.ndarray):  The position of atoms in the supercell.

    Returns:
        distance (float): The distance between atoms.
    """
    distance = ((pos_new[center_atom][0] - pos_su[second_atom][0]) ** 2 +
                (pos_new[center_atom][1] - pos_su[second_atom][1]) ** 2 +
                (pos_new[center_atom][2] - pos_su[second_atom][2]) ** 2) ** 0.5
    return distance


def cal_neighbors(ele, pos_new, ele_su, cutoff1, cutoff2, pos_su, skip_distance, ion):
    """
    Calculate the neighbor of atoms in the cell.

    Args:
        ele (list): The elements of atoms in the cell.
        pos_new (numpy.ndarray): The position of atoms in the cell.
        ele_su (list): The elements of atoms in the supercell.
        cutoff1 (float): The cutoff of thr 1st shell.
        cutoff2 (float): The cutoff of the 2nd shell.
        pos_su (numpy.ndarray):  The position of atoms in the supercell.
        skip_distance (float): Skip the neighbor distance under 0.1 (Default).
        ion (list): The selected elements for the neighbor calculation.

    Returns:
        neighbors (list): The neighbors distance of atoms in the 1st shell.
        neighbors_ele (list):  The neighbors element of atoms in the 1st shell.
        neighbors_2 (list): The neighbors distance of atoms in the 2nd shell.
        neighbors_ele_2 (list): The neighbors element of atoms in the 2nd shell.
    """
    neighbors = []
    neighbors_ele = []
    neighbors_2 = []
    neighbors_ele_2 = []
    for i in range(0, len(ele)):  # select center atoms
        distance_temp = []
        neighbors_ele_temp = []
        distance_temp_2 = []
        neighbors_ele_temp_2 = []
        if ele[i] in ion:
            pos_center = pos_new[i]
            for j in range(0, len(ele_su)):  # select neighbor atoms
                if pos_center[0] - cutoff2 <= pos_su[j][0] <= pos_center[0] + cutoff2:
                    if pos_center[1] - cutoff2 <= pos_su[j][1] <= pos_center[1] + cutoff2:
                        if pos_center[2] - cutoff2 <= pos_su[j][2] <= pos_center[2] + cutoff2:
                            if ele_su[j] in ion:
                                temp = cal_distance(i, j, pos_new, pos_su)
                                if skip_distance < temp <= cutoff1:
                                    distance_temp.append(temp)
                                    ele_temp = ele_su[j]
                                    neighbors_ele_temp.append(ele_temp)
                                elif cutoff1 < temp <= cutoff2:
                                    distance_temp_2.append(temp)
                                    ele_temp = ele_su[j]
                                    neighbors_ele_temp_2.append(ele_temp)
        neighbors.append(distance_temp)
        neighbors_ele.append(neighbors_ele_temp)
        neighbors_2.append(distance_temp_2)
        neighbors_ele_2.append(neighbors_ele_temp_2)
        # print('neighbor_shell_1',len(distance_temp), distance_temp, neighbors_ele_temp)
        # print(len(neighbors))
        # print('neighbor_shell_2',len(distance_temp_2), distance_temp_2, neighbors_ele_temp_2)
    return neighbors, neighbors_ele, neighbors_2, neighbors_ele_2


################################################################################################
# plot neighbors


def plot_neighbors(cutoff, neighbors_1, neighbors_2):
    """
    Plot the neighbor distribution of atoms in the cell.

    Args:
        cutoff (float): The maximum cutoff for the neighbor calculation.
        neighbors_1 (list): The neighbors distance of atoms in the 1st shell.
        neighbors_2 (list): The neighbors distance of atoms in the 2nd shell.

    Returns:
        save_plot_data (list): The plot data of neighbor calculation.
    """
    list_x = np.arange(0, cutoff + 0.02, 0.01)
    list_y = np.zeros(len(list_x))
    save_plot_data = [list_x]
    for i in neighbors_1:
        for j in i:
            temp = '%.2f' % j
            for k in range(0, len(list_x)):
                if '%.2f' % list_x[k] == temp:
                    list_y[k] += 1
    for i in neighbors_2:
        for j in i:
            temp = '%.2f' % j
            for k in range(0, len(list_x)):
                if '%.2f' % list_x[k] == temp:
                    list_y[k] += 1
    save_plot_data.append(list_y)
    # list_y_smooth = savgol_filter(list_y,15,10)
    plt.plot(list_x, list_y, '-', linewidth=1.5)
    # plt.plot(list_x, list_y_smooth, '-', linewidth=1.5)
    plt.xlabel('Neighbor distance (Å)', fontsize=11)
    plt.ylabel('Neighbor count', fontsize=11)
    return save_plot_data


def plot_partial_neighbors(cutoff, neighbors_1, neighbors_2, ele, neighbors_ele, neighbors_ele_2, ion1):
    """
    Plot the partial neighbor distribution of atoms in the cell.

    Args:
        cutoff (float): The maximum cutoff for the neighbor calculation.
        neighbors_1 (list): The neighbors distance of atoms in the 1st shell.
        neighbors_2 (list): The neighbors distance of atoms in the 2nd shell.
        ele (list): The element of atoms.
        neighbors_ele (list):  The neighbors element of atoms in the 1st shell.
        neighbors_ele_2 (list): The neighbors element of atoms in the 2nd shell.
        ion1 (list): The selected elements for the neighbor calculation.

    Returns:
        save_plot_data (list): The plot data of partial neighbor calculation.
    """
    # print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    list_x = np.arange(0, cutoff + 0.02, 0.01)
    save_plot_data = [list_x]
    ele_list = copy.deepcopy(ion1)
    ele_list_temp = copy.deepcopy(ele_list)
    for ii in ele_list:
        center_ele = ii
        for jj in ele_list_temp:
            second_ele = jj
            list_y = np.zeros(len(list_x))
            for i in range(len(ele)):
                if ele[i] == center_ele:
                    for j in range(len(neighbors_1[i])):
                        if neighbors_ele[i][j] == second_ele:
                            temp = '%.2f' % neighbors_1[i][j]
                            for k in range(0, len(list_x)):
                                if '%.2f' % list_x[k] == temp:
                                    list_y[k] += 1
                    for j in range(len(neighbors_2[i])):
                        if neighbors_ele_2[i][j] == second_ele:
                            temp = '%.2f' % neighbors_2[i][j]
                            for k in range(0, len(list_x)):
                                if '%.2f' % list_x[k] == temp:
                                    list_y[k] += 1
            plt.plot(list_x, list_y, '-', label=f'{center_ele}-{second_ele}', linewidth=1)
            save_plot_data.append(list_y)
        del ele_list_temp[0]
    plt.legend(frameon=False, fontsize=10.5)
    plt.xlabel('Neighbor distance (Å)', fontsize=11)
    plt.ylabel('Neighbor count', fontsize=11)
    # print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    return save_plot_data


################################################################################################
# plot rdf


def plot_rdf(cutoff, neighbors_1, neighbors_2, particle_number_sum, v):
    """
    Plot the radial distribution function of atoms in the cell.

    Args:
        cutoff (float): The maximum cutoff for the neighbor calculation.
        neighbors_1 (list): The neighbors distance of atoms in the 1st shell.
        neighbors_2 (list): The neighbors distance of atoms in the 2nd shell.
        particle_number_sum (float): The particle number of atoms in cell.
        v (float): The volume of cell.

    Returns:
        save_plot_data (list): The plot data of partial neighbor calculation.
    """
    list_x = np.arange(0, cutoff + 0.02, 0.01)
    list_y = np.zeros(len(list_x))
    for i in neighbors_1:
        for j in i:
            temp = '%.2f' % j
            for k in range(0, len(list_x)):
                if '%.2f' % list_x[k] == temp:
                    list_y[k] += 1
    for i in neighbors_2:
        for j in i:
            temp = '%.2f' % j
            for k in range(0, len(list_x)):
                if '%.2f' % list_x[k] == temp:
                    list_y[k] += 1
    list_y_final = []
    for i in range(0, len(list_x) - 1):
        temp_new = ((list_y[i] / particle_number_sum) /
                    ((4 / 3) * np.pi * (list_x[i + 1] ** 3 - list_x[i] ** 3)) /
                    (particle_number_sum / v))
        list_y_final.append(temp_new)
    list_x_final = np.arange(0, cutoff + 0.01, 0.01)
    save_plot_data = [list_x_final]
    save_plot_data.append(list_y_final)
    plt.plot(list_x_final, list_y_final, '-', linewidth=1.5)
    plt.xlabel('r (Å)', fontsize=11)
    plt.ylabel('g (Å$^{-1}$)', fontsize=11)
    return save_plot_data


def plot_partial_rdf(cutoff, neighbors_1, neighbors_2, ele, neighbors_ele, neighbors_ele_2,
                     ion1, particle_number, v):
    """
    Plot the partial radial distribution function of atoms in the cell.

    Args:
        cutoff (float): The maximum cutoff for the neighbor calculation.
        neighbors_1 (list): The neighbors distance of atoms in the 1st shell.
        neighbors_2 (list): The neighbors distance of atoms in the 2nd shell.
        ele (list): The element of atoms.
        neighbors_ele (list):  The neighbors element of atoms in the 1st shell.
        neighbors_ele_2 (list): The neighbors element of atoms in the 2nd shell.
        ion1 (list): The selected elements for the neighbor calculation.
        particle_number (list): The particle number of elements in ion group.
        v (float): The volume of cell.

    Returns:
        save_plot_data (list): The plot data of partial neighbor calculation.
    """
    # print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    list_x = np.arange(0, cutoff + 0.02, 0.01)
    list_x_final = np.arange(0, cutoff + 0.01, 0.01)
    save_plot_data = [list_x_final]
    ele_list = copy.deepcopy(ion1)
    ele_list_temp = copy.deepcopy(ele_list)
    for ii in ele_list:
        center_ele = ii
        for jj in ele_list_temp:
            second_ele = jj
            list_y = np.zeros(len(list_x))
            for i in range(len(ele)):
                if ele[i] == center_ele:
                    for j in range(len(neighbors_1[i])):
                        if neighbors_ele[i][j] == second_ele:
                            temp = '%.2f' % neighbors_1[i][j]
                            for k in range(0, len(list_x)):
                                if '%.2f' % list_x[k] == temp:
                                    list_y[k] += 1
                    for j in range(len(neighbors_2[i])):
                        if neighbors_ele_2[i][j] == second_ele:
                            temp = '%.2f' % neighbors_2[i][j]
                            for k in range(0, len(list_x)):
                                if '%.2f' % list_x[k] == temp:
                                    list_y[k] += 1
            list_y_final = []
            for j in range(len(ion1)):
                if ion1[j] == center_ele:
                    if center_ele == second_ele:
                        particle_number_2 = particle_number[j]
                    else:
                        for k in range(len(ion1)):
                            if ion1[k] == second_ele:
                                particle_number_2 = particle_number[k]
                    for i in range(0, len(list_x) - 1):
                        temp_new = ((list_y[i] / particle_number[j]) /
                                    ((4 / 3) * np.pi * (list_x[i + 1] ** 3 - list_x[i] ** 3)) /
                                    (particle_number_2 / v))
                        list_y_final.append(temp_new)
            plt.plot(list_x_final, list_y_final, '-', label=f'{center_ele}-{second_ele}', linewidth=1)
            save_plot_data.append(list_y_final)
        del ele_list_temp[0]
    plt.legend(frameon=False, fontsize=10.5)
    plt.xlabel('r (Å)', fontsize=11)
    plt.ylabel('g (Å$^{-1}$)', fontsize=11)
    # print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    return save_plot_data


################################################################################################

def plot_for_rdf(partial_neighbors, cutoff2, neighbors, neighbors_2,
                 ele, neighbors_ele, neighbors_ele_2, ion1, plot_save, cell):
    """
    plot rdf function.

    Args:
        partial_neighbors (bool): Whether calculate the partial neighbors.
        cutoff2 (float): The cutoff of the 2nd shell.
        neighbors (list): The neighbors distance of atoms in the 1st shell.
        neighbors_ele (list):  The neighbors element of atoms in the 1st shell.
        neighbors_2 (list): The neighbors distance of atoms in the 2nd shell.
        neighbors_ele_2 (list): The neighbors element of atoms in the 2nd shell.
        ele (list): The elements of atoms in the cell.
        ion1 (list): The selected elements for the PM-SRO calculation.
        plot_save (bool): Whether save the neighbor plot.
        cell (ase.cell.Cell): The axis vectors of cell.

    Returns:
        save_plot_data (list): The plot data for saving.
    """
    particle_number = []
    for i in ion1:
        temp_count = 0
        for j in range(len(ele)):
            if ele[j] == i:
                temp_count += 1
        particle_number.append(temp_count)
    particle_number_sum = np.sum(particle_number)
    v_temp = [cell[0], cell[1], cell[2]]
    v_temp = np.array(v_temp)
    v = np.linalg.det(v_temp)
    if partial_neighbors:
        # save_plot_data = plot_partial_neighbors(cutoff2, neighbors, neighbors_2, ele,
        #                                         neighbors_ele, neighbors_ele_2, ion1)
        save_plot_data = plot_partial_rdf(cutoff2, neighbors, neighbors_2, ele, neighbors_ele,
                                          neighbors_ele_2, ion1, particle_number, v)
        if plot_save:
            plt.savefig('partial_rdf', dpi=300)
            plt.show()
        else:
            plt.show()
    else:
        # save_plot_data = plot_neighbors(cutoff2, neighbors, neighbors_2)
        save_plot_data = plot_rdf(cutoff2, neighbors, neighbors_2, particle_number_sum, v)
        if plot_save:
            plt.savefig('rdf', dpi=300)
            plt.show()
        else:
            plt.show()
    return save_plot_data


def plot_for_nieighbors(partial_neighbors, cutoff2, neighbors, neighbors_2,
                        ele, neighbors_ele, neighbors_ele_2, ion1, plot_save):
    """
    plot neighbors function.

    Args:
        partial_neighbors (bool): Whether calculate the partial neighbors.
        cutoff2 (float): The cutoff of the 2nd shell.
        neighbors (list): The neighbors distance of atoms in the 1st shell.
        neighbors_ele (list):  The neighbors element of atoms in the 1st shell.
        neighbors_2 (list): The neighbors distance of atoms in the 2nd shell.
        neighbors_ele_2 (list): The neighbors element of atoms in the 2nd shell.
        ele (list): The elements of atoms in the cell.
        ion1 (list): The selected elements for the PM-SRO calculation.
        plot_save (bool): Whether save the neighbor plot.

    Returns:
        save_plot_data (list): The plot data for saving.
    """
    if partial_neighbors:
        save_plot_data = plot_partial_neighbors(cutoff2, neighbors, neighbors_2, ele,
                                                neighbors_ele, neighbors_ele_2, ion1)
        if plot_save:
            plt.savefig('partial_neighbors', dpi=300)
            plt.show()
        else:
            plt.show()
    else:
        save_plot_data = plot_neighbors(cutoff2, neighbors, neighbors_2)
        if plot_save:
            plt.savefig('neighbors', dpi=300)
            plt.show()
        else:
            plt.show()
    return save_plot_data
