import math
import numpy as np
import copy
from operator import abs


################################################################################################
# get default settings

def default_settings(ion1, cutoff1, cutoff2, skip_distance, plot_save, cal_same_pair, safe_mode):
    """
    Get the default settings from inputs.

    Args:
        ion1 (str): The selected elements for the PM-SRO calculation.
        cutoff1 (float): The cutoff of thr 1st shell.
        cutoff2 (float): The cutoff of the 2nd shell.
        skip_distance (float): Skip the neighbor distance under 0.1 (Default).
        cal_same_pair (str): Whether calculate the wcp of same elements but different center atoms.
        safe_mode (str): Whether use the supercell selection function, which can reduce the calculation time.
        plot_save (str): Whether save the neighbor plot.

    Returns:
        ion1 (list): The selected elements for the PM-SRO calculation.
        cutoff1 (float): The cutoff of thr 1st shell.
        cutoff2 (float): The cutoff of the 2nd shell.
        skip_distance (float): Skip the neighbor distance under 0.1 (Default).
        cal_same_pair (bool): Whether calculate the wcp of same elements but different center atoms.
        safe_mode (bool): Whether use the supercell selection function, which can reduce the calculation time.
        plot_save (bool): Whether save the neighbor plot.
        dual_cutoff (bool): Whether the cutoff1 and cutoff2 existed at the same time.
        single_ele (bool): Whether use the WC-SRO to calculate the only one input element.
    """
    ion1 = ion1.split(' ')
    if len(ion1) == 1:
        single_ele = True
    else:
        single_ele = False
    dual_cutoff = True
    if skip_distance is None:
        skip_distance = 0.1  # skip neighbor distance under
    if cal_same_pair is None or cal_same_pair.lower() == 'yes' or cal_same_pair.lower() == 'y':
        cal_same_pair = False  # whether calculate the wcp of same elements but different center atoms
    elif cal_same_pair.lower() == 'no' or cal_same_pair.lower() == 'n':
        cal_same_pair = True
    if safe_mode is None or safe_mode.lower() == 'no' or safe_mode.lower() == 'n':
        safe_mode = False  # whether use the supercell selection function, which can reduce the calculation time
    elif safe_mode.lower() == 'yes' or safe_mode.lower() == 'y':
        safe_mode = True
    if cutoff2 is None:
        cutoff2 = cutoff1  # only calculate the 2nd shell
        dual_cutoff = False
    if plot_save is None or plot_save.lower() == 'no' or plot_save.lower() == 'n':
        plot_save = False  # Whether save the neighbor plot. Default: No
    return ion1, cutoff1, cutoff2, skip_distance, plot_save, cal_same_pair, safe_mode, dual_cutoff, single_ele


################################################################################################
# get config

def get_axis(cell):
    """
    Get axis lengths from the cell.

    Args:
        cell (ase.cell.Cell): The axis vectors of cell.

    Returns:
        axis (list): The length of the axis.
    """
    axis = []
    for i in cell:  # calculate axis
        temp = (i[0] ** 2 + i[1] ** 2 + i[2] ** 2) ** 0.5
        axis.append(temp)
    # print(axis)
    return axis


def periodicity(cutoff, axis, atoms):
    """
    Set up a supercell for the neighbors calculation.

    Args:
        cutoff (float): The maximum cutoff for the neighbors calculation.
        axis (list): The length of the axis.
        atoms (ase.atoms.Atoms): The configuration of crystal structure.

    Returns:
        supercell (ase.atoms.Atoms): The configuration of supercell crystal structure.
    """
    if cutoff < min(axis):
        supercell = atoms * (3, 3, 3)
        # print('type 1')
        # write("ouput_type1.vasp", sort(supercell), direct=True)
    else:
        temp = math.ceil(cutoff / min(axis)) * 2 + 1
        supercell = atoms * (temp, temp, temp)
        # print('type 2')
        # write("ouput_type2.vasp", sort(supercell), direct=True)
    return supercell


def move_original_positions(pos, cell, cell_su):
    """
    Move the original cell's position to the center of the supercell.

    Args:
        pos (numpy.ndarray): The position of atoms.
        cell (ase.cell.Cell): The axis vectors of the cell.
        cell_su (ase.cell.Cell): The axis vectors of the supercell.

    Returns:
        pos_new (numpy.ndarray): The new position of the cell.
    """
    su_center_moves = []
    for i in range(0, 3):
        temp = ((cell_su[0][i] + cell_su[1][i] + cell_su[2][i]) - (cell[0][i] + cell[1][i] + cell[2][i])) * 0.5
        su_center_moves.append(temp)
    # print(su_center_moves)
    pos_new = copy.deepcopy(pos)
    for i in pos_new:
        i[0] += su_center_moves[0]
        i[1] += su_center_moves[1]
        i[2] += su_center_moves[2]
    return pos_new


def cal_plane_normal_vector(cell):
    """
    Calculate the normalized plane normal vector for axis planes.

    Args:
        cell (ase.cell.Cell): The axis vectors of the cell.

    Returns:
        plane_n_new (list): The normal vector of axis planes.
    """
    plane_list = [[cell[1], cell[2]], [cell[2], cell[0]], [cell[0], cell[1]]]
    plane_n = []
    for i in plane_list:
        temp = np.cross(i[0], i[1])
        plane_n.append(list(temp))
    # print(plane_n)
    sin_zero = []
    for i in plane_n:
        temp = []
        for j in i:
            if j != 0:
                temp.append(j)
        sin_zero.append(temp)
    # print(sin_zero)
    sin_zero_min = []
    for i in sin_zero:
        temp = []
        for j in i:
            temp.append(abs(j))
        sin_zero_min.append(min(temp))
    # print(sin_zero_min)
    plane_n_new = []
    for i in range(0, 3):
        temp = plane_n[i]
        for j in range(0, 3):
            temp[j] = temp[j] / sin_zero_min[i]
        plane_n_new.append(temp)
    # print(plane_n_new)
    return plane_n_new


def cal_distance_plane(pos_single, plane_n_single):
    """
    Calculate the distance between a vector and a plane.

    Args:
        pos_single (numpy.ndarray): The position of a single atom.
        plane_n_single (list): The normal vector of an axis plane.

    Returns:
        distance (float): The distance between an atom and a plane.
    """
    angle = cal_vectors_angle(pos_single, plane_n_single)
    distance = np.linalg.norm(pos_single) * np.cos(angle / 180 * np.pi)
    return distance


def cal_vectors_angle(pos_single, plane_n_single):
    """
    Get the angle between two vectors.

    Args:
        pos_single (numpy.ndarray): The position of a single atom.
        plane_n_single (list): The normal vector of an axis plane.

    Returns:
        angle (float): The angle between an atom and a normal vector of an axis plane.
    """
    angle_temp = np.arccos(np.dot(pos_single, plane_n_single) /
                           (np.linalg.norm(pos_single) * np.linalg.norm(plane_n_single))) / (0.5 * np.pi) * 90
    temp = '%.2f' % angle_temp
    angle = float(temp)
    # print(angle)
    return angle


def select_supercell(pos_su, ele_su, cutoff, cell, cell_su, plane_n_su):
    """
    Select useful atoms in the supercell.

    Args:
        pos_su (numpy.ndarray): The position of atoms in the supercell
        ele_su (numpy.ndarray): The element of atoms in the supercell
        cutoff (float): he maximum cutoff for the neighbors calculation.
        cell (ase.cell.Cell): The axis vectors of the cell.
        cell_su (ase.cell.Cell): The axis vectors of the supercell.
        plane_n_su (list): The normal vector of axis planes in the supercell.

    Returns:
        pos_su_new (numpy.ndarray): The selected positions in the supercell.
        ele_su_new (numpy.ndarray): The selected element in the supercell.
    """
    distance_limit = []
    pos_su_new = []
    ele_su_new = []
    axis_new = [cal_distance_plane(cell[0], plane_n_su[0]),
                cal_distance_plane(cell[1], plane_n_su[1]),
                cal_distance_plane(cell[2], plane_n_su[2])]
    axis_su_new = [cal_distance_plane(cell_su[0], plane_n_su[0]),
                   cal_distance_plane(cell_su[1], plane_n_su[1]),
                   cal_distance_plane(cell_su[2], plane_n_su[2])]
    for i in range(0, 3):
        distance_limit_temp = ((axis_su_new[i] - axis_new[i]) / 2) - cutoff
        distance_limit.append(distance_limit_temp)
    limit = [[distance_limit[0], axis_su_new[0] - distance_limit[0]],
             [distance_limit[1], axis_su_new[1] - distance_limit[1]],
             [distance_limit[2], axis_su_new[2] - distance_limit[2]]]
    # print(limit)
    distance_table = []
    for i in range(0, len(ele_su)):
        if pos_su[i][0] == pos_su[i][1] == pos_su[i][2] == 0:
            distance_table.append([0, 0, 0])
        else:
            distance_table.append([cal_distance_plane(pos_su[i], plane_n_su[0]),
                                   cal_distance_plane(pos_su[i], plane_n_su[1]),
                                   cal_distance_plane(pos_su[i], plane_n_su[2])])
        # print(distance_table[-1])
        if limit[0][0] <= distance_table[-1][0] <= limit[0][1]:
            if limit[1][0] <= distance_table[-1][1] <= limit[1][1]:
                if limit[2][0] <= distance_table[-1][2] <= limit[2][1]:
                    # print('ok')
                    pos_su_new.append(pos_su[i])
                    ele_su_new.append(ele_su[i])
    # print(len(ele_su_new), pos_su_new)
    # print(len(ele_su_new))
    return pos_su_new, ele_su_new






