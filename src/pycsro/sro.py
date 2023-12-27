import numpy as np


################################################################################################
# SRO parameters

def select_ele(ele):
    """
    Select the types of elements in the cell.

    Args:
        ele (list): The elements of atoms in the cell.

    Returns:
        ele_list (list): The types of elements in the cell.
    """
    ele_list = list(set(ele))
    return ele_list


def cal_average_density(ele, second_ele, ion):
    """
    Calculate the average density of second element.

    Args:
        ele (list): The elements of atoms in the cell.
        second_ele (str): The element type of the second atoms.
        ion (list): The selected elements for the neighbor calculation.

    Returns:
        average_density (float): The average density of the second element.
    """
    count_ele = 0
    sum_ele = 0
    for i in ele:
        if i == second_ele:
            count_ele += 1
        if i in ion:
            sum_ele += 1
    average_density = count_ele / sum_ele
    # print(second_ele, average_density)
    return average_density


def cal_neighbor_density(ele, neighbors_ele, center_ele, second_ele):
    """
    Calculate the second element local density of center atoms.

    Args:
        ele (list): The elements of atoms in the cell.
        center_ele (str): The element type of the center atoms.
        second_ele (str): The element type of the second atoms.
        neighbors_ele (list): The neighbors element of atoms.

    Returns:
        neighbor_density_mean (float): The second element local density of center atoms.
    """
    neighbor_density = []
    for i in range(0, len(ele)):
        count_sec_ele = 0
        if ele[i] == center_ele:
            for j in neighbors_ele[i]:
                if j == second_ele:
                    count_sec_ele += 1
            if len(neighbors_ele[i]) != 0:
                temp = count_sec_ele / len(neighbors_ele[i])
            else:
                temp = 0
            neighbor_density.append(temp)
    neighbor_density_mean = np.mean(neighbor_density)
    # print(neighbor_density_mean)
    return neighbor_density_mean


def cal_pm_sro(neighbor_density, average_density, center_ele, second_ele):
    """
    Calculate the PM-SRO parameter of elements pair.

    Args:
        neighbor_density (float): The second element local density of center atoms.
        average_density (float): The average density of the second element.
        center_ele (str): The element type of the center atoms.
        second_ele (str): The element type of the second atoms.

    Returns:
        pmsro (float): the PM-SRO parameter of elements pair.
    """
    if center_ele == second_ele:
        pmsro = - (neighbor_density - average_density) / (1 - average_density)
    else:
        pmsro = (neighbor_density - average_density) / (0 - average_density)
    # print(pmsro)
    print(f'| {center_ele}-{second_ele} {pmsro}')
    return pmsro


def cal_wc_sro(neighbor_density, average_density, center_ele, second_ele):
    """
    Calculate the WC-SRO parameter of elements pair.

    Args:
        neighbor_density (float): The second element local density of center atoms.
        average_density (float): The average density of the second element.
        center_ele (str): The element type of the center atoms.
        second_ele (str0: The element type of the second atoms.

    Returns:
        wcsro (float): the WC-SRO parameter of elements pair.
    """
    wcsro = 1 - neighbor_density / average_density
    print(f'| {center_ele}-{second_ele} {wcsro}')
    return wcsro

