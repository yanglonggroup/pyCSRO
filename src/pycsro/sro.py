import numpy as np
import copy


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


################################################################################################

def cal_sro_default_settings(ion1):
    """
    Get the default setting of SRO calculation.

    Args:
        ion1 (list): The calculated elements list.

    Returns:
        wc_list (list): the SRO parameter save list of the 1st shell.
        wc_list_2 (list): the SRO parameter save list of the 2nd shell.
        ele_list (list): The center elements list.
        ele_list_temp (list): The second elements list for the 1st shell calculation.
        ele_list_temp_2 (list): The second elements list for the 2nd shell calculation.
    """
    wc_list = [f'# The PM-SRO parameter for {ion1} element group in the 1st shell']
    wc_list_2 = ['# +-----------------------------------------------------------------------------',
                 f'# The PM-SRO parameter for {ion1} element group in the 2nd shell']
    ele_list = copy.deepcopy(ion1)
    ele_list_temp = copy.deepcopy(ele_list)
    ele_list_temp_2 = copy.deepcopy(ele_list)
    print(f'| The PM-SRO parameter for {ion1} element group in the 1st shell')
    return wc_list, wc_list_2, ele_list, ele_list_temp, ele_list_temp_2


def single_sro_cal_fun(ion1, ele, neighbors_ele, neighbors_ele_2, wc_list, wc_list_2, dual_cutoff):
    """
    Calculate the SRO of single element input.

    Args:
        ion1 (list): The calculated elements list.
        ele (list): The element of atoms.
        neighbors_ele (list):  The neighbors element of atoms in the 1st shell.
        neighbors_ele_2 (list): The neighbors element of atoms in the 2nd shell.
        wc_list (list): the SRO parameter save list of the 1st shell.
        wc_list_2 (list): the SRO parameter save list of the 2nd shell.
        dual_cutoff (bool): Whether the cutoff1 and cutoff2 existed at the same time.

    Returns:

        wc_list (list): the SRO parameter save list of the 1st shell.
        wc_list_2 (list): the SRO parameter save list of the 2nd shell.
    """
    center_ele = ion1[0]
    second_ele = ion1[0]
    average_density = 1
    neighbor_density = cal_neighbor_density(ele, neighbors_ele, center_ele, second_ele)
    wc_list_temp = cal_wc_sro(neighbor_density, average_density, center_ele, second_ele)
    wc_list.append(f'{center_ele}-{second_ele} {wc_list_temp}')
    print(f'+-----------------------------------------------------------------------------')
    if dual_cutoff:
        print(f'| The PM-SRO parameter for {ion1} element group in the 2nd shell')
        neighbor_density = cal_neighbor_density(ele, neighbors_ele_2, center_ele, second_ele)
        wc_list_temp = cal_wc_sro(neighbor_density, average_density, center_ele, second_ele)
        wc_list_2.append(f'{center_ele}-{second_ele} {wc_list_temp}')
        print(f'+-----------------------------------------------------------------------------')
    return wc_list, wc_list_2


def pmsro_cal_fun(wc_list, wc_list_2, ele_list, ele_list_temp, ele_list_temp_2, ele, ion1,
                  neighbors_ele, neighbors_ele_2, cal_same_pair, dual_cutoff):
    """
    Calculate the SRO of multiple elements input.

    Args:
        wc_list (list): the SRO parameter save list of the 1st shell.
        wc_list_2 (list): the SRO parameter save list of the 2nd shell.
        ele_list (list): The center elements list.
        ele_list_temp (list): The second elements list for the 1st shell calculation.
        ele_list_temp_2 (list): The second elements list for the 2nd shell calculation.
        ele (list): The element of atoms.
        ion1 (list): The calculated elements list.
        neighbors_ele (list):  The neighbors element of atoms in the 1st shell.
        neighbors_ele_2 (list): The neighbors element of atoms in the 2nd shell.
        cal_same_pair (bool): Whether calculate the wcp of same elements but different center atoms.
        dual_cutoff (bool): Whether the cutoff1 and cutoff2 existed at the same time.

    Returns:

        wc_list (list): the SRO parameter save list of the 1st shell.
        wc_list_2 (list): the SRO parameter save list of the 2nd shell.
    """
    for i in ele_list:  # calculate pm-sro parameter for 1st shell
        center_ele = i
        for j in ele_list_temp:
            second_ele = j
            average_density = cal_average_density(ele, second_ele, ion1)
            neighbor_density = cal_neighbor_density(ele, neighbors_ele, center_ele, second_ele)
            wc_list_temp = cal_pm_sro(neighbor_density, average_density, center_ele, second_ele)
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
                average_density = cal_average_density(ele, second_ele, ion1)
                neighbor_density = cal_neighbor_density(ele, neighbors_ele_2, center_ele, second_ele)
                wc_list_temp = cal_pm_sro(neighbor_density, average_density, center_ele, second_ele)
                wc_list_2.append(f'{center_ele}-{second_ele} {wc_list_temp}')
            if cal_same_pair:
                del ele_list_temp_2[0]
        print(f'+-----------------------------------------------------------------------------')
    return wc_list, wc_list_2
