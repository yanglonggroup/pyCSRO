import os
from ase.io import read


################################################################################################
# read file

def readfile(filename):
    """
    Read file with ase model.

    Args:
        filename (str): The absolute path of input file.

    Returns:
        atoms (ase.atoms.Atoms): The configuration of crystal structure.
    """
    atoms = read(filename)
    return atoms


def get_config(atoms):
    """
    Get configs from the file.

    Args:
        atoms (ase.atoms.Atoms): The configuration of crystal structure.

    Returns:
        cell (ase.cell.Cell): The axis vectors of cell.
        pos (numpy.ndarray): The position of atoms.
        ele (list): The element of atoms.
    """
    cell = atoms.get_cell()
    pos = atoms.get_positions()
    ele = atoms.get_chemical_symbols()
    return cell, pos, ele


################################################################################################
# save function

def save_wcp(save_name, list):
    """
    Save the PM-SRO parameter in a new file.

    Args:
        save_name (str): The name of the saved file.
        list (list): List to be saved
    """
    with open(f'{save_name}', 'w', encoding='utf-8') as file:
        for i in list:
            file.write(f'{i}\n')


def save_wcp_2(save_name, list):
    """
    Save the PM-SRO parameter in a old file.

    Args:
        save_name (str): The name of the saved file.
        list (list): List to be saved
    """
    with open(f'{save_name}', 'a') as file:
        for i in list:
            file.write(f'{i}\n')


def save_plot(save_name, save_plot_data):
    """
    Save the plot data in a old file.

    Args:
        save_name (str): The name of the saved file.
        save_plot_data (list): data to be saved
    """
    with open(f'{save_name}', 'a') as file:
        file.write(f'\n')
        for i in range(len(save_plot_data[0])):
            for j in range(len(save_plot_data)):
                file.write(f'{save_plot_data[j][i]} ')
            file.write(f'\n')


def print_settings(ion1, cutoff1, cutoff2, file_name, save_name, skip_distance, dual_cutoff, cal_same_pair, safe_mode):
    """
    Print and save the details of PM-SRO calculation.

    Args:
        ion1 (list): The selected elements for the PM-SRO calculation.
        cutoff1 (float): The cutoff of thr 1st shell.
        cutoff2 (float): The cutoff of the 2nd shell.
        file_name (str): The absolute path of input file.
        save_name (str): The absolute path of saved file.
        skip_distance (float): Skip the neighbor distance under 0.1 (Default).
        dual_cutoff (str): Whether the cutoff1 and cutoff2 existed at the same time.
        cal_same_pair (str): Whether calculate the wcp of same elements but different center atoms.
        safe_mode (str): Whether use the supercell selection function, which can reduce the calculation time.
    """
    if dual_cutoff:
        cutoff2_print = f'{cutoff2} Å'
    else:
        cutoff2_print = None
    if cal_same_pair:
        cal_same_pair_print = 'No'
    else:
        cal_same_pair_print = 'Yes'
    if safe_mode:
        safe_mode_print = 'Yes'
    else:
        safe_mode_print = 'No'
    print(f'+-----------------------------------------------------------------------------\n'
          f'| Element group: {ion1}\n'
          f'| Cutoff for the 1st shell: {cutoff1} Å       Cutoff for the 2nd shell: {cutoff2_print}\n'
          f'| Read file: {file_name}       Save file: {save_name}\n'
          f'| Skip neighbor distance under {skip_distance} Å\n'
          f'| Calculate same pair: {cal_same_pair_print}       Safe mode: {safe_mode_print}\n'
          f'+-----------------------------------------------------------------------------')
    if skip_distance < 0.1 or safe_mode:
        if skip_distance < 0.1:
            print(f'| WARNING！Calculateing neighbor distance under {skip_distance} may cause problems.')
        if safe_mode:
            print('| Safe Model! SRO calculation without selecting atoms in supercell.\n'
                  '| Which may use much more time...')
        print('+-----------------------------------------------------------------------------')
    settings = ['+-----------------------------------------------------------------------------',
                f'| Element group: {ion1}',
                f'| Cutoff for the 1st shell: {cutoff1} Å       Cutoff for the 2nd shell: {cutoff2_print}',
                f'| Skip neighbor distance under {skip_distance} Å',
                f'| Calculate same pair: {cal_same_pair_print}       Safe mode: {safe_mode_print}',
                '+-----------------------------------------------------------------------------']
    if save_name:
        save_wcp(os.path.abspath(save_name), settings)










