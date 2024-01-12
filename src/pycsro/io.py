import os
from ase.io import read
from pymatgen.core import Structure


################################################################################################
# read file

def readfile(file_name):
    """
    Read file with ase model.

    Args:
        file_name (str): The absolute path of input file.

    Returns:
        atoms (ase.atoms.Atoms): The configuration of crystal structure.
    """
    atoms = read(file_name)
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


def check_readfile(file_name):
    """
    Read file with ase model.

    Args:
        file_name (str): The absolute path of input file.

    Returns:
        file_name (str): The absolute path of input file.
        initial_file_name (str): The initial absolute path of input file.
    """
    initial_file_name = file_name
    try:
        atoms = read(file_name)
    except Exception:
        structure = Structure.from_file(file_name)
        structure.to(filename="_temp.cif")
        atoms = read("_temp.cif")
        file_name= "_temp.cif"
    return file_name, initial_file_name


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
        for i in range(len(save_plot_data[0])):
            for j in range(len(save_plot_data)):
                file.write(f'{save_plot_data[j][i]} ')
            file.write(f'\n')


def save_file(save_name, dual_cutoff, wc_list, wc_list_2, save_plot_data, xyz):
    """
    Save file function.

    Args:
        save_name (str): The name of the saved file.
        dual_cutoff (bool): Whether the cutoff1 and cutoff2 existed at the same time.
        wc_list (list): the SRO parameter save list of the 1st shell.
        wc_list_2 (list): the SRO parameter save list of the 2nd shell.
        save_plot_data (list): data to be saved
        xyz (bool): Whether the format of input file is XYZ.
    """
    if save_name:  # savefile
        save_wcp_2(os.path.abspath(save_name), wc_list)
        if dual_cutoff:
            save_wcp_2(os.path.abspath(save_name), wc_list_2)
        end = ['# +-----------------------------------------------------------------------------']
        save_wcp_2(os.path.abspath(save_name), end)
        rdf = ['\n',
               '#### start rdf data',
               '#r (A)  g(r) (A^{-1})']
        neighbors = ['\n',
                     '#### start neighbors data',
                     '#Neighbor distance (A)  Neighbor count']
        if xyz:
            save_wcp_2(os.path.abspath(save_name), neighbors)
        else:
            save_wcp_2(os.path.abspath(save_name), rdf)
        save_plot(save_name, save_plot_data)


################################################################################################


def print_settings(ion1, cutoff1, cutoff2, file_name, save_name, skip_distance, dual_cutoff, cal_same_pair, safe_mode,
                   xyz, initial_file_name):
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
        xyz (bool): Whether the format of input file is XYZ.
        initial_file_name (str): The initial absolute path of input file.
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
    if file_name == "_temp.cif":
        file_name = initial_file_name
    print(f'+-----------------------------------------------------------------------------\n'
          f'| Element group: {ion1}\n'
          f'| Cutoff for the 1st shell: {cutoff1} Å       Cutoff for the 2nd shell: {cutoff2_print}\n'
          f'| Read file: {file_name}       Save file: {save_name}\n'
          f'| Skip neighbor distance under {skip_distance} Å\n'
          f'| Calculate same pair: {cal_same_pair_print}       Safe mode: {safe_mode_print}\n'
          f'+-----------------------------------------------------------------------------')
    if cutoff2_print == f'{cutoff2} Å':
        cutoff2_print = f'{cutoff2} A'
    settings = ['# +-----------------------------------------------------------------------------',
                f'# | Element group: {ion1}',
                f'# | Cutoff for the 1st shell: {cutoff1} A       Cutoff for the 2nd shell: {cutoff2_print}',
                f'# | Read file: {file_name}       Save file: {save_name}',
                f'# | Skip neighbor distance under {skip_distance} A',
                f'# | Calculate same pair: {cal_same_pair_print}       Safe mode: {safe_mode_print}',
                '# +-----------------------------------------------------------------------------']
    if save_name:
        save_wcp(os.path.abspath(save_name), settings)
    if skip_distance < 0.1 or safe_mode or xyz:
        if skip_distance < 0.1:
            print(f'| WARNING！Calculateing neighbor distance under 0.1 may cause problems.')
            settings_2 = [f'# | WARNING！Calculateing neighbor distance under 0.1 may cause problems.']
            if save_name:
                save_wcp_2(os.path.abspath(save_name), settings_2)
        if safe_mode:
            print('| Safe Model! SRO calculation without selecting atoms in supercell.\n'
                  '| Which may use much more time...')
            settings_3 = ['# | Safe Model! SRO calculation without selecting atoms in supercell.',
                          '# | Which may use much more time...']
            if save_name:
                save_wcp_2(os.path.abspath(save_name), settings_3)
        if xyz:
            print('| Reading the XYZ file!')
            settings_4 = ['# | Reading the XYZ file!']
            if save_name:
                save_wcp_2(os.path.abspath(save_name), settings_4)
        print('+-----------------------------------------------------------------------------')
        end = ['# +-----------------------------------------------------------------------------']
        if save_name:
            save_wcp_2(os.path.abspath(save_name), end)