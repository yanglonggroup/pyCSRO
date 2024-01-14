#!/usr/bin/env python
import sys
import os

# if __name__ == "__main__":
#     # show output results from every test function
#     args = ["-v"]
#     # show the message output for skipped and expected failure tests
#     if len(sys.argv) > 1:
#         args.extend(sys.argv[1:])
#     print("pytest arguments: {}".format(args))
#     # call pytest and exit with the return code from pytest so that
#     # travis will fail correctly if tests fail
#     exit_res = pytest.main(args)
#     sys.exit(exit_res)


def get_path(folder, file_name):
    current_path = os.path.dirname(os.path.abspath(__file__))
    base_path = os.path.dirname(current_path)  # get into upper folder
    target_path = os.path.join(base_path, folder, file_name)
    return target_path


current_path = os.path.dirname(os.path.abspath(__file__))
base_path = os.path.dirname(current_path)  # get into upper folder
target_path = os.path.join(base_path, 'src', 'pycsro')
target_path_2 = os.path.join(base_path, 'src')
sys.path.append(target_path)  #
sys.path.append(target_path_2)  #
from pycsro.main import *


current_path = os.path.dirname(os.path.abspath(__file__))
base_path = os.path.dirname(current_path)
target_path = os.path.join(base_path, 'example', 'PbTe_pcell.vasp')
target_path_2 = os.path.join(base_path, 'example', 'Cu.xyz')


def test_main():
    ion1 = 'Pb Te'
    cutoff1 = 3.7
    ion1_test, cutoff1_test, target_path_test, sro_test = run_pycsro_pmsro(ion1, cutoff1, target_path)
    assert ion1_test == ion1.split(' ')
    assert cutoff1_test == cutoff1
    assert target_path_test == target_path


def test_wc_sro():
    ion1 = 'Pb'
    cutoff1 = 5
    ion1_test, cutoff1_test, target_path_test, sro_test = run_pycsro_pmsro(ion1, cutoff1, target_path)
    assert sro_test == 'Pb-Pb 0.0'


def test_pm_sro():
    ion1 = 'Pb Te'
    cutoff1 = 3.7
    ion1_test, cutoff1_test, target_path_test, sro_test = run_pycsro_pmsro(ion1, cutoff1, target_path)
    assert sro_test == 'Pb-Pb 1.0'


def test_file_read():
    atom_test = readfile(target_path)
    assert atom_test == read(target_path)


def test_xyz_file_read():
    atoms = readfile(target_path_2)
    cell, pos, ele = get_config(atoms)
    if cell.any():
        xyz = False
    else:
        xyz = True
    assert xyz == True


def test_check_readfile():
    file_name, initial_file_name = check_readfile(target_path)
    assert file_name == initial_file_name


def test_axis():
    atoms = readfile(target_path)
    cell, pos, ele = get_config(atoms)
    axis = get_axis(cell)
    assert axis[0] == 4.563667163090972


def test_neighbor():
    ion1 = ['Pb', 'Te']
    cutoff1 = 3.7
    cutoff2 = 5
    safe_mode = False
    skip_distance = 0.1
    cell, pos_new, ele, cell_su, pos_su, ele_su, plane_n_su = readfile_pmsro(file_name=target_path, cutoff2=cutoff1)
    pos_su_new, ele_su_new = safemode(safe_mode, pos_su, ele_su, cutoff2, cell, cell_su, plane_n_su)
    neighbors_1, neighbors_ele, neighbors_2, neighbors_ele_2 = \
        cal_neighbors(ele, pos_new, ele_su_new, cutoff1, cutoff2, pos_su_new, skip_distance, ion1)
    assert neighbors_1 == [[3.2269999981, 3.2269999981, 3.2269999981, 3.2269999981, 3.2269999981, 3.2269999981],
                           [3.2269999981, 3.2269999981, 3.2269999981, 3.2269999981, 3.2269999981, 3.2269999981]]
    assert neighbors_ele == [['Te', 'Te', 'Te', 'Te', 'Te', 'Te'], ['Pb', 'Pb', 'Pb', 'Pb', 'Pb', 'Pb']]
