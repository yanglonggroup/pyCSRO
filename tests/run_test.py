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
from pycsro.main import run_pycsro_pmsro


def test_pbte():
    ion1 = 'Pb Te'
    cutoff1 = 3
    current_path = os.path.dirname(os.path.abspath(__file__))
    base_path = os.path.dirname(current_path)
    target_path = os.path.join(base_path, 'example', 'PbTe_pcell.vasp')
    ion1_test, cutoff1_test, file_name_test = run_pycsro_pmsro(ion1, cutoff1, target_path)
    assert ion1_test == ion1.split(' ')
    assert cutoff1_test == cutoff1
    assert file_name_test == target_path
