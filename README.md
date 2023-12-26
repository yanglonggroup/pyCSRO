pyCSRO: a python package for calculating chemical short range ordering.
===================================================================================================

pyCSRO is a python package using pairwise multi-component short-range order (PM-SRO) parameters to discover the chemical short range ordering in materials, including crystalline, amorphous, and high-entropy alloys, etc.
The PM-SRO parameter was extended from the Warren-Cowley short-range order (WC-SRO) to describe the local distribution of element pairs in the multicomponent system, which is defined as:

$$\alpha_{ij}^m = \frac{p_{ij}^m - X_j}{\delta_{ij} - X_j},$$ 

where $p_{ij}^m$ is the average probability of finding a j-type atom around an i-type atom in the $m$-th shell, and $X_j$ is the overall concentration of A atoms in the system. The $\delta_{ij}$ equals 1 if $i = j$ and 0 if $i ​\neq j$.
A negative value of $\alpha_{ij}^m$ indecates the tendency of mixing of i and j atoms, whereas a positive one suggests the tendency toward segregation of i and j atoms. And the value of $\alpha_{ij}^m$ would be zero if the i and j atoms are randomly distributed.

Note that 

Required Dependencies:
------------
* Python 3.9+
* matplotlib
* numpy
* ase


Installation
------------
For python package and conda users:
```
pip install pycsro
```

Usage
--------
```
from pycsro import *
pycsro.run_pycsro_pmsro(ion1, cutoff1, filename, cutoff2, savename, skip_distance, plotsave, cal_same_pair, safe_mode, partial_neighbors)
```

- ion1: The selected elements for the PM-SRO calculation. (Required, range: elements in the structure model, type: str)

- cutoff1: The cutoff of thr 1st shell. (Required, range: positive number, type: float)

- cutoff2: The cutoff of the 2nd shell. (Cutoff2 should be equal or greater than cutoff1, range: positive number, type: float)

- file_name: The absolute path of input file. (Required, type: str)

- save_name: The absolute path of saved file. (type: str)

- skip_distance: Skip the neighbor distance under this parameter. (Default: 0.1, range: positive number, type: float)

- cal_same_pair: Whether calculate the wcp of same elements but different center atoms. (Default: Yes, range: Yes or No, type: str)

- safe_mode: Whether use the supercell selection function, which can reduce the calculation time. (Default: No, range: Yes or No, type: str)

- plot_save: Whether save the neighbor plot. (Default: No, range: Yes or No, type: str)

- partial_neighbors: Whether plot the partial neighbor distribution of atoms in the cell. (Default: No, range: Yes or No, type: str)

Attention: 
1. The calculation requires the user to input parameters of ion1, cutoff1, and filename.
2. Supported input file formats include CIF, POSCAR.
3. For compounds, cations and anions should be placed in different ion groups, and calculated separately.
4. You can try a larger cutoff value at the first time, and then adjust the cutoff value to the trough of the neighbor distribution.
5. The second shell was defined as cutoff1 to cutoff2.


Example
--------
Example models are placed at `/example/`.

The basic command for running Ni system:
```
pycsro.run_pycsro_pmsro(ion1='Ni', cutoff1=3, filename='XXX/pycsro/example/Ni.vasp')
```

The output:
```
+-----------------------------------------------------------------------------
| Element group: ['Ni']
| Cutoff for the 1st shell: 3 Å       Cutoff for the 2nd shell: None
| Read file: Ni.vasp       Save file: None
| Skip neighbor distance under 0.1 Å
| Calculate same pair: Yes       Safe mode: No
+-----------------------------------------------------------------------------
| The PM-SRO parameter for ['Ni'] element group in the 1st shell
| Ni-Ni 0.0
+-----------------------------------------------------------------------------
```

Contacts
--------

For more information, please email Prof. Long Yang at long_yang@tongji.edu.cn
 