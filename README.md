# FEM Final Project (ME 254/441) Fall'23
### Afnan Mostafa
### Dec 22, 2023 (Happy Holidays!!)

## System and Sample Meshing:
<img src="https://github.com/afnanmostafa/Finite-Elements-ME254/blob/4d8d06c460e71bf6c80f29a5d0b168ad7df5b5fe/figures/pro1.png" alt="drawing" width="200"/>
![system](https://github.com/afnanmostafa/Finite-Elements-ME254/blob/4d8d06c460e71bf6c80f29a5d0b168ad7df5b5fe/figures/pro1.png)
![meshing](https://github.com/afnanmostafa/Finite-Elements-ME254/blob/b8b44e4284e439990dc749d78de22d6751afccc0/figures/model.png)

### ./src/ subdirectory includes *functions* and *main script* for running FEM simulations in MATLAB.


1. main_script_FEM_AfnanMostafa.m	= 	main MATLAB script

2. IntegrandStiffMatQ4.m		=	function file for calculating integrand

3. JacobianMatQ4.m			=	function file (called by IntegrandStiffMatQ4.m to get Jacobian)

4. GaussQuadQ4.m			=	function file for evaluating Gauss Qudrature integration

5. globalizeStiffMat.m			=	function file for generating global stiffness matrix
		
6. nodalStress.m			=	function file for getting nodal stresses

7. show_displacements.m			=	function file for plotting contours


#### ./src/input-files/ includes Other files are input files necessary for running these codes.

### __Runtime of this code at different architectures:__

1. Fine mesh:	
* Full-int: 	**~60 seconds** (Intel 13th Gen i7, 32 GB RAM);	**~90 seconds** (Mac M2, RAM 24 GB)
* Red-int: 	**~30 seconds** (Intel 13th Gen i7, 32 GB RAM);	**~45 seconds** (Mac M2, RAM 24 GB)

2. Coarse mesh:
* Full-int: 	 **~20 seconds** (Intel 13th Gen i7, 32 GB RAM); **~30 seconds** (Mac M2, RAM 24 GB)
* Red-int: 	 **~10 seconds** (Intel 13th Gen i7, 32 GB RAM); **~15 seconds** (Mac M2, RAM 24 GB)


## Sample Output (displacement contours):
![disp-fine](https://github.com/afnanmostafa/Finite-Elements-ME254/blob/a5be71b8df82a2dcf2b68bcd44262b46dc5f66e6/figures/f-s-d.png)

