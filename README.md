
# CIgen 


This code implements the methodology in the paper:
"Solving the Schrödinger Equation in the Configuration Space with Generative Machine Learning" 
at "https://arxiv.org/abs/2208.06708"

This is a Fortran 90 implementation of a Restricted Boltmann Machine (RBM)
applied to the Configuration Interaction (CI) framework of quantum chemistry.
The main idea is to take as input an approximate CI wavefunction, use the RBM to learn the joint probability distribution
of the occupation of the spin-orbitals, and generate as output important determinants to improve the accuracy of the current wavefunction.

It is currently interfaced with Quantum Package (QP) ("URL") for the wafunction diagonalization

Contains:

	-src/ directory:
		-the Fortran sources
		-Makefile to compile the Fortran code
		-two python scripts "cigen2qp.py" and "qp2cigen.py", that assure the transformation of the wafunction representation between 
		 QP and CIgen. Those scripts also take care of the wafunction pruning and the creation of the RBM training set.
	-example/ directory:
		-test/ directory:
			-test wavefunction to test the code
		-water_631g_quantum_package/ directory:
			-bash script example to use the code with QP
			-qp_plugins/ directory: plugins added to QP to read/write the wavefunction

Dependencies:

	-To compile the code, a Fortran compiler is required (only gfortran was tested), and openmp for the parallelized pieces of code

	-To use the code in an actual situation, a code is requiered that can read/write and diagonalize a CI wavefunction.
	 For this project Quantum Package (v2) was used (https://quantum-package.readthedocs.io/en/master/index.html), 
	 together with added plugin ml_ci to read/write the QP wavefunction. Those can be found in example/water_631g_quantum_package/qp_plugins/



## Aknowledgments:

	In our implementation we would like to acknowledge the use of the following sorting routines:
		
		-function rargsort_int from "https://github.com/certik/fortran-utils/blob/master/src/sorting.f90"
		 based on code written by John E. Pask, LLNL.
		 Licence: MIT Licence


		-subroutine quicksort written by t-nisie, see "https://gist.github.com/t-nissie/479f0f16966925fa29ea"
		 Licence: GPLv3



		-routines uniinv, nearless and I_unista were written by Michel Olagnon and part of:
		 "ORDERPACK 2.0 -- Unconditional, Unique, and Partial Ranking, Sorting, and Permutation Downloadable Fortran 90 source code"
		 available at http://www.fortran-2000.com/rank/


## Author

- Basile Herzog: basile.herzog@univ-lorraine.fr
