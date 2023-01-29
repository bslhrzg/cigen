#!/bin/bash

cp dets_cisd dets	#CISD determinants of 6-31g Water at equilibrium
cp coefs_cisd coefs	#correspondind coefficients of the ground-state CISD wavefunction
cp n_det_cisd n_det	#number of determinants

bindir=../../src/

#We grep parameters from input.nml for consistency
N=`grep "N =" input.nml | awk '{print $3}'` 		#number of spin-orbitals
Ne=`grep "Ne =" input.nml | awk '{print $3}'`		#number of electrons
tol=`grep tol input.nml | awk '{print $3}'`		#pruning: determinants whose coefficient c is such that |c|^2 < 1e(-tol) are discarded
Ns=`grep Nsample input.nml | awk '{print $3}'`		#size of the training sample


#transform determinants from quantum package representation to CIgen representation
#prune
#create training sample
$bindir/qp2cigen.py $N $Ne $tol $Ns

#CIgen program: read parameters from input.nml, read determinants, train the RBM and sample new determinants
$bindir/CI_gen 

nd=$(wc -l < dets)
echo 'Total number of determinants after generation ' $nd

#transform back determinants from CIgen representation to quantum package one
$bindir/cigen2qp.py $N $Ne
