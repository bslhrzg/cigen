#!/bin/env python

##### QP ----> CIgen #####


#### CIgen Slater Determinants are decimal integers ####
#### whose corresponding binary representation are  ####
#### the spin-orbital occupations, for example with ####
#### 4 electrons among 6 spin-orbitals, the HF      ####
#### Slater Determinant is 111100 wtih decimal rep  ####
#### 60.                                            ####
#### Now in Quantum Package, they use a decimal rep ####
#### as well but separate spin up from spin down    ####
#### so this would be 110  and 110                  ####
#### They use a reversed order (011,011) = (11,11)  ####
#### giving (3,3) in decimal representation         ####

import numpy as np
import sys
    
def QP2CIgen(a,b,N) : 
    ########################################################
    #### This function transform the determinant (a,b)  ####
    #### from QP into the form suitable for CIgen       ####
    ########################################################

    def decabtocbs(a,b,N):
        # take the decimal rep of spin up and down (a,b)
        # create the corresponding binary string by merging them
        # reverse the order by construction
        c=''
        n=int(N/2)
        for k in range(n):
            ca = a & 1 #is odd ?
            cb = b & 1
            a, b = a // 2, b // 2
            c+='1' if ca else '0'
            c+= '1' if cb else '0'
        return c
    
    bsc=decabtocbs(a,b,N)
    
    c=int(bsc,2)

    return c
    
def CIgen2QP(c,N) :
    ########################################################
    #### This function transform the determinant (c)    ####
    #### from CIgen into the form (a,b) suitable for QP ####
    ########################################################
    
    
    def decctoabbs(c,N):
        # take the dec=''imal rep of determinant (c)
        # create the corresponding binary strings a,b for
        # spin up and down
        a,b='',''
        for k in range(int(N/2)):
            ac = c & 1 #is odd ?
            c = c // 2
            a+='1' if ac else '0'
            
            ab = c & 1 
            c = c // 2
            b+='1' if ab else '0'
        return a,b
    
    bsa,bsb=decctoabbs(c, N)
    a=int(bsa,2)
    b=int(bsb,2)
    
    return b,a

def dec2bin(c,N):
    ########################################################
    #### This function transform the determinant (c)    ####
    #### from CIgen into a binary array for the training####
    ########################################################    
   
     
    b=[int(x) for x in bin(c)[2:]] #2: because of the '0b' 

    while len(b) < N : b=[0]+b #need to pad low unoccupdied orbitals else 001 is written 1
    return b


def create_training_sample(N,Ns,dets,coefs) :
    ########################################################
    #### Create a training sample dets_s of binary array####
    #### from dets and coefs, while excludiing HF       ####
    ######################################################## 
    
    cc=coefs*coefs # Probability
    cc=cc[1:] # exclude Hartree-Fock from training
    cc=cc/np.sum(cc) 

    #random dets with replacement and probability coefs*coefs
    ids=np.random.choice(np.arange(len(cc)),size=Ns,p=cc)
    
    dets_s=[]
    for i in range(Ns) :
      dets_s.append(dec2bin(dets[ids[i]+1],N))

    return np.array(dets_s)

#********************************************************************************

N=int(sys.argv[1])
Ne=int(sys.argv[2])

tol=float(sys.argv[3])
tol=10**(-1.0*tol)
Nsample=int(sys.argv[4])


n_det_o=int(np.loadtxt('n_det'))
qp_det=np.loadtxt('dets') #[[127,127],[..,..], ..]
c=np.loadtxt('coefs')

# Pruning
print("-------------------------------------------------------------------------")
print('number of determinants before pruning : ', len(c))
print("Pruning of the wavefunction for coefficients c with |c|^2 < ",tol)
for i in range(len(c)) : 
  if c[i]**2<tol : 
    break
c=c[:i]
qp_det=qp_det[:i,:]

# Transform QP dets to CIgen dets
qp_det=qp_det.reshape((i*2))
D=np.array([QP2CIgen(int(qp_det[i]),int(qp_det[i+1]),N) for i in np.arange(0,len(qp_det),2)])


print('number of determinants after pruning : ', len(D))
np.savetxt('dets',D,fmt='%i')
np.savetxt('coefs',c,fmt='%24.15E')

# Create training sample
dets_s=create_training_sample(N,Nsample,D,c)
np.savetxt('dets_sample',dets_s,fmt='%i')
