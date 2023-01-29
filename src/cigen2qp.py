#!/bin/env python

##### CIgen ----> QP #####


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

#***********************************************************************

N=int(sys.argv[1])
Ne=int(sys.argv[2])

dets=np.loadtxt('dets',dtype = np.int64)

# Transform CIgen dets to QP dets
D=[CIgen2QP(c,N) for c in dets]

dd=np.array(D)
dd=dd.reshape((int(len(D)),2))


# QP can produce few doublons when n_det ~ 1e6
# this supress those or QP crashes
d_=np.unique(dd,axis=0)
d_=d_.reshape(np.shape(d_)[0]*2)
np.savetxt('dets', d_, fmt='% 4d')

# QP has been compiled with a max n_det of 2e7
# nd is the new n_det size after CIgen has generated new dets
nd=min(int(d_.size/2),int(2e7))

# The new determinants get zero coefficient value before diagonalization
c=np.loadtxt('coefs') 
if nd > len(c) :
 zero=[0.0 for i in range(nd-len(c))]
 COEF=list(c)+zero
else : 
 COEF=list(c)


np.savetxt('coefs', COEF, fmt='%24.15E')
np.savetxt('n_det',[nd],fmt='% 4d')
