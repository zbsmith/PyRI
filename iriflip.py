#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 19:01:58 2019

@author: zsmith
"""

import os, sys
import numpy as np
import scipy
import spacepy as sp
import geopack.geopack as gp
import os.path
import datetime
#from datetime import parser
import string
import math
import re

def FACSR(UVFAC,F107,F107A):
    '''The Schumann-Runge factors are scaled according to F10.7
    from Torr et al. GRL 1980 p6063'''
    
    I = int(0)
    LSR = int(0)
    UVFAC = np.array(60, dtype = np.float64)
    
    #............. Schumann-Runge scaling
    SRFLUX = np.array([2.4,1.4,.63,.44,.33,.17,.12,.053], dtype=np.float64)
    #...... first two SRA and SRB values out of order in Marsha's paper
    SRA = np.array([25.5,20.7,13.2,11.6,11.3,7.86,7.68,4.56], dtype = np.float64)
    SRB = np.array([222.,129.,53.4,36.0,25.0,11.3,6.35,2.05], dtype = np.float64)
    #----  Test to see if need to scale - see DATRD2 subroutine      
    #  !IF(NINT(UVFAC(58)).EQ.-1.OR.NINT(UVFAC(58)).EQ.-3) THEN
    
    for I in range(38,50+1):
        LSR=I-37
        UVFAC[I]=1.0
        if(LSR <= 8):
            UVFAC[I]=(SRA[LSR]*1.0E7*F107+SRB[LSR]*1.0E9)/SRFLUX[LSR]/1.0E11
    return

def FACEUV(UVFAC,F107,F107A):
    '''This routine uses the EUV scaling from Richards et al.[1994]
       The EUVAC flux model is based on the F74113 solar reference
       spectrum and Hinteregger's scaling factors. This subroutine
       just provides the scaling factors as a function of the proxy
       (F107+F107A)/2'''
      
    I = int(0)
    UVFAC = np.array(59)
    HFG200 = np.array(37)

    HFG200 = np.array([2.202,1.855,2.605,3.334,1.333,17.522,4.176,4.0,1.4,
                       3.694,1.791,5.385,1.889,1.899,3.427,2.051,1.392,1.619,
                       1.439,2.941,1.399,2.416,1.512,1.365,1.570,1.462,2.537,
                       1.393,1.572,1.578,1.681,1.598,1.473,1.530,1.622,1.634,
                       1.525], dtype = np.float64)
    #Test to see if need to scale - see DATRD2 subroutine
    if ((np.rint(UVFAC[58]) == -1) or (np.rint(UVFAC[58]) == -3)):
        #... EUV scaling
        F107AV = (F107+F107A)*0.5
        for I in range(1,37+1):
            A=(HFG200[I]-1)/120.0
            B=1-A*80.0
            UVFAC[I]=A*F107AV+B
            if (UVFAC[I] < 0.8): UVFAC[I]=0.8
# 50      CONTINUE
#      ENDIF
    return

def SCHUMN(J,Z,ZO2,COLUMN,SCHUPR,SCHUHT):
    '''production of o(1d) by schumann-runge bands.  The fluxes are from 
    Torr et al. GRL 1980 p6063. Scaling is done using UVFAC which may be 
    set according to F10.7 cm flux may be done in FACEUV'''
    
    global UVFAC,EUV

    SRFLUX = np.array([2.4,1.4,.63,.44,.33,.17,.12,.053], dtype = np.float64)
    SRXS = np.array([.5,1.5,3.4,6,10,13,15,12], dtype = np.float64)
    SRLAM = np.array([1725,1675,1625,1575,1525,1475,1425,1375], 
                     dtype = np.float64)
    #lmax=# of lambdas in sub. primpr: schuht=heating: schupr=o(1d) prod
    LMAX=37

    for LSR in range(1,8+1):
        #... photoabsorption cross section
        SRXSCT=1.0E-18*SRXS[LSR]
        HSRX=SRXSCT*COLUMN[2]
        if (HSRX > 70): HSRX=70
        #.. attentuated solar flux
        FLD=UVFAC[LMAX+LSR]*1.0E+11*SRFLUX[LSR]*np.exp(-HSRX)
        #... neutral heating SCHUHT and photodissociation rate SCHUPR
        SCHUHT=SCHUHT+1.24E+4*(FLD*SRXSCT)*ZO2/SRLAM[LSR]
        SCHUPR=SCHUPR+FLD*SRXSCT
        #!IF(JTI.EQ.0) WRITE(24,90) LSR,SRXSCT,FLD,SCHUPR,COLUMN(2),FLD,
        #!>   UVFAC(LMAX+LSR)
    #505  CONTINUE
    
    SCHUPR=ZO2*SCHUPR
    #90   FORMAT(2X,I5,1P,9E9.1)
#    JTI=1
    return

def PROBO2(ISW,L,ZLAM,PROB,JPTS):
    '''... o2 branching ratios are taken from kirby et al table d
       ... columns 4 & 9 are combined. columns 5,6,7&8 are combined'''
      
    A = np.array(5, dtype = np.float64)
    B = np.array(5, dtype = np.float64)
    PROB = np.array((3,6,37), dtype = np.float64)
    IPTS = 20
    X = np.array([0,20,304.,323.,454.,461.,504.,537.,556.,573.,584.,598.
                  ,610.,637.,645.,662.,684.,704.,720.,737.,774.,1026.], dtype = np.float64)
    Y = np.array([[.365,.374,.432,.435,.384,.345,.356,.365,.306,.23,.235,.245,.34,.27,.482,.675,.565,.565,1.,1.]
                  [.205,.21,.243,.245,.27,.29,.23,.27,.33,.295,.385,.35,.305,.385,.518,.325,.435,.435,0.,0.]
                  [.125,.124,.12,.12,.126,.13,.225,.216,.21,.375,.305,.37,.33,.345,0.,0.,0.,0.,0.,0.]
                  [.055,.167,.11,.105,.194,.234,.189,.149,.155,.103,.075,.036,.025,0.,0.,0.,0.,0.,0.,0.]
                  [.25,.125,.095,.95,.026,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]], dtype = np.float64)
    
    #... if zlam is too big set equal to x(max)
    #... if zlam is outside range of data values set equal to max or min
    if (ZLAM > X[20]):
        YLAM = X[20]
    elif (ZLAM <= X[1]):
        YLAM=X[1]+1.0E-3
    else:
        YLAM=ZLAM

#       DO 10 I=1,IPTS
#C kjh 6/22/92   NOTE:  I realize the following statement is strange
#C   looking, but its purpose is to prevent the CRAY compiler from
#C   vectorizing this loop.  (Which it does incorrectly).
#      if(i.eq.25)write(6,*)' '
#      IF(YLAM.GT.X(I).AND.YLAM.LE.X(I+1))  GO TO 20
# 10   CONTINUE
    I = 1
    while (((YLAM > X[I]) and (YLAM <= X[I+1])) and (I < IPTS)):
        I += 1
    
    #20   SUM=0.0
    SUM = 0.0
    
    for J in range(1,JPTS+1):
        A[J] = (Y[I+1,J]-Y[I,J])/(X[I+1]-X[I])
        B[J] = Y[I,J]-A[J]*X[I]
        SUM += A[J]*YLAM+B[J]
    #30   CONTINUE
    
    for J in range(1,JPTS+1):
        PROB[2,J,L] = (A[J]*YLAM+B[J])/SUM
    #40   CONTINUE
    
    return

def PROBN2(ISW,L,ZLAM,PROB,JPTS):
    '''the n2 probabilities are taken from kirby et al tables b and c
    the yield of n+ is determined first then the remaining portion
    of the cross section is distributed amongst the n2+ ion states
    (x,a,b). the dissociation yield is divided between the 3 higher
    energy states according to wight et al. j.phys. b, 1976
    the 2 other states of kirby et al are not included'''
    
    A = np.array(6, dtype = np.float64)
    B = np.array(6, dtype = np.float64)
#    PROB(3,6,37),SUM,YIELD,YLAM,ZLAM
    IPTS = 14
    X = np.arrau([0.,50.,210.,240.,280.,300.,332.,428.,500.,600.,660.,660.01,
                  720.,747.,796.], dtype = np.float64)
    Y = np.arrau([[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
                  [0,.32,.32,.32,.32,.32,.3,.46,.404,.308,.308,.308,.42,1.,1.],
                  [0,.55,.55,.55,.55,.55,.52,.46,.506,.589,.589,.692,.58,.0,.0],
                  [0,.13,.13,.13,.13,.13,.12,.08,.09,.103,.103,0.0,0.0,0.0,0.0],
                  [0,0.0,0.0,0.0,.05,.1,.15,.83,1.,.0,.0,.0,.0,.0,.0],
                  [0,.0,.0,.0,.3,.4,.79,.17,.0,.0,.0,.0,.0,.0,.0],
                  [0,1.,1.,1.,.65,.5,.06,.0,.0,.0,.0,.0,.0,.0,.0]],dtype = np.float64)
    
    #if zlam is too big set equal to x(max)
    YLAM=ZLAM
    #Prevent divide by zero
    if (ZLAM > X[14]): YLAM=X[14]-1
    if (ZLAM < X[1]): YLAM=X[1]+1
    
    YIELD=0.0
    #determine yield of n+, and store in prob array
    YLDISS(1,YLAM,YIELD)
    
#    DO 10 I=1,IPTS
#C kjh 6/22/92   NOTE:  I realize the following statement is strange
#C   looking, but its purpose is to prevent the CRAY compiler from
#C   vectorizing this loop.  (Which it does incorrectly).
#      if(i.eq.25)write(6,*)' '
#      IF(YLAM.GT.X(I).AND.YLAM.LE.X(I+1))  GO TO 20
# 10   CONTINUE
    I = 0
    while ((YLAM > X[I]) and (YLAM <= X[I+1])) and (I <= IPTS):
        I += 1

# 20   SUM=0.0
    SUM=0.0
    #fit straight line between points
    for J in range(1,JPTS+1):
        A[J]=(Y[I+1,J]-Y[I,J])/(X[I+1]-X[I])
        B[J]=Y[I,J]-A[J]*X[I]
    #30   CONTINUE
    #determine probabilities of n2+ states
    for J in range(1,JPTS+1):
        if (J <= 3): PROB[3,J,L]=(A[J]*YLAM+B[J])*(1-YIELD)
        if (J > 3): PROB[3,J,L]=(A[J]*YLAM+B[J])*YIELD
        SUM += PROB[3,J,L]
    #40   CONTINUE
    
    if (SUM == 0.0): return
    #normalise probabilities
    PROB[3,:,:] /= SUM
#      DO 50 J=1,JPTS
# 50   PROB(3,J,L)=PROB(3,J,L)/SUM
    return

def YLDISS(ISW,ZLAM,YIELD):
    '''determination of dissociative yield of n+, refer to kirby et al
    page 66 and table b'''
    
    IPTS = 9
    X = np.array([0,50.,210.,240.,302.,387.,477.,496.,509.,2000.])
    Y = np.array([0,.36,.36,.346,.202,.033,.041,.024,0.0,0.0])

#       DO 10 I=1,IPTS
#C kjh 6/22/92   NOTE:  I realize the following statement is strange
#C   looking, but its purpose is to prevent the CRAY compiler from
#C   vectorizing this loop.  (Which it does incorrectly).
#      if(i.eq.25)write(6,*)' '
#      IF(ZLAM.GE.X(I).AND.ZLAM.LT.X(I+1))  GO TO 20
# 10   CONTINUE
    I = 1
    while (((ZLAM >= X[I]) and (ZLAM < X[I+1])) and (I < IPTS)):
        I += 1
        
    #20   IF(ZLAM.GT.387.AND.ZLAM.LT.477) GO TO 40
    if (ZLAM > 387) and (ZLAM < 477):
    #... parabolic interpolation see formula page 66 kirby et al
    #40   YIELD=.0329+8.13E-6*(ZLAM-442)**2
        YIELD=.0329+8.13E-6*(ZLAM-442)**2
    else:
    #... linear interpolation
        YIELD=(ZLAM-X[I])/(X[I+1]-X[I])*(Y[I+1]-Y[I])+Y[I]
#      GO TO 30
# 30   RETURN
    return (YIELD)

def PROBS(ISW,PROB,ZLAM,LMAX,NNI):
    '''program for finding branching ratios (probabilities for various ion
    and molecular states) of o,o2,n2
    ---refs--- m. torr et al grl 1979 page 771, kirby et al atomic data
    and nuclear tables 1979 23,page 63'''
    
    NNI = np.array(3, dtype = np.int16)
    PROB = np.array((3,6,37), dtype = np.float64)
    ZLAM = np.array(37, dtype = np.float64)
    #coefficients of o ionization cross sections from torr et al table 2
    YO = np.array([[.19,.486,.952,1.311,1.539,1.77,1.628,1.92,1.925,2.259,2.559,2.523,3.073,3.34,3.394,3.421,3.65,3.92,3.62,3.61,3.88,4.25,5.128,4.89,6.739,4.0,3.89,3.749,5.091,3.498,4.554,1.315,0.0,0.0,0.0,0.0,0.0],
                   [.206,.529,1.171,1.762,2.138,2.62,2.325,2.842,2.849,3.446,3.936,3.883,4.896,5.37,5.459,5.427,5.67,6.02,5.91,6.17,6.29,6.159,11.453,6.57,3.997,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
                   [.134,.345,.768,1.144,1.363,1.63,1.488,1.92,1.925,2.173,2.558,2.422,2.986,3.22,3.274,3.211,3.27,3.15,3.494,3.62,3.23,2.956,0.664,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
                   [.062,.163,.348,.508,.598,.71,.637,.691,.693,.815,.787,.859,.541,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
                   [.049,.13,.278,.366,.412,.35,.383,.307,.308,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]],dtype = np.float64)

    #....... production of o states from torr et al table 2 (yo array)
    #....... need to reverse order of yo to correspond with lambda
    for L in range(1,LMAX+1):
        LL=LMAX+1-L
        SUM=YO[1,LL]+YO[2,LL]+YO[3,LL]+YO[4,LL]+YO[5,LL]
        for I in range(1,5+1):
            PROB[1,I,L]=0.0
            if (SUM != 0.0): PROB[1,I,L]=YO[I,LL]/SUM
    
    #call separate subroutines for o2 and n2 probabilities
    for L in range(1,LMAX+1):
        PROBO2(1,L,ZLAM[L],PROB,NNI[2])
        PROBN2(1,L,ZLAM[L],PROB,NNI[3])
    
    if (ISW == 0): return
#      WRITE(17,95)
# 95   FORMAT(/5X,' Photoionization branching ratios for O, O2, N2'
#     > ,/3X,'Lam    4S   2D   2P   4P   2P*   -   X2   a+A  b4   B2 '
#     > ,'  dis   -  X2   A2   B2   C2   F2   2s')
#      DO 50 L=1,LMAX
# 50   WRITE(17,90) ZLAM(L),((PROB(IS,J,L),J=1,6),IS=1,3)
# 90   FORMAT(F8.2,22F5.2)
    return

def PARAMS(ISW,LMAX):
    '''this program determines the cross sections, solar fluxes, and
    given in m. torr et al g.r.l 1979 p771, table 2 and 3. but
    the longer wavelengths come first'''
    
#      Integer, Intent(In) :: ISW
#      Integer, Intent(INOUT) :: LMAX

#      INTEGER  I,IN,IS,J,L,LAMAX,NNI
#      REAL EUV,FFAC, SIGABS,SIGION,TPOT,UVFAC,ZFLUX,ZLAM
    global ZFLUX,SIGABS,ZLAM,SIGION,TPOT,NNI,LAMAX,UVFAC,EUV
    
    SIGABS = np.array((3,37), dtype=np.float64)
    ZLAM = np.array(37, dtype=np.float64)
    SIGION = np.array((3,37), dtype=np.float64)
    TPOT = np.array((3,10), dtype=np.float64)
    NNI = np.array(3, dtype=np.float64)
    
    #ionization potentials for o,o2 and n2 see kirby et al note the
    #o2 2(pi)u and 2(sigma-)u , and n2 f2(sigma)u pots are guesses
    #the sixth n2 potential is for dissociation
    X3 = np.array([13.6,16.9,18.6,28.5,40.,0.0,12.1,16.5,18.2,20.,25.,0.,15.6,
                   16.7,18.8,25.,29.,37.], dtype=np.float64)
    #wavelength data. average is taken for bands
    ZLX = np.array([1025.,1031.91,1025.72,975.,977.02,925.,875.,825.,775.,
                    789.36,770.41,765.15,725.,703.36,675.,625.,629.73,609.76,
                    575.,584.33,554.31,525.,475.,465.22,425.,375.,368.07,325.,
                    303.78,303.31,275.,284.15,256.3,225.,175.,125.,75.], dtype=np.float64)
    #fluxes from table 3. these are for 74113. just replace this data
    #for other years in the table. note!!!! flux doubled for lambda<250
    #shortest wavelenghts have been tripled
    ZFX = np.array([2.4665,2.1,3.5,1.4746,4.4,3.,3.537,1.625,.758,.702,.26,
                    .17,.141,.36,.23,.342,1.59,.53,.357,1.27,.72,.452,.285,
                    .29,.383,.314,.65,.965,6.9,.8,1.679,.21,.46,3.1,4.8,.45,
                    1.2], dtype=np.float64)
    #absorption cross sections -- o first ,o2, then n2
    X1 = np.array([0.0,0.0,0.0,0.0,0.0,1.315,4.554,3.498,5.091,3.749,3.89,4,
                   10.736,11.46,17.245,13.365,13.4,13.4,13.024,13.09,12.59,
                   12.059,12.127,11.93,11.496,9.687,9.84,8.693,7.7,7.68,6.461,
                   7.08,6.05,5.202,3.732,1.839,.73,1.346,1.0,1.63,21.108,
                   18.73,12.817,8.562,16.631,22.145,26.668,18.91,20.8,28.535,
                   27.44,21.919,26.017,32.06,28.07,26.61,22.79,26.04,24.606,
                   23.101,21.91,20.31,18.118,18.32,17.438,16.81,16.8,14.387,
                   15.79,13.37,10.9,7.509,3.806,1.316,0.0,0.0,0.0,50.988,2.24,
                   9.68,20.249,16.992,33.578,16.487,14.18,120.49,24.662,26.54,
                   31.755,23.339,23.37,22.79,22.787,22.4,24.13,24.501,23.471,
                   23.16,21.675,16.395,16.91,13.857,11.7,11.67,10.493,10.9,
                   10.21,8.392,4.958,2.261,0.72], dtype=np.float64)
    #ionization cross sections 
    X2 = np.array([0.0,0.0,0.0,0.0,0.0,1.315,4.554,3.498,5.091,3.749,3.89,4,
                   10.736,11.46,17.245,13.365,13.4,13.4,13.024,13.09,12.59,
                   12.059,12.127,11.93,11.496,9.687,9.84,8.693,7.7,7.68,6.461,
                   7.08,6.05,5.202,3.732,1.839,.73,.259,0.0,1.05,13.94,15.54,
                   9.374,5.494,6.413,10.597,10.191,8.47,11.72,23.805,23.75,
                   21.306,24.937,31.1,26.39,26.61,22.79,26.04,24.606,23.101,
                   21.91,20.31,18.118,18.32,17.438,16.81,16.8,14.387,15.79,
                   13.37,10.9,7.509,3.806,1.316,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                   0.0,14.274,8.86,8.5,65.8,15.06,25.48,29.235,23.339,23.37,
                   22.79,22.787,22.4,24.13,24.501,23.471,23.16,21.675,16.395,
                   16.91,13.857,11.7,11.67,10.493,10.9,10.21,8.392,4.958,
                   2.261,0.72], dtype=np.float64)

    NNI[1]=5
    NNI[2]=5
    NNI[3]=6
    LMAX=37
#      IF(ISW.NE.0) WRITE(17,95)
# 95   FORMAT(/5X,'EUV fluxes, Photoabsorption, and Photoionization ',
#     >  'Cross sections',
#     > /4X,'I',5X,'lam',5X,'flux',4X,'sigaOX',3X,'sigaO2'
#     > ,3X,'sigaN2',3X,'sigiOX',3X,'sigiO2',3X,'sigiN2',3X,'UVfac')

    for L in range(1,LMAX+1):
        ZLAM[L]=ZLX[L]
        FFAC=UVFAC[LMAX+1-L]
        if (ZFX[L] < 100): ZFLUX[L]=ZFX[L]*1.0E+9*FFAC
        #..- setting up ionization potentials
        if (L <= 6):
            TPOT[1,L] = X3[L]
            TPOT[2,L] = X3[6+L]
            TPOT[3,L] = X3[12+L]
        #..- setting up cross sections
        for IS in range(1,3+1):
            IN=LMAX*(IS-1)+L
            SIGABS[IS,L] = X1[IN]*1.0E-18
            SIGION[IS,L]=X2[IN]*1.0E-18
            if (SIGABS[IS,L] < SIGION[IS,L]): SIGABS[IS,L] = SIGION[IS,L]
        #10   CONTINUE
        
#      IF(ISW.EQ.0) GO TO 20
#      WRITE(17,90) L,ZLAM(L),ZFLUX(L),(SIGABS(I,L),I=1,3)
#     > ,(SIGION(I,L),I=1,3),FFAC
    #20   CONTINUE
    
    #IF(ISW.EQ.0) RETURN
#      WRITE(17,94)
# 94   FORMAT(/5X,' Ionization potentials for O, O2, N2'
#     > ,/2X,'4S   2D   2P   4P   2P*  -   X2   a+A  b4   B2   dis  -'
#     > ,'  X2   A2   B2   C2   F2   2s')
# 60   WRITE(17,91) ((TPOT(I,J),J=1,6),I=1,3)

    return
# 90   FORMAT(1X,I4,F9.2,1P,22E9.2)
# 91   FORMAT(22F5.1)

