# -*- coding: latin-1 -*-

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

# iridreg.for, version number can be found at the end of this comment.
#-----------------------------------------------------------------------
#
# This file contains the D-region models of Friedrich and Torkar (2001)
# (subroutine F00 and block data statement).
# The subroutine DRegion of Danilov et al. (1995) was moved to IRIFUN],
# because of consistent problems of some Fortran compilers wit the long
# BLOCK DATA statement. 
#
# !!!USER NOTE!!! If your compiler has problems with this subroutine you 
# can compile IRI without this file. But you first have to comment out  
# the following two line in IRISUB: 
#            call F00(HEIGHT,LATI,DAYNR,XHI,F107D,EDENS,IERROR)
#            if(ierror.eq.0.or.ierror.eq.2) outf(1,kk)=edens
#
#-----------------------------------------------------------------------
# Corrections/Version Numbers:
#-Version-mm/dd/yy-description (person reporting correction)
# 2001.01 05/07/01 initial version
# 2001.02 07/11/01 new version of F00 (as provided by K. Torkar)
# 2002.01 28/10/02 replace TAB/6 blanks, PARAMETER () (D. Simpson)
# 2007.00 05/18/07 Release of IRI-2007
# 2012.00 12/29/11 Release of IRI-2012; no change in iridreg.for
# 2012.00 01/18/12 Moved subroutine DRegion (Danilov model) to IRIFUN
#-----------------------------------------------------------------------
#
#
NHGT=81
NLAT=5
NMON=12
NZEN=12
NF107=3
#
global EDEN,TABHE,TABLA,TABMO,TABZA,TABFL
EDEN = np.full((NHGT+1,NLAT+1,NMON+1,NZEN+1,NF107+1),0)

TABHE = np.full(NHGT+1,0)
TABLA = np.full(NLAT+1,0)
TABMO = np.full(NMON+1,0)
TABZA = np.full(NZEN+1,0)
TABFL = np.full(NF107+1,0)
#
#      INTEGER I,L
#
#     altitudes in km
TABHE = [0,60.,61.,62.,63.,64.,65.,66.,67.,68.,69.,
         70.,71.,72.,73.,74.,75.,76.,77.,78.,79.,
         80.,81.,82.,83.,84.,85.,86.,87.,88.,89.,
         90.,91.,92.,93.,94.,95.,96.,97.,98.,99.,
         100.,101.,102.,103.,104.,105.,106.,107.,108.,109.,
         110.,111.,112.,113.,114.,115.,116.,117.,118.,119.,
         120.,121.,122.,123.,124.,125.,126.,127.,128.,129.,
         130.,131.,132.,133.,134.,135.,136.,137.,138.,139.,140.]
#
#     latitudes in degree
TABLA = [0,0.,15.,30.,45.,60.]
#
#     months
TABMO = [0,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.]
#
#     solar zenith angles in degree
TABZA = [0,0.,30.,45.,60.,75.,80.,85.,90.,95.,100.,130.,180.]
#
#     log10(F10.7) for 75,130,200 Jy
TABFL = [0,1.87506, 2.11394, 2.30103]
#
#     log10 electron densities, ordered as
#     I,J,K,L,M = Height,Latitude,Month,Zenithangle,F10.7
#     8 heights in each line
#     12 zenith angles in each DATA statement
#     innermost loop: J (5 latitudes)
#     next loop:      K (12 months)
#     next loop:      M (3 F10.7-fluxes)
#     outermost loop: I (11 groups of heights)
                
def F00(HGT,GLAT1,IDAY,ZANG,F107T,EDENS,IERROR):
    '''---------------------------------------------------------------------
        PURPOSE:
            THIS SUBROUTINE COMPUTES "FIRI" ELECTRON DENSITIES
            
            COMMON BLOCK REQUIRED:
                REAL EDEN,TABHE,TABLA,TABMO,TABZA,TABFL
                COMMON/FIRCOM/EDEN(81,5,12,12,3)],
                1              TABHE(81),TABLA(5),TABMO(12),TABZA(12),TABFL(3)
                
            ARRAY EDEN contains LOG10(tabulated electron density],
            ordered in (height,latitude,month,zenithangle,f107)
            Quantity      Minimum  Maximum  Number of steps
            Height        60       140      81
            Latitude       0        60       5
            Month          1        12      12
            Zenith angle   0       180      12
            F10.7         75       200       3
            
        PARAMETERS:
            HGT   height in km (input, REAL)
            GLAT1 latitude in degrees, north is positive (input, REAL)
            IDAY  day of the year (input, INTEGER)
            ZANG  solar zenith angle in degrees (input, REAL)
            F107T 10.7cm flux in Ja (input, REAL)
            EDENS model electron density in m**-3 (output, REAL)
            IERROR  Error code (INTEGER, output)
            
        Error code
         0         no error
         1         model undefined for given combination of input
                   parameters, output is set to zero
         2         input parameters outside valid range, output is invalid
         3         both error conditions detected, output is zero
         
        USAGE
          CALL F00(HGT,GLAT1,IDAY,ZANG,F107T,EDENS,IERROR)
          
        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
            none
            
        Reference: Friedrich, M., Torkar, K. FIRI: a semiempirical model of the
                    lower ionosphere. 
                   J. Geophys. Res. 106 (A10), 21409Ð21418, 2001.
        WRITTEN BY K. TORKAR, IWF GRAZ
        Klaus.Torkar@oeaw.ac.at
        
        LAST MODIFICATION:  06.07.2001
        
        VERSION: 1.1
        
       ------------------------------------------------------------------'''
#
#      REAL HGT,GLAT1,ZANG,F107T,EDENS,F107L
#      INTEGER IDAY,IERROR
#
    NHGT=81
    NLAT=5
    NMON=12
    NZEN=12
    NF107=3
#
#      REAL EDEN,TABHE,TABLA,TABMO,TABZA,TABFL
    global EDEN,TABHE,TABLA,TABMO,TABZA,TABFL
#      COMMON/FIRCOM/EDEN(81,5,12,12,3)],
#     1              TABHE(81),TABLA(5),TABMO(12),TABZA(12),TABFL(3)
#
#      INTEGER MON,I,J,L,M,ISTEPJ,I1,I2,J1,J2,K1,K2,L1,L2,M1,M2
#      INTEGER TABM(12)
    TABM = np.zeros(12)
    EDENI = np.zeros([2,2,2,2])
    EDENIJ = np.zeros([2,2,2])
    EDENIJK = np.zeros([2,2])
    EDENIJKL = np.zeros(2)
#      REAL STEPJ,DAY1,H1,DEG1,XHI1,FLX1,EL
#
    TABM = np.array([0,31,59,90,120,151,181,212,243,273,304,334])
    STEPJ = 15.0
    ISTEPJ = 15
#
#     INDICES:
#     I=HEIGHT, J=LATITUDE, K=MONTH, L=ZANG, M=F10.7
#
#     CHECK INPUT
#
    IERROR=0
    F107L=np.log10(np.min(1000.0,np.max(1.0,F107T)))
    if (((((((((HGT < TABHE[1]) or (HGT > TABHE[NHGT])) or (GLAT1 > TABLA[NLAT])) or (GLAT1 < -TABLA[NLAT])) or (IDAY < 1)) or (IDAY > 366)) or (ZANG < TABZA[1])) or (ZANG > TABZA[NZEN])) or (F107L < TABFL[1])) or (F107L > TABFL[NF107]): 
        IERROR=2
#
#     assume height table is in 1 km steps from 60 to 140 km
    I=np.min(NHGT-1,int(HGT)-59)
    if (I < 1): I=1
    H1=HGT-TABHE[I]
    I1=I
    I2=I+1
#
#     assume latitude table is in 15 deg steps from 0 to 60 deg
    J=np.max(1,np.min(NLAT-1,int(np.abs(GLAT1))/ISTEPJ))
    DEG1=(np.abs(GLAT1)-TABLA[J])/STEPJ
    J1=J
    J2=J+1
#
#     assume month table is given for each month
    MON=12
    while (TABM[MON] > IDAY):
        MON=MON-1
    
    DAY1=float(IDAY-TABM[MON]-15)/30.0
    if (DAY1 < 0.0): MON=MON-1
    if (MON > 1) and (MON < 11):
        K1=MON
        K2=MON+1
    else:
        K1=12
        K2=1
#
#     assume zenith angle table has 12 entries between 0 and 180 deg
    L = 2
    while (ZANG >= TABZA[L]) and (L < NZEN):
        L += 1  #loop without break
#    for L in range(2,NZEN-1):
#        if (ZANG < TABZA[L]): break #GOTO 1
    if (ZANG >= TABZA[L]): L=NZEN
#1     L=L-1
    L=L-1
    L1=L
    L2=L+1
    XHI1=(ZANG-TABZA[L1])/(TABZA[L2]-TABZA[L1])
#
#     assume solar activity table has 3 entries
    F107L=np.min(TABFL[3],np.max(TABFL[1],F107L))
    if (F107L < TABFL[NF107-1]):
        M1=1
        M2=2
    else:
        M1=2
        M2=3
    FLX1=(F107L-TABFL[M1])/(TABFL[M2]-TABFL[M1])
#
#     ADJUST SOUTHERN LATITUDES TO NORTH AND MONTH+6
#
    if (GLAT1 < 0.0):
        K1=K1+6
        if (K1 > 12): K1=K1-12
        K2=K2+6
        if (K2 > 12): K2=K2-12
#
#     EDEN(hgt,lat,mon,zang,f107)
#          I   J   K   L    M
#
    for M in range(M1,M2):
#
        MH=M+1-M1
#       INTERPOLATE IN HEIGHT I
        for L in range(L1,L2):
            if (EDEN[I1,J1,K1,L,M] == 0.0) or (EDEN[I2,J1,K1,L,M] == 0.0) or (EDEN[I1,J2,K1,L,M] == 0.0) or (EDEN[I2,J2,K1,L,M] == 0.0) or (EDEN[I1,J1,K2,L,M] == 0.0) or (EDEN[I2,J1,K2,L,M] == 0.0) or (EDEN[I1,J2,K2,L,M] == 0.0) or (EDEN[I2,J2,K2,L,M] == 0.0):
                 EDENS=0.0
                 IERROR=IERROR+1
                 return
          
            if (HGT < TABHE[1]):
                EDENI[1,1,L+1-L1,MH] = EDEN[I1,J1,K1,L,M]
                EDENI[2,1,L+1-L1,MH] = EDEN[I1,J2,K1,L,M]
                EDENI[1,2,L+1-L1,MH] = EDEN[I1,J1,K2,L,M]
                EDENI[2,2,L+1-L1,MH] = EDEN[I1,J2,K2,L,M]
            elif (HGT > TABHE[NHGT]):
                EDENI[1,1,L+1-L1,MH] = EDEN[I2,J1,K1,L,M]
                EDENI[2,1,L+1-L1,MH] = EDEN[I2,J2,K1,L,M]
                EDENI[1,2,L+1-L1,MH] = EDEN[I2,J1,K2,L,M]
                EDENI[2,2,L+1-L1,MH] = EDEN[I2,J2,K2,L,M]
            else:
                EDENI[1,1,L+1-L1,MH] = EDEN[I1,J1,K1,L,M] + H1*(EDEN[I2,J1,K1,L,M] - EDEN[I1,J1,K1,L,M])
                EDENI[2,1,L+1-L1,MH] = EDEN[I1,J2,K1,L,M] + H1*(EDEN[I2,J2,K1,L,M] - EDEN[I1,J2,K1,L,M])
                EDENI[1,2,L+1-L1,MH] = EDEN[I1,J1,K2,L,M] + H1*(EDEN[I2,J1,K2,L,M] - EDEN[I1,J1,K2,L,M])
                EDENI[2,2,L+1-L1,MH] = EDEN[I1,J2,K2,L,M] + H1*(EDEN[I2,J2,K2,L,M] - EDEN[I1,J2,K2,L,M])        
#
#       INTERPOLATE IN LATITUDE J
        for L in range(1,2):
            if (np.abs(GLAT1) > TABLA[NLAT]):
                EDENIJ[1,L,MH] = EDENI[2,1,L,MH]
                EDENIJ[2,L,MH] = EDENI[2,2,L,MH]
            else:
                EDENIJ[1,L,MH] = EDENI[1,1,L,MH] + DEG1*(EDENI[2,1,L,MH] - EDENI[1,1,L,MH])
                EDENIJ[2,L,MH] = EDENI[1,2,L,MH] + DEG1*(EDENI[2,2,L,MH] - EDENI[1,2,L,MH])
#
#       INTERPOLATE IN MONTH K
        EDENIJK[1,MH] = EDENIJ[1,1,MH] + DAY1*(EDENIJ[2,1,MH] - EDENIJ[1,1,MH])
        EDENIJK[2,MH] = EDENIJ[1,2,MH] + DAY1*(EDENIJ[2,2,MH] - EDENIJ[1,2,MH])
#
#       INTERPOLATE IN ZENITH ANGLE L
        EDENIJKL[MH]=EDENIJK[1,MH]+XHI1*(EDENIJK[2,MH]-EDENIJK[1,MH])
#
    EL=EDENIJKL[1]+FLX1*(EDENIJKL[2]-EDENIJKL[1])
#
    EDENS = math.pow(10.,EL)
#
    return (EDENS, IERROR)

def arrayIndexes(line):
    cutLine = line.lstrip('      DATA ((EDEN(').rstrip(')/')
    indices = cutLine.split(')')[0].split(',')
    return (indices)

def variables(line):
    indices = arrayIndexes(line)
    numInd = len(indices)
    variables = []
    element = []
    for i in range(0,numInd):
        temp = indices.pop()
        if temp.isalpha():
            variables.append(temp)
            element.append(i)
            del temp
    variables.reverse()
#        element.reverse()
    return (variables, element)

def valRange (line, index):
    firstInstance = line.find(index)
    linePart = line[(firstInstance+1):]
    secondInstance = linePart.find(index)
    linePart = linePart[(secondInstance+1):]
    equal = linePart.find('=')
    comma = linePart.find(',')
    numStr = linePart[(equal+1):comma]
    closePar = linePart.find(')')
    numStp = linePart[(comma+1):closePar]
    return (numStr, numStp)

#
#
#    BLOCK DATA
#
#     PURPOSE:
#     DEFINES TABLES OF FIRI(2000) IN
#     ARRAY EDEN(height,latitude,month,zenithangle,f107)
#     Quantity      Minimum  Maximum  Number of steps
#     Height        60       140      81
#     Latitude       0        60       5
#     Month          1        12      12
#     Zenith angle   0       180      12
#     F10.7         75       200       3
#
#     WRITTEN BY K. TORKAR, IWF GRAZ
#     Klaus.Torkar@oeaw.ac.at
#
#     LAST MODIFICATION:  01.09.2000
#
#     VERSION: 1.1
#
#     ------------------------------------------------------------------
#
#     EDEN(hgt,lat,mon,zang,f107)
#          I   J   K   L    M
f = open('fortranData.dat')
p = re.compile('^c', re.IGNORECASE)
line = f.readline()
while p.search(line) == None:        
    indices = arrayIndexes(line)
    baseIndex = np.ones(len(indices), dtype=int)
    for i in range(0,len(indices)):
        if indices[i].isnumeric(): baseIndex[i] = int(indices[i])
    var, ele = variables(line)
    star = np.zeros(len(var), dtype=int)
    stp = np.zeros(len(var), dtype=int)
    for ind in range(0,len(var)):
        vR = valRange(line, var[ind])
        if vR[0].isalpha():
            star[ind] = int(eval(vR[0]))
        else:
            star[ind] = int(vR[0])
        if vR[1].isalpha():
            stp[ind] = int(eval(vR[1]))
        else:
            stp[ind] = int(vR[1])
        del vR

    for u in range(int(star[1]),int(stp[1])+1):
        line = f.readline()
        line = line.strip(line[0:(line.find('*')+1)])
        for v in range(int(star[0]),int(stp[0]+1)):
            index = baseIndex
            index[ele[0]-1] = int(v)
            index[ele[1]-1] = int(u)
            EDEN[index] = np.float64(line[0:(line.find(',')-1)])

            line = line[(line.find(',')+1):]
    line = f.readline()