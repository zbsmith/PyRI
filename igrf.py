import numpy as np
#import scipy
#import spacepy as sp
#import geopack.geopack as gp
#import os.path
import datetime
#from datetime import parser
#import string
import math
#-----------------------------------------------------------------------        
#
# Subroutines to compute IGRF parameters for IRI and all functions and 
# subroutines required for this computation, including:
# 	IGRF_SUB, IGRF_DIP, FINDB0, SHELLG, STOER, FELDG, FELDCOF, GETSHC, 
# 	INTERSHC, EXTRASHC, GEODIP, fmodip
#
# CGM coordinates : GEOCGM01, OVL_ANG, CGMGLA, CGMGLO, DFR1DR, 
#   AZM_ANG, MLTUT, MFC, FTPRNT, GEOLOW, CORGEO, GEOCOR, SHAG, RIGHT, 
#   IGRF, RECALC, SPHCAR, BSPCAR, GEOMAG, MAGSM, SMGSM
#
# MLT: CLCMLT, DPMTRX
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Required i/o units:  
#  KONSOL= 6 Program messages (used when jf(12)=.true. -> konsol)
#  KONSOL=11 Program messages (used when jf(12)=.false. -> MESSAGES.TXT)
#
#     COMMON/iounit/konsol,mess is used to pass the value of KONSOL from 
#     IRISUB to IRIFUN and IGRF. If mess=false then messages are turned off.
#     
#  UNIT=14 IGRF/GETSHC: IGRF coeff. (DGRF%%%%.DAT or IGRF%%%%.DAT, %%%%=year)
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Corrections:
# 11/01/91 SHELLG: lowest starting point for B0 search is 2  
#  1/27/92 Adopted to IGRF-91 coeffcients model
#  2/05/92 Reduce variable names: INTER(P)SHC,EXTRA(P)SHC,INITI(ALI)ZE
#  8/08/95 Updated to IGRF-45-95; new coeff. DGRF90, IGRF95, IGRF95S
#  5/31/00 Updated to IGRF-45-00; new coeff.: IGRF00, IGRF00s
#-Version-mm/dd/yy-Description (Person reporting the correction)
# 2000.01 05/07/01 initial version
# 2000.02 07/11/01 replace feldi(xi,h) by feldi (P. Wilkinson)
# 2000.02 07/11/01 variables EGNR, AGNR,OGNR not used (P. Wilkinson)
# 2000.01 10/28/02 replace TAB/6 blanks, enforce 72/line (D. Simpson)
# 2000.02 11/08/02 change unit for coefficients to 14
# 2000.03 06/05/03 correct DIPL computation (V. Truhlik)
# 2005.00 04/25/05 CALL FELDI and DO 1111 I=1,7 (Alexey Petrov)
# 2005.01 11/10/05 added igrf_dip and geodip (MLAT) 
# 2005.02 11/10/05 FELDCOF: updated to IGRF-10 version
# 2005.03 12/21/06 GH2(120) -> GH2(144)
# 2007.00 05/18/07 Release of IRI-2007
# 2007.08 07/30/09 SHELLG,STOER,FELDG,FELDCOF: NMAX=13; H/G-arrays(195) 
# 2007.10 02/26/10 FELDCOF: updated to IGRF-11; DGRF05, IGRF10, IGRF10S
# 2007.11 04/27/10 RECALC: updated to IGRF-11
# 2007.11 04/27/10 Make all arrays(195) to arrays(196) 
# 2007.11 04/27/10 FELDCOF: corrected Filmod and also IGRF10.DAT
# 2007.11 04/29/10 New files dgrf%%%%.asc; new GETSHC; char*12 to 13
#
# 2012.00 10/05/11 IRI-2012: bottomside B0 B1 model (SHAMDB0D, SHAB1D),
# 2012.00 10/05/11    bottomside Ni model (iriflip.for), auroral foE
# 2012.00 10/05/11    storm model (storme_ap), Te with PF10.7 (elteik),
# 2012.00 10/05/11    oval kp model (auroral_boundary), IGRF-11(igrf.for), 
# 2012.00 10/05/11    NRLMSIS00 (cira.for), CGM coordinates, F10.7 daily
# 2012.00 10/05/11    81-day 365-day indices (apf107.dat), ap->kp (ckp),
# 2012.00 10/05/11    array size change jf(50) outf(20,1000), oarr(100).
# 2012.02 12/17/12 igrf_dip: Add magnetic declination as output parameter
# 2014.01 07/20/14 igrf_dip,FTPRNT,RECALC: ASIN(x): abs(x)>1.0 x=sign(1.,x)
# 2014.02 07/24/14 COMMON/iounit: added 'mess' 
# 2015.01 02/10/15 Updating to IGRF-12 (2015)
# 2015.01 07/12/15 use mess,konsol in IGRF and RECALC
# 2015.02 08/23/15 initialization of Earth constants moved to IRI_SUB
# 2015.03 10/14/15 CLCMLT,DPMTRX <--- IRIFUN.FOR
# 2015.03 10/14/15 RECALC: update with IGRF-12 until 2020
# 2015.03 10/14/15 IGRF_SUB,_DIP: move CALL FELDCOF to IRISUB.FOR
# 2015.03 10/14/15 FELDCOF,SHELLG: DIMO to COMMON/IGRF1/
# 2016.01 02/17/16 GEODIP: add PI to CONST
# 2017.01 07/07/17 IGRF: updated with newest 2010, 2015, 2015s coeff.
# 2019.01 09/01/19 First Python Version
#-----------------------------------------------------------------------        
# 

t0 =  datetime.datetime(1970,1,1) #T0 for UTC
MA = 0
IYR = 0
G = np.array(67)
H = np.array(67)
REC = np.array(67)

def igrf_sub(xlat,xlong,year,height,xl,icode,dipl,babs):
    '''-----------------------------------------------------------------------        
        INPUT:
                xlat      geodetic latitude in degrees
                xlong     geodetic longitude in degrees
                year      decimal year (ycear+(month-0.5)/12.0-0.5 or
                          year+day-of-year/365 or ../366 if leap year)
                height    height in km
        OUTPUT:
                xl        L value
                icode      =1  L is correct; 
                           =2  L is not correct;
                           =3  an approximation is used
                dipl      dip latitude in degrees
                babs      magnetic field strength in Gauss
       ----------------------------------------------------------------------'''
    global UMR,PI
    
    bnorth = 0.0
    beast = 0.0
    bdown = 0.0
    bab1 = 0.0
    lati=xlat
    longi=xlong
    #      CALL FELDCOF(YEAR,DIMO)
    FELDG(lati,longi,height,bnorth,beast,bdown,babs)
    SHELLG(lati,longi,height,xl,icode,bab1)
    dipl = math.atan(bdown/2.0/np.sqrt(bnorth*bnorth+beast*beast))/UMR
    return
#
#
def igrf_dip(xlat,xlong,year,height,dec,dip,dipl,ymodip):
    '''-----------------------------------------------------------------------
        INPUT:
                xlat      geodetic latitude in degrees
                xlong     geodetic longitude in degrees
                year      decimal year (year+month/12.0-0.5 or
                          year+day-of-year/365 or ../366 if leap year)
                height    height in km
        OUTPUT:
                dec       magnetic declination in degrees
                dip       magnetic inclination (dip) in degrees
                dipl      dip latitude in degrees
                ymodip    modified dip latitude = asin{dip/sqrt[dip^2+cos(LATI)]}
        -----------------------------------------------------------------------'''
    global UMR,PI
    
    bnorth = 0.0
    beast = 0.0
    bdown = 0.0
    babs = 0.0
    xlati = xlat
    xlongi = xlong
    h = height
    #      CALL FELDCOF(YEAR,DIMO)
    FELDG(xlati,xlongi,h,bnorth,beast,bdown,babs)
    decarg=beast/np.sqrt(beast*beast+bnorth*bnorth)
    if (np.abs(decarg) > 1.): decarg=1. * decarg/np.abs(decarg)
    dec=math.asin(decarg)
    bdba=bdown/babs
    if (np.abs(bdba) > 1.): bdba=1. * bdba/np.abs(bdba)
    dip=math.asin(bdba)
    dipdiv=dip/np.sqrt(dip*dip+np.cos(xlati*UMR))
    if (np.abs(dipdiv) > 1.): dipdiv=1. * dipdiv/np.abs(dipdiv)
    smodip=math.asin(dipdiv)
    #       DIPL1=ATAN(0.5*TAN(DIP))/UMR
    dipl=math.atan(bdown/2.0/np.sqrt(bnorth*bnorth+beast*beast))/UMR
    ymodip=smodip/UMR
    dec=dec/UMR
    dip=dip/UMR
    return
#
#
# SHELLIG.FOR
#
# 11/01/91 SHELLG: lowest starting point for B0 search is 2  
#  1/27/92 Adopted to IGRF-91 coeffcients model
#  2/05/92 Reduce variable-names: INTER(P)SHC,EXTRA(P)SHC,INITI(ALI)ZE
#  8/08/95 Updated to IGRF-45-95; new coeff. DGRF90, IGRF95, IGRF95S
#  5/31/00 Updated to IGRF-45-00; new coeff.: IGRF00, IGRF00s
#  3/24/05 Updated to IGRF-45-10; new coeff.: IGRF05, IGRF05s
#  4/25/05 ENTRY FELDI(XI,H) and  DO 1111 I=1,7 [Alexey Petrov]
#  7/22/09 SHELLG: NMAX=13 for DGRF00 and IGRF05; H/G-arrays(195)
#  2/26/10 FELDCOF: Updated IGRF45-15; new coeff: DGRF05, IGRF10, IGRF10S
#  4/29/10 H/H-arrays(196); FELDCOF: corrected IGRF00 and ..00S
#  4/29/10 Change to new files dgrf%%%%.asc; new GETSHC; char*12 to 13
#
#*********************************************************************
#  SUBROUTINES SHELLG, STOER, FELDG, FELDCOF, GETSHC,                *
#       INTERSHC, EXTRASHC                                           *
#*********************************************************************
#*********************************************************************
#
#
#
def SHELLG(GLAT,GLON,ALT,FL,ICODE,B0):
    '''SUBROUTINE SHELLG(GLAT,GLON,ALT,DIMO,FL,ICODE,B0)
       -----------------------------------------------------------------------
           CALCULATES L-VALUE FOR SPECIFIED GEODAETIC COORDINATES, ALTITUDE
           AND GEMAGNETIC FIELD MODEL.
           REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTER, INTERNAL NOTE
           NO. 67, 1970.
           G. KLUGE, COMPUTER PHYSICS COMMUNICATIONS 3, 31-35, 1972
       -----------------------------------------------------------------------
          CHANGES (D. BILITZA, NOV 87):
              - USING CORRECT DIPOL MOMENT I.E.,DIFFERENT COMMON/MODEL/
              - USING IGRF EARTH MAGNETIC FIELD MODELS FROM 1945 TO 1990
              09/07/22 NMAX=13 for DGRF00 and IGRF05; H/G-arrays(195)
       -----------------------------------------------------------------------
           INPUT:  ENTRY POINT SHELLG
             GLAT  GEODETIC LATITUDE IN DEGREES (NORTH)
             GLON  GEODETIC LONGITUDE IN DEGREES (EAST)
             ALT   ALTITUDE IN KM ABOVE SEA LEVEL
           ENTRY POINT SHELLC
             V(3)  CARTESIAN COORDINATES IN EARTH RADII (6371.2 KM)
                     X-AXIS POINTING TO EQUATOR AT 0 LONGITUDE
                     Y-AXIS POINTING TO EQUATOR AT 90 LONG.
                     Z-AXIS POINTING TO NORTH POLE
             DIMO     DIPOL MOMENT IN GAUSS (NORMALIZED TO EARTH RADIUS)
             COMMON
             X(3)    NOT USED
             H(144)  FIELD MODEL COEFFICIENTS ADJUSTED FOR SHELLG
       -----------------------------------------------------------------------
          OUTPUT: FL           L-VALUE
                  ICODE        =1 NORMAL COMPLETION
                               =2 UNPHYSICAL CONJUGATE POINT (FL MEANINGLESS)
                               =3 SHELL PARAMETER GREATER THAN LIMIT UP TO
                                  WHICH ACCURATE CALCULATION IS REQUIRED;
                                  APPROXIMATION IS USED.
                  B0           MAGNETIC FIELD STRENGTH IN GAUSS
       -----------------------------------------------------------------------'''
       
    global X,H
    X = np.array(4)
#      COMMON/FIDB0/     SP	/CONST/UMR,PI      
    global UMR,PI
    global ERA,AQUAD,BQUAD,DIMO
#
#-- RMIN, RMAX ARE BOUNDARIES FOR IDENTIFICATION OF ICODE=2 AND 3
#-- STEP IS STEP SIZE FOR FIELD LINE TRACING
#-- STEQ IS STEP SIZE FOR INTEGRATION
# 
    RMIN = 0.05
    RMAX = 1.01
    STEP = 0.20
    STEQ = 0.03
    BEQU=1.0E10
    #*****ENTRY POINT  SHELLG  TO BE USED WITH GEODETIC CO-ORDINATES
    RLAT=GLAT*UMR
    CT=math.sin(RLAT)
    ST=math.cos(RLAT)
    D=math.sqrt(AQUAD-(AQUAD-BQUAD)*CT*CT)
    X[1]=(ALT+AQUAD/D)*ST/ERA
    X[3]=(ALT+BQUAD/D)*CT/ERA
    RLON=GLON*UMR
    X[2]=X[1]*math.sin(RLON)
    X[1]=X[1]*math.cos(RLON)
    #*****CONVERT TO DIPOL-ORIENTED CO-ORDINATES                     
    U = np.array([0.3511737,-0.9148385,-0.1993679,+0.9335804,+0.3583680,
                  0.0000000,+0.0714471,-0.1861260,+0.9799247])
    RQ=1./(X[1]*X[1]+X[2]*X[2]+X[3]*X[3])     #9
    R3H=math.sqrt(RQ*math.sqrt(RQ))                                      
    P[1,2]=(X[1]*U[1,1]+X[2]*U[2,1]+X[3]*U[3,1])*R3H
    P[2,2]=(X[1]*U[1,2]+X[2]*U[2,2])*R3H
    P[3,2]=(X[1]*U[1,3]+X[2]*U[2,3]+X[3]*U[3,3])*RQ
    #*****FIRST THREE POINTS OF FIELD LINE
    STEP=-STEP * P[3,2]/math.fabs(P[3,2])
    STOER(P[1,2],BQ2,R2)
    B0=math.sqrt(BQ2)
    P[1,3]=P[1,2]+0.5*STEP*P[4,2]
    P[2,3]=P[2,2]+0.5*STEP*P[5,2]
    P[3,3]=P[3,2]+0.5*STEP
    STOER(P[1,3],BQ3,R3)
    P[1,1]=P[1,2]-STEP*(2.*P[4,2]-P[4,3])
    P[2,1]=P[2,2]-STEP*(2.*P[5,2]-P[5,3])
    P[3,1]=P[3,2]-STEP
    STOER(P[1,1],BQ1,R1)
    P[1,3]=P[1,2]+STEP*(20.*P[4,3]-3.*P[4,2]+P[4,1])/18.
    P[2,3]=P[2,2]+STEP*(20.*P[5,3]-3.*P[5,2]+P[5,1])/18.
    P[3,3]=P[3,2]+STEP
    STOER(P[1,3],BQ3,R3)
    #*****INVERT SENSE IF REQUIRED
    if (BQ3 > BQ1):  #GOTO2
        STEP=-STEP
        R3=R1
        BQ3=BQ1
        for I in range(1,7):
            ZZ=P[I,1]
            P[I,1]=P[I,3]
            P[I,3]=ZZ
        #*****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH
    if (BQ1 < BEQU):
        BEQU=BQ1
        IEQU=1
    if (BQ2 < BEQU):
        BEQU=BQ2
        IEQU=2
    if (BQ3 < BEQU):
        BEQU=BQ3
        IEQU=3
    #*****INITIALIZATION OF INTEGRATION LOOPS
    STEP12=STEP/12.
    STEP2=STEP+STEP
    STEQ=STEQ * STEP/math.fabs(STEP)
    FI=0.
    ICODE=1
    ORADIK=0.
    OTERM=0.
    STP=R2*STEQ
    Z=P[3,2]+STP
    STP=STP/0.75
    P[8,1]=STEP2*(P[1,1]*P[4,1]+P[2,1]*P[5,1])
    P[8,2]=STEP2*(P[1,2]*P[4,2]+P[2,2]*P[5,2])
    #*****MAIN LOOP (FIELD LINE TRACING
    for N in range(3,3333):
        #*****CORRECTOR (FIELD LINE TRACING)
        P[1,N]=P[1,N-1]+STEP12*(5.*P[4,N]+8.*P[4,N-1]-P[4,N-2])
        P[2,N]=P[2,N-1]+STEP12*(5.*P[5,N]+8.*P[5,N-1]-P[5,N-2])
        #*****PREPARE EXPANSION COEFFICIENTS FOR INTERPOLATION
        #*****OF SLOWLY VARYING QUANTITIES
        P[8,N]=STEP2*(P[1,N]*P[4,N]+P[2,N]*P[5,N])
        C0=P[1,N-1]**2+P[2,N-1]**2
        C1=P[8,N-1]
        C2=(P[8,N]-P[8,N-2])*0.25
        C3=(P[8,N]+P[8,N-2]-C1-C1)/6.0
        D0=P[6,N-1]
        D1=(P[6,N]-P[6,N-2])*0.5
        D2=(P[6,N]+P[6,N-2]-D0-D0)*0.5
        E0=P[7,N-1]
        E1=(P[7,N]-P[7,N-2])*0.5
        E2=(P[7,N]+P[7,N-2]-E0-E0)*0.5
        #*****INNER LOOP (FOR QUADRATURE)
        goto10 = False
        goto30 = False
        T=(Z-P[3,N-1])/STEP  #4
        while True:
            if (T < 1.): break   #GOTO5
            HLI=0.5*(((C3*T+C2)*T+C1)*T+C0)
            ZQ=Z*Z
            R=HLI+math.sqrt(HLI*HLI+ZQ)
            if (R <= RMIN):      #GOTO30
                goto30 = True
                break
            RQ=R*R
            FF=math.sqrt(1.+3.*ZQ/RQ)
            RADIK=B0-((D2*T+D1)*T+D0)*R*RQ*FF
            if ((R-RMAX) > 0): 
                ICODE=2
                RADIK=RADIK-12.*(R-RMAX)**2
            if (RADIK+RADIK <= ORADIK):
                goto10 = True
                break
            TERM=math.sqrt(RADIK)*FF*((E2*T+E1)*T+E0)/(RQ+ZQ)
            FI=FI+STP*(OTERM+TERM)
            ORADIK=RADIK
            OTERM=TERM
            STP=R*STEQ
            Z=Z+STP
        #*****PREDICTOR (FIELD LINE TRACING)
        if not(goto10) and not(goto30):
            P[1,N+1]=P[1,N]+STEP12*(23.*P[4,N]-16.*P[4,N-1]+5.*P[4,N-2])  #5
            P[2,N+1]=P[2,N]+STEP12*(23.*P[5,N]-16.*P[5,N-1]+5.*P[5,N-2])
            P[3,N+1]=P[3,N]+STEP
            STOER(P[1,N+1],BQ3,R3)
            #*****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH
            if (BQ3 < BEQU):
                IEQU=N+1
                BEQU=BQ3
            #3     CONTINUE
    if not(goto30):
        if (IEQU < 2): IEQU=2      #10
        SP[1]=P[1,IEQU-1]
        SP[2]=P[2,IEQU-1]
        SP[3]=P[3,IEQU-1]
        if (ORADIK >= 1E-15):     #GOTO11
            FI=FI+STP/0.75*OTERM*ORADIK/(ORADIK-RADIK)
            #
            #-- The minimal allowable value of FI was changed from 1E-15 to 1E-12,
            #-- because 1E-38 is the minimal allowable arg. for ALOG in our envir.
            #-- D. Bilitza, Nov 87.
            #
        FI=0.5*math.fabs(FI)/math.sqrt(B0)+1E-12          #11
        #
        #*****COMPUTE L FROM B AND I.  SAME AS CARMEL IN INVAR.  
        #
        #-- Correct dipole moment is used here. D. Bilitza, Nov 87.
        #
        DIMOB0=DIMO/B0
        arg1=math.log(FI)
        arg2=math.log(DIMOB0)
        #      arg = FI*FI*FI/DIMOB0
        #      if(abs(arg).gt.88.0) arg=88.0
        XX=3*arg1-arg2
        if (XX > 23.0):       #GOTO 776
            GG=XX-3.0460681E0 
        elif (XX > 11.7):     #GOTO 775   
            GG=(((((2.8212095E-8*XX-3.8049276E-6)*XX+2.170224E-4)*XX-6.7310339E-3)*XX+1.2038224E-1)*XX-1.8461796E-1)*XX+2.0007187E0
        elif (XX > +3.0):     #GOTO 774
            GG=((((((((6.3271665E-10*XX-3.958306E-8)*XX+9.9766148E-07)*XX-1.2531932E-5)*XX+7.9451313E-5)*XX-3.2077032E-4)*XX+2.1680398E-3)*XX+1.2817956E-2)*XX+4.3510529E-1)*XX+6.222355E-1
        elif (XX > -3.0):     #GOTO 773
            GG=((((((((2.6047023E-10*XX+2.3028767E-9)*XX-2.1997983E-8)*XX-5.3977642E-7)*XX-3.3408822E-6)*XX+3.8379917E-5)*XX+1.1784234E-3)*XX+1.4492441E-2)*XX+4.3352788E-1)*XX+6.228644E-1
        elif (XX > -22.):     #GOTO 772
            GG=((((((((-8.1537735E-14*XX+8.3232531E-13)*XX+1.0066362E-9)*XX+8.1048663E-8)*XX+3.2916354E-6)*XX+8.2711096E-5)*XX+1.3714667E-3)*XX+1.5017245E-2)*XX+4.3432642E-1)*XX+6.2337691E-1
        else:                 # 771
            GG=3.33338E-1*XX+3.0062102E-1
        FL=EXP(ALOG((1.+EXP(GG))*DIMOB0)/3.0)   #777
        return
    #*****APPROXIMATION FOR HIGH VALUES OF L.
    ICODE=3        #30
    T=-P[3,N-1]/STEP
    FL=1./(math.fabs(((C3*T+C2)*T+C1)*T+C0)+1E-15)
    return                                                          
#
def SHELLC(V,FL,B0):
    '''*****ENTRY POINT  SHELLC  TO BE USED WITH CARTESIAN CO-ORDINATES'''
    X[1]=V[1]
    X[2]=V[2]
    X[3]=V[3]
    #*****CONVERT TO DIPOL-ORIENTED CO-ORDINATES                     
    U = np.array([0.3511737,-0.9148385,-0.1993679,+0.9335804,+0.3583680,
                  0.0000000,+0.0714471,-0.1861260,+0.9799247])
    RQ=1./(X[1]*X[1]+X[2]*X[2]+X[3]*X[3])     #9
    R3H=math.sqrt(RQ*math.sqrt(RQ))                                      
    P[1,2]=(X[1]*U[1,1]+X[2]*U[2,1]+X[3]*U[3,1])*R3H
    P[2,2]=(X[1]*U[1,2]+X[2]*U[2,2])*R3H
    P[3,2]=(X[1]*U[1,3]+X[2]*U[2,3]+X[3]*U[3,3])*RQ
    #*****FIRST THREE POINTS OF FIELD LINE
    STEP=-STEP * P(3,2)/math.fabs(P(3,2))
    STOER(P[1,2],BQ2,R2)
    B0=math.sqrt(BQ2)
    P[1,3]=P[1,2]+0.5*STEP*P[4,2]
    P[2,3]=P[2,2]+0.5*STEP*P[5,2]
    P[3,3]=P[3,2]+0.5*STEP
    STOER(P[1,3],BQ3,R3)
    P[1,1]=P[1,2]-STEP*(2.*P[4,2]-P[4,3])
    P[2,1]=P[2,2]-STEP*(2.*P[5,2]-P[5,3])
    P[3,1]=P[3,2]-STEP
    STOER(P[1,1],BQ1,R1)
    P[1,3]=P[1,2]+STEP*(20.*P[4,3]-3.*P[4,2]+P[4,1])/18.
    P[2,3]=P[2,2]+STEP*(20.*P[5,3]-3.*P[5,2]+P[5,1])/18.
    P[3,3]=P[3,2]+STEP
    STOER(P[1,3],BQ3,R3)
    #*****INVERT SENSE IF REQUIRED
    if (BQ3 > BQ1):  #GOTO2
        STEP=-STEP
        R3=R1
        BQ3=BQ1
        for I in range(1,7):
            ZZ=P[I,1]
            P[I,1]=P[I,3]
            P[I,3]=ZZ
        #*****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH
    if (BQ1 < BEQU):
        BEQU=BQ1
        IEQU=1
    if (BQ2 < BEQU):
        BEQU=BQ2
        IEQU=2
    if (BQ3 < BEQU):
        BEQU=BQ3
        IEQU=3
    #*****INITIALIZATION OF INTEGRATION LOOPS
    STEP12=STEP/12.
    STEP2=STEP+STEP
    STEQ=STEQ * STEP/math.fabs(STEP)
    FI=0.
    ICODE=1
    ORADIK=0.
    OTERM=0.
    STP=R2*STEQ
    Z=P[3,2]+STP
    STP=STP/0.75
    P[8,1]=STEP2*(P[1,1]*P[4,1]+P[2,1]*P[5,1])
    P[8,2]=STEP2*(P[1,2]*P[4,2]+P[2,2]*P[5,2])
    #*****MAIN LOOP (FIELD LINE TRACING
    for N in range(3,3333):
        #*****CORRECTOR (FIELD LINE TRACING)
        P[1,N]=P[1,N-1]+STEP12*(5.*P[4,N]+8.*P[4,N-1]-P[4,N-2])
        P[2,N]=P[2,N-1]+STEP12*(5.*P[5,N]+8.*P[5,N-1]-P[5,N-2])
        #*****PREPARE EXPANSION COEFFICIENTS FOR INTERPOLATION
        #*****OF SLOWLY VARYING QUANTITIES
        P[8,N]=STEP2*(P[1,N]*P[4,N]+P[2,N]*P[5,N])
        C0=P[1,N-1]**2+P[2,N-1]**2
        C1=P[8,N-1]
        C2=(P[8,N]-P[8,N-2])*0.25
        C3=(P[8,N]+P[8,N-2]-C1-C1)/6.0
        D0=P[6,N-1]
        D1=(P[6,N]-P[6,N-2])*0.5
        D2=(P[6,N]+P[6,N-2]-D0-D0)*0.5
        E0=P[7,N-1]
        E1=(P[7,N]-P[7,N-2])*0.5
        E2=(P[7,N]+P[7,N-2]-E0-E0)*0.5
        #*****INNER LOOP (FOR QUADRATURE)
        goto10 = False
        goto30 = False
        T=(Z-P[3,N-1])/STEP  #4
        while True:
            if (T < 1.): break   #GOTO5
            HLI=0.5*(((C3*T+C2)*T+C1)*T+C0)
            ZQ=Z*Z
            R=HLI+math.sqrt(HLI*HLI+ZQ)
            if (R <= RMIN):      #GOTO30
                goto30 = True
                break
            RQ=R*R
            FF=math.sqrt(1.+3.*ZQ/RQ)
            RADIK=B0-((D2*T+D1)*T+D0)*R*RQ*FF
            if ((R-RMAX) > 0): 
                ICODE=2
                RADIK=RADIK-12.*(R-RMAX)**2
            if (RADIK+RADIK <= ORADIK):
                goto10 = True
                break
            TERM=math.sqrt(RADIK)*FF*((E2*T+E1)*T+E0)/(RQ+ZQ)
            FI=FI+STP*(OTERM+TERM)
            ORADIK=RADIK
            OTERM=TERM
            STP=R*STEQ
            Z=Z+STP
        #*****PREDICTOR (FIELD LINE TRACING)
        if not(goto10) and not(goto30):
            P[1,N+1]=P[1,N]+STEP12*(23.*P[4,N]-16.*P[4,N-1]+5.*P[4,N-2])  #5
            P[2,N+1]=P[2,N]+STEP12*(23.*P[5,N]-16.*P[5,N-1]+5.*P[5,N-2])
            P[3,N+1]=P[3,N]+STEP
            STOER(P[1,N+1],BQ3,R3)
            #*****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH
            if (BQ3 < BEQU):
                IEQU=N+1
                BEQU=BQ3
            #3     CONTINUE
    if not(goto30):
        if (IEQU < 2): IEQU=2      #10
        SP[1]=P[1,IEQU-1]
        SP[2]=P[2,IEQU-1]
        SP[3]=P[3,IEQU-1]
        if (ORADIK >= 1E-15):     #GOTO11
            FI=FI+STP/0.75*OTERM*ORADIK/(ORADIK-RADIK)
            #
            #-- The minimal allowable value of FI was changed from 1E-15 to 1E-12,
            #-- because 1E-38 is the minimal allowable arg. for ALOG in our envir.
            #-- D. Bilitza, Nov 87.
            #
        FI=0.5*math.fabs(FI)/math.sqrt(B0)+1E-12          #11
        #
        #*****COMPUTE L FROM B AND I.  SAME AS CARMEL IN INVAR.  
        #
        #-- Correct dipole moment is used here. D. Bilitza, Nov 87.
        #
        DIMOB0=DIMO/B0
        arg1=math.log(FI)
        arg2=math.log(DIMOB0)
        #      arg = FI*FI*FI/DIMOB0
        #      if(abs(arg).gt.88.0) arg=88.0
        XX=3*arg1-arg2
        if (XX > 23.0):       #GOTO 776
            GG=XX-3.0460681E0 
        elif (XX > 11.7):     #GOTO 775   
            GG=(((((2.8212095E-8*XX-3.8049276E-6)*XX+2.170224E-4)*XX-6.7310339E-3)*XX+1.2038224E-1)*XX-1.8461796E-1)*XX+2.0007187E0
        elif (XX > +3.0):     #GOTO 774
            GG=((((((((6.3271665E-10*XX-3.958306E-8)*XX+9.9766148E-07)*XX-1.2531932E-5)*XX+7.9451313E-5)*XX-3.2077032E-4)*XX+2.1680398E-3)*XX+1.2817956E-2)*XX+4.3510529E-1)*XX+6.222355E-1
        elif (XX > -3.0):     #GOTO 773
            GG=((((((((2.6047023E-10*XX+2.3028767E-9)*XX-2.1997983E-8)*XX-5.3977642E-7)*XX-3.3408822E-6)*XX+3.8379917E-5)*XX+1.1784234E-3)*XX+1.4492441E-2)*XX+4.3352788E-1)*XX+6.228644E-1
        elif (XX > -22.):     #GOTO 772
            GG=((((((((-8.1537735E-14*XX+8.3232531E-13)*XX+1.0066362E-9)*XX+8.1048663E-8)*XX+3.2916354E-6)*XX+8.2711096E-5)*XX+1.3714667E-3)*XX+1.5017245E-2)*XX+4.3432642E-1)*XX+6.2337691E-1
        else:                 # 771
            GG=3.33338E-1*XX+3.0062102E-1
        FL=EXP(ALOG((1.+EXP(GG))*DIMOB0)/3.0)   #777
        return
    #*****APPROXIMATION FOR HIGH VALUES OF L.
    ICODE=3        #30
    T=-P[3,N-1]/STEP
    FL=1./(math.fabs(((C3*T+C2)*T+C1)*T+C0)+1E-15)
    return                                                          
#
def STOER(P,BQ,R):
    '''*******************************************************************
       * SUBROUTINE USED FOR FIELD LINE TRACING IN SHELLG                *
       * CALLS ENTRY POINT FELDI IN GEOMAGNETIC FIELD SUBROUTINE FELDG   *
          09/07/22 NMAX=13 for DGRF00 and IGRF05; H/G-arrays(195)
       *******************************************************************'''

    global XI,H
    #*****XM,YM,ZM  ARE GEOMAGNETIC CARTESIAN INVERSE CO-ORDINATES
    ZM=P[3]
    FLI=P[1]*P[1]+P[2]*P[2]+1E-15
    R=0.5*(FLI+math.sqrt(FLI*FLI+(ZM+ZM)**2))
    RQ=R*R
    WR=math.sqrt(R)
    XM=P[1]*WR
    YM=P[2]*WR
    #*****TRANSFORM TO GEOGRAPHIC CO-ORDINATE SYSTEM
    U = np.array([[+0.3511737,-0.9148385,-0.1993679],
                  [+0.9335804,+0.3583680,+0.0000000],
                  [+0.0714471,-0.1861260,+0.9799247]])
    XI[1]=XM*U[1,1]+YM*U[1,2]+ZM*U[1,3]
    XI[2]=XM*U[2,1]+YM*U[2,2]+ZM*U[2,3]
    XI[3]=XM*U[3,1]+ZM*U[3,3]
    #*****COMPUTE DERIVATIVES
    #      CALL FELDI(XI,H)
    FELDI
    Q=H[1]/RQ
    DX=H[3]+H[3]+Q*XI[1]
    DY=H[4]+H[4]+Q*XI[2]
    DZ=H[2]+H[2]+Q*XI[3]
    #*****TRANSFORM BACK TO GEOMAGNETIC CO-ORDINATE SYSTEM
    DXM=U[1,1]*DX+U[2,1]*DY+U[3,1]*DZ
    DYM=U[1,2]*DX+U[2,2]*DY
    DZM=U[1,3]*DX+U[2,3]*DY+U[3,3]*DZ
    DR=(XM*DXM+YM*DYM+ZM*DZM)/R
    #*****FORM SLOWLY VARYING EXPRESSIONS
    P[4]=(WR*DXM-0.5*P[1]*DR)/(R*DZM)
    P[5]=(WR*DYM-0.5*P[2]*DR)/(R*DZM)
    DSQ=RQ*(DXM*DXM+DYM*DYM+DZM*DZM)
    BQ=DSQ*RQ*RQ
    P[6]=math.sqrt(DSQ/(RQ+3.*ZM*ZM))
    P[7]=P[6]*(RQ+ZM*ZM)/(RQ*DZM)
    return                                                              
#
#
def FELDG(GLAT,GLON,ALT,BNORTH,BEAST,BDOWN,BABS):
    '''-----------------------------------------------------------------------
        CALCULATES EARTH MAGNETIC FIELD FROM SPHERICAL HARMONICS MODEL
        REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTRE, INTERNAL NOTE 61, 
                                                                       1970.
       -----------------------------------------------------------------------
       CHANGES (D. BILITZA, NOV 87):
           - FIELD COEFFICIENTS IN BINARY DATA FILES INSTEAD OF BLOCK DATA
           - CALCULATES DIPOL MOMENT
       09/07/22 NMAX=13 for DGRF00 and IGRF05; H/G-arrays(195)
       -----------------------------------------------------------------------
       INPUT:  GLAT  GEODETIC LATITUDE IN DEGREES (NORTH)
               GLON  GEODETIC LONGITUDE IN DEGREES (EAST)
               ALT   ALTITUDE IN KM ABOVE SEA LEVEL
               
       COMMON BLANK AND ENTRY POINT FELDI ARE NEEDED WHEN USED
            IN CONNECTION WITH L-CALCULATION PROGRAM SHELLG.
            
            COMMON /MODEL/ AND /IGRF1/
               UMR     = ATAN(1.0)*4./180.   <DEGREE>*UMR=<RADIANT>
               ERA     EARTH RADIUS FOR NORMALIZATION OF CARTESIAN 
                       COORDINATES (6371.2 KM)
               AQUAD, BQUAD   SQUARE OF MAJOR AND MINOR HALF AXIS OF 
                              EARTH ELLIPSOID AS RECOMMENDED BY INTERNAT. 
                              ASTRONOMICAL UNION (6378.160, 6356.775 KM).
                      NMAX    MAXIMUM ORDER OF SPHERICAL HARMONICS
                      TIME    YEAR (DECIMAL: 1973.5) FOR WHICH MAGNETIC 
                              FIELD IS TO BE CALCULATED
                      G(M)    NORMALIZED FIELD COEFFICIENTS (SEE FELDCOF)
                              M=NMAX*(NMAX+2)
       -----------------------------------------------------------------------
       OUTPUT: BABS   MAGNETIC FIELD STRENGTH IN GAUSS
           BNORTH, BEAST, BDOWN   COMPONENTS OF THE FIELD WITH RESPECT
                            TO THE LOCAL GEODETIC COORDINATE SYSTEM, WITH AXIS
                            POINTING IN THE TANGENTIAL PLANE TO THE NORTH, EAST
                            AND DOWNWARD.
       -----------------------------------------------------------------------'''
    V = np.zeros(3)
    B = np.zeros(3)
    global XI,H,NMAX,TIME,G,NAME,ERA,AQUAD,BQUAD,DIMO,UMR,PI#
#-- IS RECORDS ENTRY POINT
#
#*****ENTRY POINT  FELDG  TO BE USED WITH GEODETIC CO-ORDINATES         
    IS=1
    RLAT=GLAT*UMR
    CT=math.sin(RLAT)
    ST=math.cos(RLAT)
    D=math.sqrt(AQUAD-(AQUAD-BQUAD)*CT*CT)
    RLON=GLON*UMR
    CP=math.cos(RLON)
    SP=math.sin(RLON)
    ZZZ=(ALT+BQUAD/D)*CT/ERA
    RHO=(ALT+AQUAD/D)*ST/ERA
    XXX=RHO*CP
    YYY=RHO*SP
    RQ=1./(XXX*XXX+YYY*YYY+ZZZ*ZZZ)
    XI[1]=XXX*RQ
    XI[2]=YYY*RQ
    XI[3]=ZZZ*RQ
    
    IHMAX=NMAX*NMAX+1    #20    
    LAST=IHMAX+NMAX+NMAX
    IMAX=NMAX+NMAX-1
    for I in range(IHMAX,LAST): H[I]=G[I]     #8
    for K in [1,3,2]:   #6
        I=IMAX
        IH=IHMAX
        while True:
            IL=IH-I     #1
            F=2./np.float64(I-K+2)
            X=XI[1]*F
            Y=XI[2]*F
            Z=XI[3]*(F+F)
            I=I-2
            if (I-1) <=0:     #5
                H[IL]=G[IL]+Z*H[IH]+2.*(X*H[IH+1]+Y*H[IH+2])
                IH=IL
            elif (I == K): #4
                H[IL+2]=G[IL+2]+Z*H[IH+2]+X*H[IH+4]-Y*(H[IH+3]+H[IH])
                H[IL+1]=G[IL+1]+Z*H[IH+1]+Y*H[IH+4]+X*(H[IH+3]-H[IH])
                H[IL]=G[IL]+Z*H[IH]+2.*(X*H[IH+1]+Y*H[IH+2])
                IH=IL
            else : #2
                for M in [3,I,2]:
                    H[IL+M+1]=G[IL+M+1]+Z*H[IH+M+1]+X*(H[IH+M+3]-H[IH+M-1])-Y*(H[IH+M+2]+H[IH+M-2])
                    H[IL+M]=G[IL+M]+Z*H[IH+M]+X*(H[IH+M+2]-H[IH+M-2])+Y*(H[IH+M+3]+H[IH+M-1])
                H[IL+2]=G[IL+2]+Z*H[IH+2]+X*H[IH+4]-Y*(H[IH+3]+H[IH])
                H[IL+1]=G[IL+1]+Z*H[IH+1]+Y*H[IH+4]+X*(H[IH+3]-H[IH])
                H[IL]=G[IL]+Z*H[IH]+2.*(X*H[IH+1]+Y*H[IH+2])
                IH=IL
            if (I < K): #5 #GOTO 1
                break
        if (IS == 3): return
    S=.5*H[1]+2.*(H[2]*XI[3]+H[3]*XI[1]+H[4]*XI[2])
    T=(RQ+RQ)*math.sqrt(RQ)
    BXXX=T*(H[3]-S*XXX)
    BYYY=T*(H[4]-S*YYY)
    BZZZ=T*(H[2]-S*ZZZ)
    if (IS != 2): #GOTO 7
        BABS=math.sqrt(BXXX*BXXX+BYYY*BYYY+BZZZ*BZZZ)
        BEAST=BYYY*CP-BXXX*SP
        BRHO=BYYY*SP+BXXX*CP
        BNORTH=BZZZ*ST-BRHO*CT
        BDOWN=-BZZZ*CT-BRHO*ST
        return
    B[1]=BXXX    #7
    B[2]=BYYY
    B[3]=BZZZ
    return                                                          
#
def FELDI():                                                       
    '''-----------------------------------------------------------------------
        CALCULATES EARTH MAGNETIC FIELD FROM SPHERICAL HARMONICS MODEL
        REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTRE, INTERNAL NOTE 61, 
                                                                       1970.
       -----------------------------------------------------------------------
       CHANGES (D. BILITZA, NOV 87):
           - FIELD COEFFICIENTS IN BINARY DATA FILES INSTEAD OF BLOCK DATA
           - CALCULATES DIPOL MOMENT
       09/07/22 NMAX=13 for DGRF00 and IGRF05; H/G-arrays(195)
       -----------------------------------------------------------------------
       COMMON BLANK AND ENTRY POINT FELDI ARE NEEDED WHEN USED
            IN CONNECTION WITH L-CALCULATION PROGRAM SHELLG.
            
            COMMON /MODEL/ AND /IGRF1/
               UMR     = ATAN(1.0)*4./180.   <DEGREE>*UMR=<RADIANT>
               ERA     EARTH RADIUS FOR NORMALIZATION OF CARTESIAN 
                       COORDINATES (6371.2 KM)
               AQUAD, BQUAD   SQUARE OF MAJOR AND MINOR HALF AXIS OF 
                              EARTH ELLIPSOID AS RECOMMENDED BY INTERNAT. 
                              ASTRONOMICAL UNION (6378.160, 6356.775 KM).
                      NMAX    MAXIMUM ORDER OF SPHERICAL HARMONICS
                      TIME    YEAR (DECIMAL: 1973.5) FOR WHICH MAGNETIC 
                              FIELD IS TO BE CALCULATED
                      G(M)    NORMALIZED FIELD COEFFICIENTS (SEE FELDCOF)
                              M=NMAX*(NMAX+2)
       -----------------------------------------------------------------------
       OUTPUT: BABS   MAGNETIC FIELD STRENGTH IN GAUSS
           BNORTH, BEAST, BDOWN   COMPONENTS OF THE FIELD WITH RESPECT
                            TO THE LOCAL GEODETIC COORDINATE SYSTEM, WITH AXIS
                            POINTING IN THE TANGENTIAL PLANE TO THE NORTH, EAST
                            AND DOWNWARD.
       -----------------------------------------------------------------------'''
    V = np.zeros(3)
    B = np.zeros(3)
    global XI,H,NMAX,TIME,G,NAME,ERA,AQUAD,BQUAD,DIMO,UMR,PI
    
    IS=3                                                              
    IHMAX=NMAX*NMAX+1    #20    
    LAST=IHMAX+NMAX+NMAX
    IMAX=NMAX+NMAX-1
    for I in range(IHMAX,LAST): H[I]=G[I]     #8
    for K in [1,3,2]:   #6
        I=IMAX
        IH=IHMAX
        while True:
            IL=IH-I     #1
            F=2./np.float64(I-K+2)
            X=XI[1]*F
            Y=XI[2]*F
            Z=XI[3]*(F+F)
            I=I-2
            if (I-1) <=0:     #5
                H[IL]=G[IL]+Z*H[IH]+2.*(X*H[IH+1]+Y*H[IH+2])
                IH=IL
            elif (I == K): #4
                H[IL+2]=G[IL+2]+Z*H[IH+2]+X*H[IH+4]-Y*(H[IH+3]+H[IH])
                H[IL+1]=G[IL+1]+Z*H[IH+1]+Y*H[IH+4]+X*(H[IH+3]-H[IH])
                H[IL]=G[IL]+Z*H[IH]+2.*(X*H[IH+1]+Y*H[IH+2])
                IH=IL
            else : #2
                for M in [3,I,2]:
                    H[IL+M+1]=G[IL+M+1]+Z*H[IH+M+1]+X*(H[IH+M+3]-H[IH+M-1])-Y*(H[IH+M+2]+H[IH+M-2])
                    H[IL+M]=G[IL+M]+Z*H[IH+M]+X*(H[IH+M+2]-H[IH+M-2])+Y*(H[IH+M+3]+H[IH+M-1])
                H[IL+2]=G[IL+2]+Z*H[IH+2]+X*H[IH+4]-Y*(H[IH+3]+H[IH])
                H[IL+1]=G[IL+1]+Z*H[IH+1]+Y*H[IH+4]+X*(H[IH+3]-H[IH])
                H[IL]=G[IL]+Z*H[IH]+2.*(X*H[IH+1]+Y*H[IH+2])
                IH=IL
            if (I < K): #5 #GOTO 1
                break
        if (IS == 3): return
    S=.5*H[1]+2.*(H[2]*XI[3]+H[3]*XI[1]+H[4]*XI[2])
    T=(RQ+RQ)*math.sqrt(RQ)
    BXXX=T*(H[3]-S*XXX)
    BYYY=T*(H[4]-S*YYY)
    BZZZ=T*(H[2]-S*ZZZ)
    if (IS != 2): #GOTO 7
        BABS=math.sqrt(BXXX*BXXX+BYYY*BYYY+BZZZ*BZZZ)
        BEAST=BYYY*CP-BXXX*SP
        BRHO=BYYY*SP+BXXX*CP
        BNORTH=BZZZ*ST-BRHO*CT
        BDOWN=-BZZZ*CT-BRHO*ST
        return
    B[1]=BXXX    #7
    B[2]=BYYY
    B[3]=BZZZ
    return                                                          
#
def FELDC(V,B):                                                  
    '''-----------------------------------------------------------------------
        CALCULATES EARTH MAGNETIC FIELD FROM SPHERICAL HARMONICS MODEL
        REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTRE, INTERNAL NOTE 61, 
                                                                       1970.
       -----------------------------------------------------------------------
       CHANGES (D. BILITZA, NOV 87):
           - FIELD COEFFICIENTS IN BINARY DATA FILES INSTEAD OF BLOCK DATA
           - CALCULATES DIPOL MOMENT
       09/07/22 NMAX=13 for DGRF00 and IGRF05; H/G-arrays(195)
       -----------------------------------------------------------------------
       INPUT:
               V(3)  CARTESIAN COORDINATES IN EARTH RADII (6371.2 KM)
                       X-AXIS POINTING TO EQUATOR AT 0 LONGITUDE
                       Y-AXIS POINTING TO EQUATOR AT 90 LONG.
                       Z-AXIS POINTING TO NORTH POLE
                       
       COMMON BLANK AND ENTRY POINT FELDI ARE NEEDED WHEN USED
            IN CONNECTION WITH L-CALCULATION PROGRAM SHELLG.
            
            COMMON /MODEL/ AND /IGRF1/
               UMR     = ATAN(1.0)*4./180.   <DEGREE>*UMR=<RADIANT>
               ERA     EARTH RADIUS FOR NORMALIZATION OF CARTESIAN 
                       COORDINATES (6371.2 KM)
               AQUAD, BQUAD   SQUARE OF MAJOR AND MINOR HALF AXIS OF 
                              EARTH ELLIPSOID AS RECOMMENDED BY INTERNAT. 
                              ASTRONOMICAL UNION (6378.160, 6356.775 KM).
                      NMAX    MAXIMUM ORDER OF SPHERICAL HARMONICS
                      TIME    YEAR (DECIMAL: 1973.5) FOR WHICH MAGNETIC 
                              FIELD IS TO BE CALCULATED
                      G(M)    NORMALIZED FIELD COEFFICIENTS (SEE FELDCOF)
                              M=NMAX*(NMAX+2)
       -----------------------------------------------------------------------
       OUTPUT: B            MAGNETIC FIELD STRENGTH IN GAUSS
       -----------------------------------------------------------------------'''

    global XI,H,NMAX,TIME,G,NAME,ERA,AQUAD,BQUAD,DIMO,UMR,PI
    IS=2
    XXX=V[1]
    YYY=V[2]
    ZZZ=V[3]
    RQ=1./(XXX*XXX+YYY*YYY+ZZZ*ZZZ)
    XI[1]=XXX*RQ
    XI[2]=YYY*RQ
    XI[3]=ZZZ*RQ
    
    IHMAX=NMAX*NMAX+1    #20    
    LAST=IHMAX+NMAX+NMAX
    IMAX=NMAX+NMAX-1
    for I in range(IHMAX,LAST): H[I]=G[I]     #8
    for K in [1,3,2]:   #6
        I=IMAX
        IH=IHMAX
        while True:
            IL=IH-I     #1
            F=2./np.float64(I-K+2)
            X=XI[1]*F
            Y=XI[2]*F
            Z=XI[3]*(F+F)
            I=I-2
            if (I-1) <=0:     #5
                H[IL]=G[IL]+Z*H[IH]+2.*(X*H[IH+1]+Y*H[IH+2])
                IH=IL
            elif (I == K): #4
                H[IL+2]=G[IL+2]+Z*H[IH+2]+X*H[IH+4]-Y*(H[IH+3]+H[IH])
                H[IL+1]=G[IL+1]+Z*H[IH+1]+Y*H[IH+4]+X*(H[IH+3]-H[IH])
                H[IL]=G[IL]+Z*H[IH]+2.*(X*H[IH+1]+Y*H[IH+2])
                IH=IL
            else : #2
                for M in [3,I,2]:
                    H[IL+M+1]=G[IL+M+1]+Z*H[IH+M+1]+X*(H[IH+M+3]-H[IH+M-1])-Y*(H[IH+M+2]+H[IH+M-2])
                    H[IL+M]=G[IL+M]+Z*H[IH+M]+X*(H[IH+M+2]-H[IH+M-2])+Y*(H[IH+M+3]+H[IH+M-1])
                H[IL+2]=G[IL+2]+Z*H[IH+2]+X*H[IH+4]-Y*(H[IH+3]+H[IH])
                H[IL+1]=G[IL+1]+Z*H[IH+1]+Y*H[IH+4]+X*(H[IH+3]-H[IH])
                H[IL]=G[IL]+Z*H[IH]+2.*(X*H[IH+1]+Y*H[IH+2])
                IH=IL
            if (I < K): #5 #GOTO 1
                break
        if (IS == 3): return
    S=.5*H[1]+2.*(H[2]*XI[3]+H[3]*XI[1]+H[4]*XI[2])
    T=(RQ+RQ)*math.sqrt(RQ)
    BXXX=T*(H[3]-S*XXX)
    BYYY=T*(H[4]-S*YYY)
    BZZZ=T*(H[2]-S*ZZZ)
    if (IS != 2): #GOTO 7
        BABS=math.sqrt(BXXX*BXXX+BYYY*BYYY+BZZZ*BZZZ)
        BEAST=BYYY*CP-BXXX*SP
        BRHO=BYYY*SP+BXXX*CP
        BNORTH=BZZZ*ST-BRHO*CT
        BDOWN=-BZZZ*CT-BRHO*ST
        return
    B[1]=BXXX    #7
    B[2]=BYYY
    B[3]=BZZZ
    return                                                          
#
def FELDCOF(YEAR):
    '''-----------------------------------------------------------------------
            DETERMINES COEFFICIENTS AND DIPOL MOMENT FROM IGRF MODELS
        
        INPUT:  YEAR    DECIMAL YEAR FOR WHICH GEOMAGNETIC FIELD IS TO
                        BE CALCULATED
         				COMMON/IGRF1/ERAD,AQUAD,BQUAD,DIMO /CONST/UMR,PI
        OUTPUT:         COMMON/MODEL/NMAX,TIME,GH1,FIL1
         				COMMON/DIPOL/GHI1,GHI2,GHI3
                         
        THE GEOMAGNETIC DIPOL MOMENT (DIMO) IN GAUSS (NORMALIZED TO EARTH'S
        RADIUS) AT THE TIME (YEAR) IS COMPUTED BUT NOT USED.
        
        05/31/2000 updated to IGRF-2000 version (###)
        03/24/2000 updated to IGRF-2005 version (###)
        07/22/2009 NMAX=13 for DGRF00 and IGRF05; H/G-arrays(195)
        02/26/2010 update to IGRF-11 (2010) (###)
        10/05/2011 added COMMON/DIPOL/ for MLT computation in DPMTRX (IRIFUN)
        02/10/2015 update to IGRF-12 (2015) (###)
       -----------------------------------------------------------------------'''
    
    global NMAX,TIME,GH1,FIL1
    global ERAD,AQUAD,BQUAD,DIMO,UMR,PI
    global GHI1,GHI2,GHI3

    # ### FILMOD, DTEMOD array-size is number of IGRF maps
    GH1 = np.zeros(196)
    GH2 = np.zeros(196)
    # ### updated coefficient file names and corresponding years
    FILMOD = np.array(['dgrf1945.dat','dgrf1950.dat','dgrf1955.dat',
                       'dgrf1960.dat','dgrf1965.dat','dgrf1970.dat','dgrf1975.dat',
                       'dgrf1980.dat','dgrf1985.dat','dgrf1990.dat','dgrf1995.dat',
                       'dgrf2000.dat','dgrf2005.dat','dgrf2010.dat','igrf2015.dat',
                       'igrf2015s.dat'])
    DTEMOD = np.array([1945., 1950., 1955., 1960., 1965.,
                       1970., 1975., 1980., 1985., 1990., 1995., 2000.,2005.,
                       2010., 2015., 2020.])
    #
    # ### numye is number of IGRF coefficient files minus 1
    #
    NUMYE=15
    #
    #  IS=0 FOR SCHMIDT NORMALIZATION   IS=1 GAUSS NORMALIZATION
    #  IU  IS INPUT UNIT NUMBER FOR IGRF COEFFICIENT SETS
    #
    IU = 14
    IS = 0
    #-- DETERMINE IGRF-YEARS FOR INPUT-YEAR
    TIME = YEAR
    IYEA = int(YEAR/5.)*5
    L = (IYEA - 1945)/5 + 1
    if (L < 1): L=1
    if (L > NUMYE): L=NUMYE
    DTE1 = DTEMOD[L]
    FIL1 = FILMOD[L]
    DTE2 = DTEMOD[L+1]
    FIL2 = FILMOD[L+1]
    #-- GET IGRF COEFFICIENTS FOR THE BOUNDARY YEARS
    GETSHC (FIL1, NMAX1, ERAD, GH1, IER)
    if (IER != 0): return
    GETSHC (FIL2, NMAX2, ERAD, GH2, IER)
    if (IER != 0): return
    #-- DETERMINE IGRF COEFFICIENTS FOR YEAR
    if (L <= NUMYE-1):
        INTERSHC (YEAR, DTE1, NMAX1, GH1, DTE2,NMAX2, GH2, NMAX, GHA)
    else:
        EXTRASHC (YEAR, DTE1, NMAX1, GH1, NMAX2,GH2, NMAX, GHA)
    #-- DETERMINE MAGNETIC DIPOL MOMENT AND COEFFIECIENTS G
    F0=0.0
    for J in range (1,3):
        F = GHA[J] * 1.E-5
        F0 = F0 + F * F
    DIMO = np.sqrt(F0)
    GHI1=GHA[1]
    GHI2=GHA[2]
    GHI3=GHA[3]
    
    GH1[1] = 0.0
    I=2
    F0=1.E-5
    if (IS == 0): F0=-F0
    SQRT2=np.sqrt(2.)
    
    for N in range(1,NMAX):
        X = N
        F0 = F0 * X * X / (4.0 * X - 2.0)
        if (IS == 0): F0 = F0 * (2.0 * X - 1.0) / X
        F = F0 * 0.50
        if (IS == 0): F = F * SQRT2
        GH1[I] = GHA[I-1] * F0
        I = I+1
        for M in range(1,N):
            F = F * (X + M) / (X - M + 1.0)
            if (IS == 0): F = F * np.sqrt((X - M + 1.0) / (X + M))
            GH1[I] = GHA[I-1] * F
            GH1[I+1] = GHA[I] * F
            I=I+2
    return
#
#
def GETSHC (fileName, NMAX, ERAD, GH, IER):                                                                                           
# ===============================================================               
#       Reads spherical harmonic coefficients from the specified     
#       file into an array.                                          
#       Input:                                                       
#           IU    - Logical unit number                              
#           FSPEC - File specification                               
#       Output:                                                      
#           NMAX  - Maximum degree and order of model                
#           ERAD  - Earth's radius associated with the spherical     
#                   harmonic coefficients, in the same units as      
#                   elevation                                        
#           GH    - Schmidt quasi-normal internal spherical          
#                   harmonic coefficients                            
#           IER   - Error number: =  0, no error                     
#                                 = -2, records out of order         
#                                 = FORTRAN run-time error number    
# ===============================================================               
                                                                                
# not needed in Python    global konsol,mess
# ---------------------------------------------------------------               
#       Open coefficient file. Read past first header record.        
#       Read degree and order of model and Earth's radius.           
# ---------------------------------------------------------------               
    IER = 0
    f = open(fileName,'r')
    trash = f.readline()
    data = f.readline()
    temp = data.split(' ')
    NMAX = int(temp[1])
    ERAD = float(temp[3])
    xmyear = float(temp[4])
    nm = NMAX*(NMAX+2)
    GH = np.zeros(nm+1)
    temp = f.read()
    data = temp.split('\n')
    for i in range(1,nm+1): GH[i] = float(data[i-1])
#    IER =  ??? work needs to be done 
    f.close
    return
#
#
def INTERSHC (DATE, DTE1, NMAX1, GH1, DTE2, NMAX2, GH2, NMAX, GH):              
# ===============================================================               
#                                                                               
#       Version 1.01                                                 
#                                                                               
#       Interpolates linearly, in time, between two spherical        
#       harmonic models.                                             
#                                                                               
#       Input:                                                       
#           DATE  - Date of resulting model (in decimal year)        
#           DTE1  - Date of earlier model                            
#           NMAX1 - Maximum degree and order of earlier model        
#           GH1   - Schmidt quasi-normal internal spherical          
#                   harmonic coefficients of earlier model           
#           DTE2  - Date of later model                              
#           NMAX2 - Maximum degree and order of later model          
#           GH2   - Schmidt quasi-normal internal spherical          
#                   harmonic coefficients of later model             
#                                                                               
#       Output:                                                      
#           GH    - Coefficients of resulting model                  
#           NMAX  - Maximum degree and order of resulting model      
#                                                                               
#       A. Zunde                                                     
#       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225    
#                                                                               
# ===============================================================
               
# ---------------------------------------------------------------               
#       The coefficients (GH) of the resulting model, at date        
#       DATE, are computed by linearly interpolating between the     
#       coefficients of the earlier model (GH1), at date DTE1,       
#       and those of the later model (GH2), at date DTE2. If one     
#       model is smaller than the other, the interpolation is        
#       performed with the missing coefficients assumed to be 0.     
# ---------------------------------------------------------------
    GH = np.zeros(L+1)
    FACTOR = (DATE - DTE1) / (DTE2 - DTE1)
    if (NMAX1 == NMAX2):
        K = NMAX1 * (NMAX1 + 2)
        NMAX = NMAX1
    elif (NMAX1 > NMAX2):
        K = NMAX2 * (NMAX2 + 2)
        L = NMAX1 * (NMAX1 + 2)
        for I in range(K + 1, L): GH[I] = GH1[I] + FACTOR * (-GH1[I])
        NMAX = NMAX1
    else:
        K = NMAX1 * (NMAX1 + 2)
        L = NMAX2 * (NMAX2 + 2)
        for I in range(K + 1, L): GH[I] = FACTOR * GH2[I]
        NMAX = NMAX2
    for I in range(1, K): GH[I] = GH1[I] + FACTOR * (GH2[I] - GH1[I])
    return
#
#
def EXTRASHC (DATE, DTE1, NMAX1, GH1, NMAX2, GH2, NMAX, GH):                    
# ===============================================================               
#                                                                               
#       Version 1.01                                                   
#                                                                               
#       Extrapolates linearly a spherical harmonic model with a        
#       rate-of-change model.                                          
#                                                                               
#       Input:                                                         
#           DATE  - Date of resulting model (in decimal year)          
#           DTE1  - Date of base model                                 
#           NMAX1 - Maximum degree and order of base model             
#           GH1   - Schmidt quasi-normal internal spherical            
#                   harmonic coefficients of base model                
#           NMAX2 - Maximum degree and order of rate-of-change         
#                   model                                              
#           GH2   - Schmidt quasi-normal internal spherical            
#                   harmonic coefficients of rate-of-change model      
#                                                                               
#       Output:                                                        
#           GH    - Coefficients of resulting model                    
#           NMAX  - Maximum degree and order of resulting model        
#                                                                               
#       A. Zunde                                                       
#       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225      
#                                                                               
# ===============================================================               
# ---------------------------------------------------------------               
#       The coefficients (GH) of the resulting model, at date          
#       DATE, are computed by linearly extrapolating the coef-         
#       ficients of the base model (GH1), at date DTE1, using          
#       those of the rate-of-change model (GH2), at date DTE2. If      
#       one model is smaller than the other, the extrapolation is      
#       performed with the missing coefficients assumed to be 0.       
# ---------------------------------------------------------------               
                                                                                
    FACTOR = (DATE - DTE1)
    
    if (NMAX1 == NMAX2):
        K = NMAX1 * (NMAX1 + 2)
        NMAX = NMAX1
    elif (NMAX1 > NMAX2):
        K = NMAX2 * (NMAX2 + 2)
        L = NMAX1 * (NMAX1 + 2)
        for I in range(K + 1, L): GH[I] = GH1[I]
        NMAX = NMAX1
    else:
        K = NMAX1 * (NMAX1 + 2)
        L = NMAX2 * (NMAX2 + 2)
        for I in range(K + 1, L): GH[I] = FACTOR * GH2[I]
        NMAX = NMAX2
    for I in range(1, K): GH[I] = GH1[I] + FACTOR * GH2[I]
    
    return
#
#
def GEODIP(IYR,SLA,SLO,DLA,DLO,J):
    '''Calculates dipole geomagnetic coordinates from geocentric coordinates
       or vice versa.
                     J=0           J=1
		INPUT:     J,SLA,SLO     J,DLA,DLO
		OUTPUT:     DLA,DLO       SLA,SLO'''

#  Last revision: November 2005 (Vladimir Papitashvili)
#  The code is modifed from GEOCOR written by V.Popov and V.Papitashvili
#  in mid-1980s. 

    global UMR,PI 

#  Earth's radius (km) RE = 6371.2

#  The radius of the sphere to compute the coordinates (in Re)
#        RH = (RE + HI)/RE
    R = 1.
    
    if (j <= 0):
        COL = (90.- SLA)*UMR
        RLO = SLO*UMR
        SPHCAR(R,COL,RLO,X,Y,Z,1)
        GEOMAG(X,Y,Z,XM,YM,ZM,1,IYR)
        SPHCAR(RM,TH,PF,XM,YM,ZM,-1)
        SZM = ZM
        DLO = PF/UMR
        DCO = TH/UMR
        DLA = 90.- DCO
    else:
        COL = (90.- DLA)*UMR
        RLO = DLO*UMR
        SPHCAR(R,COL,RLO,XM,YM,ZM,1)
        GEOMAG(X,Y,Z,XM,YM,ZM,-1,IYR)
        SPHCAR(RM,TH,PF,X,Y,Z,-1)
        SZM = ZM
        SLO = PF/UMR
        SCO = TH/UMR
        SLA = 90.- SCO
    return
# 
# 
def fmodip(xlat):
    
    global xlong,year
    
    igrf_dip(xlat,xlong,year,300.,dec,dip,dipl,ymodip)
    fmodip=ymodip

    return
#
#
def GEOCGM01(ICOR,IYEAR,HI,DAT,PLA,PLO):
    '''*********************************************************************
       Version 2011 for GEO-CGM.FOR    (good through 2015)      January 2011
       Version 2005 for GEO-CGM.FOR    (good through 2010)     November 2005
       Nov 11, 2005  IGRF and RECALC are is modified to the IGRF-10 model
                     and extended back to 1900 using the DGRF coeffcients
       Apr 11, 2001  GEOLOW is modified to account for interpolation of
                     CGM meridians near equator across the 360/0 boundary
       AUTHORS:
           Natalia E. Papitashvili (WDC-B2, Moscow, Russia, now at NSSDC,
                         NASA/Goddard Space Flight Center, Greenbelt, Maryland)
           Vladimir O. Papitashvili (IZMIRAN, Moscow, Russia, now at SPRL,
                                             University of Michigan, Ann Arbor)
       Conributions from Boris A. Belov and Vladimir A. Popov (both at
                     IZMIRAN), Therese Moretto (DMI, DSRI, now at NSF)
                         Freddy Christiansen (DMI, DSRI), and 
                         Scott Boardsen (NASA/GSFC).
       
        The original version of this code is described in the brochure by
        N.A. Tsyganenko, A.V. Usmanov, V.O. Papitashvili, N.E. Papitashvili,
        and V.A. Popov, Software for computations of geomagnetic field and
        related coordinate systems, Soviet Geophys. Committ., Moscow, 58 pp.,
        1987. A number of subroutines from the revised GEOPACK-96 software
        package developed by Nikolai A. Tsyganenko and Mauricio Peredo are
        utilized in this code with some modifications (see full versions of
        GEOPACK packages on http://www-spof.gsfc.nasa.gov/Modeling/geopack.html).'''

#  This code consists of the main subroutine GEOCGM01, five functions
#  (OVL_ANG, CGMGLA, CGMGLO, DFRIDR, and AZM_ANG), eigth new and revised
#  subroutines from the above-mentioned brochure (MLTUT, MFC, FTPRNT,
#  GEOLOW, CORGEO, GEOCOR, SHAG, and RIGHT), and 9 subroutines from
#  GEOPACK-96 (IGRF, SPHCAR, BSPCAR, GEOMAG, MAGSM, SMGSM, RECALC, SUN)

#  =====================================================================

#  Input parameters:
#     ICOR = +1    geo to cgm
#            -1    cgm to geo
#     IYEAR= year
#     HI   = altitude in km
#  Input/Output parameters:
#     DAT(1,i)=slar geocentric latitude (input/output if icor=+1/-1)
#     DAT(2,i)=slor geocentric longitude (input/output if icor=+1/-1)
#     DAT(3,i)=clar CGM latitude (input/output if icor=-1/+1)
#     DAT(4,i)=clor CGM longitude (input/output if icor=-1/+1)
#  Output parameters:
#     DAT(5,i)=rbm apex of the magnetic field line in Re (Re=6371.2 km)
#            (this parameter approximately equals the McIlwain L-value)
#     DAT(6,i)=btr IGRF Magnetic field H (nT)
#     DAT(7,i)=brr IGRF Magnetic field D (deg)
#     DAT(8,i)=ovl oval_angle as the azimuth to "magnetic north":
#                + east in Northern Hemisphere
#                + west in Southern Hemisphere
#     DAT(9,i)=azm meridian_angle as the azimuth to the CGM pole:
#                + east in Northern Hemisphere
#                + west in Southern Hemisphere
#     DAT(10,i)=utm magnetic local time (MLT) midnight in UT hours
#     		 i=1	for the start point
#     		 i=2	for the conjugate point of the start point (slac, sloc)
#			 i=3    for the footprint at 1-Re of the start point (slaf,slof)
#			 i=4    for the conjugate footprint at 1-Re of the start point
#     PLA(1)	geocentric latitude of the CGM pole in the Northern hemisphere
#     PLO(1)	geocentric longitude of the CGM pole in the Northern hemisphere
#     PLA(2)	geocentric latitude of the CGM pole in the Southern hemisphere
#     PLO(2)	geocentric longitude of the CGM pole in the Southern hemisphere
#     PLA(3)	geoce lati CGM North pole at the Earth's surface 1-Re or zero alt.
#     PLO(3)	geoce long CGM North pole at the Earth's surface 1-Re or zero alt.
#     PLA(4)	geoce lati CGM South pole at the Earth's surface 1-Re or zero alt.
#     PLO(4)	geoce long CGM South pole at the Earth's surface 1-Re or zero alt.
#
# In program:
#     dla  = dipole latitude
#     dlo  = dipole longitude

#  =====================================================================

#      COMMON /C1/ AA(27),II(2),BB(8)
    global IYR, NM
#      COMMON /RZ/ RH
    PLA = np.zeros(4)
    PLO = np.zeros(4)
#  Year (for example, as for Epoch 1995.0 - no fraction of the year)

    IYR = iyear

#  Earth's radius (km)

    RE = 6371.2

#  NM is the number of harmonics

    NM = 10

#  The radius of the sphere to compute the coordinates (in Re)

    RH = (RE + HI)/RE

#  Correction of latitudes and longitudes if they are entered beyond of
#  the limits (this actually does not affect coordinate calculations
#  but the oval/meridian angles and MLT midnight cannot be computed)

    if (DAT[1,1] > 90.): DAT[1,1] =  180. - DAT[1,1]
    if (DAT[1,1] < -90.): DAT[1,1] = -180. - DAT[1,1]
    if (DAT[3,1] > 90.): DAT[3,1] =  180. - DAT[3,1]
    if (DAT[3,1] < -90.): DAT[3,1] = -180. - DAT[3,1]
    
    if (DAT[2,1] > 360.): DAT[2,1] = DAT[2,1] - 360.
    if (DAT[2,1] < -360.): DAT[2,1] = DAT[2,1] + 360.
    if (DAT[4,1] > 360.): DAT[4,1] = DAT[4,1] - 360.
    if (DAT[4,1] < -360.): DAT[4,1] = DAT[4,1] + 360.
#  Computation of CGM coordinates from geocentric ones at high- and
#  middle latitudes

    if (ICOR == 1):
        SLAR = DAT[1,1]
        SLOR = DAT[2,1]
        if (math.fabs(SLAR) == 90.): SLOR = 360.
        GEOCOR(SLAR,SLOR,RH,DLA,DLO,CLAR,CLOR,PMR)
        DAT[3,1] = CLAR
        DAT[4,1] = CLOR
    else:
#  Computation of geocentric coordinates from CGM ones at high- and
#  middle latitudes
        CLAR = DAT[3,1]
        CLOR = DAT[4,1]
        if (math.fabs(CLAR) == 90.): CLOR = 360.
        CORGEO(SLAR,SLOR,RH,DLA,DLO,CLAR,CLOR,PMR)
        DAT[1,1] = SLAR
        DAT[2,1] = SLOR

#  PMI is L-shell parameter for the magnetic field line; limit to 16 Re
    if (PMR >= 16.): PMR = 999.99
    DAT[5,1] = PMR
#  Check if CGM_Lat has been calculated, then go for the conjugate point
    if (CLAR > 999.):
#  CGM_Lat has NOT been calculated, call GEOLOW for computation of the
#  CGM coordinates at low latitudes using the CBM approach (see the
#  reference in GEOLOW)
        GEOLOW(SLAR,SLOR,RH,CLAR,CLOR,RBM,SLAC,SLOC)
        DAT[3,1] = CLAR
        DAT[4,1] = CLOR
        if (RBM >= 16.): RBM = 999.99
        DAT[5,1] = RBM
#  Conjugate point coordinates at low latitudes

#          WRITE(STR,'(2F6.2)') SLAC,SLOC
#          READ (STR,'(2F6.2)') SLAC,SLOC
        DAT[1,2] = SLAC
        DAT[2,2] = SLOC
        GEOCOR(SLAC,SLOC,RH,DAA,DOO,CLAC,CLOC,RBM)
        if (CLAC > 999.):
            GEOLOW(SLAC,SLOC,RH,CLAC,CLOC,RBM,SLAL,SLOL)
        DAT[3,2] = CLAC
        DAT[4,2] = CLOC
        DAT[5,2] = RBM
    else:
#  Computation of the magnetically conjugated point at high- and
#  middle latitudes
        CLAC = -CLAR
        CLOC =  CLOR
        DAT[3,2] = CLAC
        DAT[4,2] = CLOC
        CORGEO(SLAC,SLOC,RH,DAA,DOO,CLAC,CLOC,PMC)
        DAT[1,2] = SLAC
        DAT[2,2] = SLOC
        if (PMC >= 16.): PMC = 999.99
        DAT[5,2] = PMC
    
#  Same RBM for footprints as for the starting and conjugate points
    DAT[5,3] = DAT[5,1]
    DAT[5,4] = DAT[5,2]
#  Calculation of the magnetic field line footprint at the
#  Earth's surface for the starting point

    if (RH > 1.) and (CLAR < 999.) and (CLAR < 999.):
        FTPRNT(RH,SLAR,SLOR,CLAR,CLOR,ACLAR,ACLOR,SLARF,SLORF,1.)
        DAT[1,3] = SLARF
        DAT[2,3] = SLORF
        DAT[3,3] = ACLAR
        DAT[4,3] = ACLOR
#  and for the conjugate point
        FTPRNT(RH,SLAC,SLOC,CLAC,CLOC,ACLAC,ACLOC,SLACF,SLOCF,1.)
        DAT[1,4] = SLACF
        DAT[2,4] = SLOCF
        DAT[3,4] = ACLAC
        DAT[4,4] = ACLOC
    else:
        for i in range(1,4+1):
            for j in range(3,4+1): DAT[i,j] = 999.99
    
#  Computation of geocentric coordinates of the North or South CGM
#  poles for a given year at the altitude RH and Earth's surface (1-Re)

    CORGEO(PLAN,PLON,RH,DAA,DOO, 90.,360.,PMP)
    PLAN1 = PLAN
    PLON1 = PLON
    
    CORGEO(PLAS,PLOS,RH,DAA,DOO,-90.,360.,PMP)
    PLAS1 = PLAS
    PLOS1 = PLOS
    
    if (RH > 1.):
        CORGEO(PLAN1,PLON1,1.,DAA,DOO, 90.,360.,PMP)
        CORGEO(PLAS1,PLOS1,1.,DAA,DOO,-90.,360.,PMM)
    
    if (CLAR > 0.):
        PLA[1] = PLAS
        PLO[1] = PLOS
    else:
        PLA[1] = PLAN
        PLO[1] = PLON
    if (ACLAR > 0.):
        PLA[3] = PLAS1
        PLO[3] = PLOS1
    else:
        PLA[3] = PLAN1
        PLO[3] = PLON1
    if (CLAC < 0.):
        PLA[2] = PLAS
        PLO[2] = PLOS
    else:
        PLA[2] = PLAN
        PLO[2] = PLON
    if (ACLAC < 0.):
        PLA[4] = PLAS1
        PLO[4] = PLOS1
    else:
        PLA[4] = PLAN1
        PLO[4] = PLON1
        
    for j in range(1,4+1):
        DAT[6,j] = 99999.
        DAT[7,j] = 999.99
        DAT[8,j] = 99999.
        DAT[9,j] = 999.99
        DAT[10,j] = 999.99
        DAT[11,j] =  99.99
    icount = 2
    if (RH > 1.): icount = 4
    RJ = RH
    for j in range(1,icount+1):
        if (j > 2): RJ = 1.
        PLAJ = PLA[j]
        PLOJ = PLO[j]
        SLAJ = DAT[1,j]
        SLOJ = DAT[2,j]
        CLAJ = DAT[3,j]
        CLOJ = DAT[4,j]
#  Computation of the IGRF components
        MFC(SLAJ,SLOJ,RJ,BTR,BFR,BRR)
        DAT[6,j] = BTR
        DAT[7,j] = BFR
        DAT[8,j] = BRR
#  Computation of the oval_angle (OVL) between the tangents to
#  geographic and CGM latitudes at a given point (the code is slightly
#  modified from the source provided by Therese Morreto in 1994). Note
#  that rotation of OVL on 90 deg anticlockwise provides the azimuth
#  to the local "magnetic" north (south) measured from the local
#  geographic meridian. The OVL_ANG can be calculated only at middle
#  and high latitudes where CGM --> GEO is permitted.

        OVL = OVL_ANG(SLAJ,SLOJ,CLAJ,CLOJ,RJ)
        DAT[9,j] = OVL
#  Computation of the meridian_angle (AZM) between the geographic
#  meridian and direction (azimuth along the great-circle arc) to
#  the North (South) CGM pole

        AZM = AZM_ANG(SLAJ,SLOJ,CLAJ,PLAJ,PLOJ)
        DAT[10,j] = AZM
#  Computation of the MLT midnight (in UT)
        MLTUT(SLAJ,SLOJ,CLAJ,PLAJ,PLOJ,UT)
        DAT[11,j] = UT
#  End of loop j = 1,icount
    return
#
#
def OVL_ANG(sla,slo,cla,clo,rr):
    '''*********************************************************************
       This function returns an estimate at the given location of the angle
       (oval_angle) between the directions (tangents) along the constant
       CGM and geographic latitudes by utilizing the function DFRIDR from
       Numerical Recipes for FORTRAN.

       This angle can be taken as the azimuth to the local "magnetic" north
       (south) if the eastward (westward) tangent to the local CGM latitude
       points south (north) from the local geographic latitude.
       
       Written by Therese Moretto in August 1994 (revised by V. Papitashvili
       in January 1999).
       *********************************************************************'''
    #external cgmgla,cgmglo,dfridr
    global clat,cr360,cr0,rh

#  Ignore points which nearly coincide with the geographic or CGM poles
#  within 0.01 degree in latitudes; this also takes care if SLA or CLA
#  are dummy values (e.g., 999.99)

    if (math.abs(sla) >= 89.99) or (math.fabs(cla) >= 89.99) or (math.fabs(sla) < 30.):
        OVL_ANG = 999.99
        return
#  Initialize values for the cgmglo and cgmgla functions

    rh = rr
    clat = cla
    cr360 = False
    cr0 = False
#  Judge if SLO may be crossing the 360-0 limit. If geocentric
#  longitude of the location is larger than 270 deg, then cr360 is
#  set "true"; if it is less than 90 deg, then cr0 is set "true".

    if (slo >= 270.): cr360 = True
    if (slo <= 90.): cr0 = True
#  An initial stepsize (in degrees)

    step = 10.
#  Note that in the near-pole region the functions CGMGLA and CGMGLO
#  could be called from DFRIDR with the CGM latitudes exceeded 90 or
#  -90 degrees (e.g., 98 or -98) when STEP is added or subtracted to a
#  given CGM latitude (CLA). This does not produce discontinuities in
#  the functions because GEOCOR calculates GEOLAT smoothly for the
#  points lying behind the pole (e.g., as for 82 or - 82 deg. in the
#  above-mentioned example). However, it could be discontinuity in
#  GEOLON if |GEOLAT| = 90 deg. - see CGMGLO for details.

    hom = dfridr(cgmgla,clo,step,err1)
    denom = dfridr(cgmglo,clo,step,err2)
    denom = denom*cos(sla*0.017453293)
    OVL_ANG = -math.atan2(hom,denom)
    OVL_ANG = OVL_ANG*57.2957751
    return
#
#
def cgmgla(clon):
    '''*********************************************************************
       This function returns the geocentric latitude as a function of CGM
       longitude with the CGM latitude held in common block CGMGEO.
       Essentially this function just calls the subroutine CORGEO.
       *********************************************************************'''
    global cclat,cr360,cr0,rh
    
    rr = rh
    if (clon > 360.): clon = clon - 360.
    if (clon > 0.): clon = clon + 360.
    CORGEO(geolat,geolon,rr,dla,dlo,cclat,clon,pmi)
    cgmgla = geolat
    
    return
#
#
def cgmglo(clon):
    '''*********************************************************************
        Same as the function CGMGLA but this returns the geocentric
        longitude. If cr360 is true, geolon+360 deg is returned when geolon
        is less than 90 deg. If cr0 is true, geolon-360 deg is returned
        when geolon is larger than 270 degrees.
        *********************************************************************'''
    global cclat,cr360,cr0,rh

    rr = rh
    if(clon > 360.): clon = clon - 360.
    if(clon < 0.): clon = clon + 360.
#   1   continue
    CORGEO(geolat,geolon,rr,dla,dlo,cclat,clon,pmi)
    while (math.fabs(geolat) >= 89.99):
#  Geographic longitude geolon could be any number (e.g., discontinued)
#  when geolat is the geographic pole

        clon = clon - 0.01
	    
        if (cr360 and (geolon <= 90.)):
            cgmglo = geolon + 360.
        else:
            if (cr0 and (geolon >= 270.)):
                cgmglo = geolon - 360.
            else:
                cgmglo = geolon
        CORGEO(geolat,geolon,rr,dla,dlo,cclat,clon,pmi)
    return
#
#
def DFRIDR(func,x,h,err):
    '''**********************************************************************
       Numerical Recipes Fortran 77 Version 2.07
       Copyright (c) 1986-1995 by Numerical Recipes Software
       **********************************************************************'''
    CON=1.4
    CON2=CON*CON
    BIG=1.0E30
    NTAB=10
    SAFE=2.
#    EXTERNAL func
    global konsol,mess
    a = np.zeros([NTAB,NTAB])
    if (h == 0.):
#        if (mess) write(konsol,100) 
#100       FORMAT('h must be nonzero in dfridr')
        return
    hh = h
    a[1,1] = (func(x+hh)-func(x-hh))/(2.0*hh)
    err = BIG
    for i in range(2,NTAB+1):
        hh = hh/CON
        a[1,i] = (func(x+hh)-func(x-hh))/(2.0*hh)
        fac = CON2
        for j in range(2,i+1):
            a[j,i] = (a[j-1,i]*fac-a[j-1,i-1])/(fac-1.)
            fac = CON2*fac
            errt = np.max([math.fabs(a[j,i]-a[j-1,i]),math.fabs(a[j,i]-a[j-1,i-1])])
            if (errt <= err):
                err = errt
                dfridr = a[j,i]
        if (math.fabs(a[i,i]-a[i-1,i-1]) >= SAFE*err): return
    return
#
#
def AZM_ANG(sla,slo,cla,pla,plo):
    '''*********************************************************************
        Computation of an angle between the north geographic meridian and
        direction to the North (South) CGM pole: positive azimuth is
        measured East (West) from geographic meridian, i.e., the angle is
        measured between the great-circle arc directions to the geographic
        and CGM poles. In this case the geomagnetic field components in
        XYZ (NEV) system can be converted into the CGM system in both
        hemispheres as:
                           XM = X cos(alf) + Y sin(alf)
                           YM =-X sin(alf) + Y cos(alf)
        
        Written by V. O. Papitashvili in mid-1980s; revised in February 1999
        
        Ignore points which nearly coincide with the geographic or CGM poles
        within 0.01 degree in latitudes; this also takes care if SLA or CLA
        are dummy values (e.g., 999.99)
       *********************************************************************'''
       
    if (math.fabs(sla) >= 89.99) or (math.fabs(cla) >= 89.99):
        AZM_ANG = 999.99
        return
    sp = 1.
    ss = 1.
    if (np.abs(sp)*np.sign(pla)) != (np.abs(ss)*np.sign(cla)):
        print('WARNING - The CGM pole PLA and station CLAT are not in the same hemisphere: AZM_ANG is incorrect!')
#        write(7,2) pla,cla
#   2    format(/
#     +  'WARNING - The CGM pole PLA = ',f6.2,' and station CLAT = ',
#     +  f6.2,' are not in the same hemisphere: AZM_ANG is incorrect!')
      
    RAD = 0.017453293
    
    am = (90. - np.abs(pla))*rad
    if (np.abs(sp)*np.sign(pla) == np.abs(ss)*np.sign(sla)):
        cm = (90. - np.abs(sla))*rad
    else:
        cm = (90. + np.abs(sla))*rad
    if (sla >= 0.):
        bet = (plo - slo)*rad
    else:
        bet = (slo - plo)*rad
    sb = math.sin(bet)
    st = math.sin(cm)/math.tan(am) - math.cos(cm)*math.cos(bet)
    alfa = math.atan2(sb,st)
    AZM_ANG = alfa/rad

    return
#
#
def MLTUT(SLA,SLO,CLA,PLA,PLO,UT):
    '''*********************************************************************
        Calculates the MLT midnight in UT hours
        Definition of the MLT midnight (MLTMN) here is different from the
        approach described elsewhere. This definition does not take into
        account the geomagnetic meridian of the subsolar point which causes
        seasonal variations of the MLTMN in UT time. The latter approach is
        perfectly applicable to the dipole or eccentric dipole magnetic
        coordinates but it fails with the CGM coordinates because there are
        forbidden areas near the geomagnetic equator where CGM coordinates
        cannot be calculated by definition [e.g., Gustafsson et al., JATP,
        54, 1609, 1992].
        In this code the MLT midnight is defined as location of a given point
        on (or above) the Earth's surface strictly behind the North (South)
        CGM pole in such the Sun, the pole, and the point are lined up.
        This approach was originally proposed and coded by Boris Belov
        sometime in the beginning of 1980s; here it is slightly edited by
        Vladimir Papitashvili in February 1999.
        Ignore points which nearly coincide with the geographic or CGM poles
        within 0.01 degree in latitudes; this also takes care if SLA or CLA
        are dummy values (e.g., 999.99)
       *********************************************************************'''
       
    if (np.abs(sla) >= 89.99) or (np.abs(cla) >= 89.99):
        UT = 99.99
        return
    TPI = 6.283185307
    RAD = 0.017453293
    sp = 1.
    ss = 1.
    if (np.abs(sp)*np.sign(pla) != np.abs(ss)*np.sign(cla)):
        print('WARNING - The CGM pole PLA and station CLAT are not in the same hemisphere: MLTMN is incorrect!')
#        write(7,2) pla,cla
#   2    format(/
#     +  'WARNING - The CGM pole PLA = ',f6.2,' and station CLAT = ',
#     +  f6.2,' are not in the same hemisphere: MLTMN is incorrect!')

#  Solve the spherical triangle

    QQ = PLO*RAD
    CFF = 90. - np.abs(PLA)
    CFF = CFF*RAD
    if (CFF < 0.0000001): CFF=0.0000001
    if (np.abs(sp)*np.sign(pla) == np.abs(ss)*np.sign(sla)):
        CFT = 90. - np.abs(SLA)
    else:
        CFT = 90. + np.abs(SLA)
    CFT = CFT*RAD
    if (CFT < 0.0000001): CFT=0.0000001
    
    QT = SLO*RAD
    A = math.sin(CFF)/math.sin(CFT)
    Y = A*math.sin(QQ) - math.sin(QT)
    X = math.cos(QT) - A*math.cos(QQ)
    UT = math.atan2(Y,X)
    
    if (UT < 0.): UT = UT + TPI
    QQU = QQ + UT
    QTU = QT + UT
    BP = math.sin(CFF)*math.cos(QQU)
    BT = math.sin(CFT)*math.cos(QTU)
    UT = UT/RAD
    UT = UT/15.
    if (BP >= BT): #GOTO 10
        
        if (UT < 12.): UT = UT + 12.
        if (UT > 12.): UT = UT - 12.
    #10  CONTINUE
    
    return
#
#
def MFC(SLA,SLO,R,H,D,Z):
    '''*********************************************************************
        Computation of the IGRF magnetic field components
        Extracted as a subroutine from the earlier version of GEO-CGM.FOR
        V. Papitashvili, February 1999
        *********************************************************************'''
    global NM,IYR

#  This takes care if SLA or CLA are dummy values (e.g., 999.99)

    if (sla >= 999.):
        X = 99999.
        Y = 99999.
        Z = 99999.
        H = 99999.
        D = 999.99
        I = 999.99
        F = 99999.
        return
#  Computation of all geomagnetic field components
    RLA = (90.-SLA)*0.017453293
    RLO = SLO*0.017453293
    IGRF(IYR,NM,R,RLA,RLO,BR,BT,BF)
    X = -BT
    Y =  BF
    Z = -BR
    h = np.sqrt(X**2+Y**2)
    D = 57.2957751*math.atan2(Y,X)
    I = 57.2957751*math.atan2(Z,h)
    F = np.sqrt(h**2+Z**2)
    return
#
#
def FTPRNT(RH,SLA,SLO,CLA,CLO,ACLA,ACLO,SLAF,SLOF,RF):
    '''*********************************************************************
        Calculation of the magnetic field line footprint at the Earth's
        (or any higher) surface.
        Extracted as a subroutine from the earlier version of GEO-CGM.FOR by
        V. Papitashvili in February 1999 but then the subroutine was revised
        to obtain the Altitude Adjusted CGM coordinates. The AACGM approach
        is proposed by Kile Baker of the JHU/APL, see their World Wide Web
        site http://sd-www.jhuapl.edu/RADAR/AACGM/ for details.
        If RF = 1-Re (i.e., at the Earth's surface), then the footprint
        location is defined as the Altitude Adjusted (AA) CGM coordinates
        for a given point (ACLA, ACLO).
        
        If RF = 1.xx Re (i.e., at any altitude above or below the starting
        point), then the conjunction between these two points can be found
        along the field line.
       *********************************************************************'''
    global NM, IYR

#  This takes care if SLA or CLA are dummy values (e.g., 999.99)
    if (sla > 999.) or (cla > 999) or (RF == RH):
        ACLA = 999.99
        ACLO = 999.99
        SLAF = 999.99
        SLOF = 999.99
        return
#  Defining the Altitude Adjusted CGM coordinates for a given point
    COL = (90. - CLA)*0.017453293
    SN2 = (math.sin(COL))**2
    DECARG = np.sqrt((SN2*RF)/RH)
    if (np.abs(DECARG) > 1.): DECARG = np.sign(DECARG)*1.0
    ACOL = math.asin(DECARG)
    ACLA = 90. - ACOL*57.29577951
    if (CLA > 0.): ACLA = -ACLA
    ACLO = CLO
    
    CORGEO(SLAF,SLOF,RF,DLAF,DLOF,ACLA,ACLO,PMIF)
    
    if (SLAF < 999.): return
#  Tracing the magnetic field line down to the Earth's surface at low
#  latitudes if CORGEO failed to calculate geocentric coordinates SLAF
#  and SLOF

    if (SN2 < 0.0000001): SN2 = 0.0000001
    RL = RH/SN2
    FRAC = 0.03/(1.+3./(RL-0.6))
#  Checking direction of the magnetic field-line, so the step along
#  the field-line will go down, to the Earth surface

    if (CLA >= 0.): FRAC = -FRAC
    DS = RH*FRAC
#    250   CONTINUE

#  Start from an initial point
    R = RH
    RSLA = (90. - SLA)*0.0174533
    RSLO = SLO*0.0174533
    SPHCAR(R,RSLA,RSLO,XF,YF,ZF,1)
    RF1 = R
    XF1 = XF
    YF1 = YF
    ZF1 = ZF
    
#    255   CALL SHAG(XF,YF,ZF,DS)
    SHAG(XF,YF,ZF,DS)
    RR = np.sqrt(XF**2+YF**2+ZF**2)
    while True:
        if (RR > RH):
            R = RH
            RSLA = (90. - SLA)*0.0174533
            RSLO = SLO*0.0174533
            SPHCAR(R,RSLA,RSLO,XF,YF,ZF,1)
            RF1 = R
            XF1 = XF
            YF1 = YF
            ZF1 = ZF
            SHAG(XF,YF,ZF,DS)
            RR = np.sqrt(XF**2+YF**2+ZF**2)
            DS = -DS
            XF = XF1
            YF = YF1
            ZF = ZF1
        elif (RR > RF):
            SHAG(XF,YF,ZF,DS)
            RR = np.sqrt(XF**2+YF**2+ZF**2)
            RF1 = RR
            XF1 = XF
            YF1 = YF
            ZF1 = ZF
        else:
            break    
#    if (RR > RH):
#        DS = -DS
#        XF = XF1
#        YF = YF1
#        ZF = ZF1
#        GOTO 250
#    if (RR > RF):
#        RF1 = RR
#        XF1 = XF
#        YF1 = YF
#        ZF1 = ZF
#        GOTO 255
#    else:
    DR1 = np.abs(RF1 - RF)
    DR0 = np.abs( RF - RR)
    DR10 = DR1 + DR0
    if (DR10 != 0.):
        DS = DS*(DR1/DR10)
        SHAG(XF1,YF1,ZF1,DS)
    SPHCAR(RR,SLAF,SLOF,XF1,YF1,ZF1,-1)
    SLAF = 90. - SLAF*57.29578
    SLOF = SLOF*57.29578
    return
#
#
def GEOLOW(SLAR,SLOR,RH,CLAR,CLOR,RBM,SLAC,SLOC):
    '''*********************************************************************
       Calculates CGM coordinates from geocentric ones at low latitudes
       where the DGRF/IGRF magnetic field lines may never cross the dipole
       equatorial plane and, therefore, the definition of CGM coordinates
       becomes invalid.
       
       The code is written by Natalia and Vladimir Papitashvili as a part
       of the earlier versions of GEO-CGM.FOR; extracted as a subroutine by
       V. Papitashvili in February 1999.
       
       Apr 11, 2001  GEOLOW is modified to account for interpolation of
                     CGM meridians near equator across the 360/0 boundary

       See the paper by  Gustafsson, G., N. E. Papitashvili, and V. O.
       Papitashvili, A revised corrected geomagnetic coordinate system for
       Epochs 1985 and 1990 [J. Atmos. Terr. Phys., 54, 1609-1631, 1992]
       for detailed description of the B-min approach utilized here.
       *********************************************************************'''
       
    global NM,IYR
    
    BC = np.zeros(2)
    ARLAT = np.zeros(181)
    ARLON = np.zeros(181)
    
#  This takes care if SLA is a dummy value (e.g., 999.99)

    if (slar > 999.):
        CLAR = 999.99
        CLOR = 999.99
        SLAC = 999.99
        SLOC = 999.99
        RBM = 999.99
        return
    #endif

#  HH is an error (nT) to determine B-min along the magnetic field line

    DHH = 0.5

#  Filling the work arrays of CGM latitudes and longitudes with 999.99
#  Note that at certain geocentric longitudes in the very near-equator
#  region no "geomagnetic equator" can be defined at all.

    for J in range(61,121):
        ARLAT[J] = 999.99
        ARLON[J] = 999.99
    #ENDDO

    SLO = SLOR

    NDIR=0

#  Finding the geomagnetic equator as a projection of the B-min point
#  found for the field lines started from the last latitude in each
#  hemisphere where the CGM coordinates were obtained from geocentric
#  ones (GEO --> CGM). First the CGM coordinates are calculated in the
#  Northern (NDIR=0) and then in the Southern hemispheres (NDIR=1)

#  53     IF(NDIR.EQ.0) THEN
    while (NDIR == 0):

#  Program works from 30 deg. latitude down to the geographic equator
#  in the Northern Hemisphere

        for JC in range(61,91):
            SLA = 90.-(JC-1)
            GEOCOR(SLA,SLO,RH,DAA,DOO,CLA,CLO,PMM)
            if (CLA > 999.):
                NDIR=1
                #GOTO 53
                break
            #ENDIF
            ARLAT[JC] = CLA
            ARLON[JC] = CLO
        #ENDDO
        NDIR=1
        #GOTO 53
    else:
#  Program works from -30 deg. latitude down to the geographic equator
#  in the Southern Hemisphere
        for JC in range(121,92,-1):
            SLA = 90.-(JC-1)
            GEOCOR(SLA,SLO,RH,DAA,DOO,CLA,CLO,PMM)
            if (CLA > 999.):
                NDIR=0
                #GOTO 57
                break
            #ENDIF
            ARLAT[JC] = CLA
            ARLON[JC] = CLO
        #ENDDO
        if (CLA <= 999.): NDIR=0
    #ENDIF
#  57   CONTINUE

#  Finding last geographic latitudes along SLO where CGM coordinates
#  can be calculated

    n999=0
    ndir=0
    for jc in range(61,121):
        if (arlat[jc] > 999.):
            if (ndir == 0):
                jcn = jc - 1
                rnlat = arlat[jcn]
                rnlon = arlon[jcn]
                ndir = 1
                n999 = 1
             #endif
        #endif
        if (arlat[jc] < 999.):
            if (ndir == 1):
                jcs = jc
                rslat = arlat[jc]
                rslon = arlon[jc]
                ndir = 0
                #goto 59
                break
            #endif
        #endif
    #enddo
# 59     continue

#  If there is no points with 999.99 found along the SLO meridian,
#  then the IHEM loop will start from 3; otherwise it starts from 1

    if (n999 == 0):
        ih = 3
#        goto 31
    else:
        ih = 1
    #endif

#  Interpolation of the appropriate CGM longitudes between last
#  geocentric latitudes along SLO where CGM coordinates were defined
# (modified by Freddy Christiansen of DMI to account for interpolation
#  across the 360/0 boundary - April 11, 2001)

# this block of code should never be executed n999 ==0 jumps to 31.
#    if (n999 == 0):
#        rdel = jcs - jcn
#        if (rdel == 0.):
#            delon = 0.
#        else:
#            if (rslon > 270.) and (rnlon < 90.):
#                delon = (rslon - (rnlon + 360.))/rdel
#            else:
#                if (rslon < 90.) and (rnlon > 270.):
#                    delon = (rslon - (rnlon - 360.))/rdel
#                else:
#                    delon = (rslon - rnlon)/rdel
#                #endif
#            #endif
#        #endif
#        for jc in range(jcn+1,jcs-1):
#            arlon[jc] = rnlon + delon*(jc-jcn)
#            if (arlon[jc] < 0.): arlon[jc] = arlon[jc] + 360.
        #enddo
    
    #31   continue

#  Finding the CGM equator at SLO on the sphere with radius RH

    NOBM = 0
    for ihem in range(ih,3):
        RM = RH
#  Defining the real equator point from the Northern Hemisphere

        if (ihem == 1):
            CLA = rnlat
            SLA = 90. - (jcn - 1.)
            SLAN = SLA
        #endif
# Defining the real equator point from the Southern Hemisphere

        if (ihem == 2):
            CLA = rslat
            SLA = 90. - (jcs - 1)
            SLAS = SLA
        #endif
#  Defining the apex of the current magnetic field line

        if (ihem == 3):
            CLA = 0.
            SLA = SLAR
        #endif
#  Here CLA is used only to calculate FRAC

        COL = (90. - CLA)*0.017453293
        SLM = (90. - SLA)*0.017453293
        SLL = SLO*0.017453293
        IGRF(IYR,NM,RM,SLM,SLL,BR,BT,BF)
        SZ = -BR
        SPHCAR(RM,SLM,SLL,XGEO,YGEO,ZGEO,1)
        BM = np.sqrt(BR*BR + BT*BT + BF*BF)
        XBM = XGEO
        YBM = YGEO
        ZBM = ZGEO
        
        RL = 1./(np.sin(COL))**2
        FRAC = 0.03/(1. + 3./(RL - 0.6))
        if (SZ <= 0.): FRAC = -FRAC
        DSD = RL*FRAC
        DS = DSD

        #5      CONTINUE
        while True:
#  Keep two consequently computed points to define B-min

            #DO 7 I = 1,2
            for I in [1,2]:
                DD = DS
                SHAG(XGEO,YGEO,ZGEO,DD)
                #11         IF(I.NE.1) GOTO 9
                if (I != 1):
                    XBM1 = XGEO
                    YBM1 = YGEO
                    ZBM1 = ZGEO
                    RBM1 = np.sqrt(XBM1**2 + YBM1**2 + ZBM1**2)
                    #9         CONTINUE
                    SPHCAR(RM,SLM,SLL,XGEO,YGEO,ZGEO,-1)
                    IGRF(IYR,NM,RM,SLM,SLL,BR,BT,BF)

#  Go and compute the conjugate point if no B-min was found at this
#  magnetic field line (could happen at very near geomagnetic equator)

                    if (RM < RH):
                        NOBM = 1
                        break
                        #GOTO 77
                    #endif
                    BC[I] = np.sqrt(BR*BR + BT*BT + BF*BF)
                    #7      CONTINUE
            if (RM < RH): break

            B2 = BC[1]
            B3 = BC[2]
            if (BM > B2) and (B2 < B3): #GO TO 15
                BB3 = np.abs(B3 - B2)
                BB2 = np.abs(BM - B2)
                if (BB2 < DHH) and (BB3 < DHH): break #GO TO 21
            if (BM >= B2) and (B2 < B3): #GO TO 17
                BM = BM
                XGEO = XBM
                YGEO = YBM
                ZGEO = ZBM
                DS = DS/2.
                continue
            if (BM > B2) and (B2.LE.B3): #GO TO 17
                BM = BM
                XGEO = XBM
                YGEO = YBM
                ZGEO = ZBM
                DS = DS/2.
                continue
            BM = BC[1]
            XGEO = XBM1
            YGEO = YBM1
            ZGEO = ZBM1
            XBM = XBM1
            YBM = YBM1
            ZBM = ZBM1
        #GOTO 5
#   15      BB3 = ABS(B3 - B2)
#           BB2 = ABS(BM - B2)
#           IF(BB2.LT.DHH.AND.BB3.LT.DHH) GO TO 21
#   17      BM = BM
#           XGEO = XBM
#           YGEO = YBM
#           ZGEO = ZBM
#           DS = DS/2.
#           GOTO 5
        #21      CONTINUE
        
        if (RM >= RH):
            SPHCAR(RBM1,RLA,RLO,XBM1,YBM1,ZBM1,-1)
            RLA = 90. - RLA*57.2957751
            RLO = RLO*57.2957751
            
            if (ihem == 1): rlan = rla
            if (ihem == 2): rlas = rla

#  Computation of the magnetically conjugate point at low latitudes

            #54      continue
#   			print*,ihem
            if (ihem == 3):
                RBM = RBM1
                RM = RBM
                DS = DSD
            #55         continue
            SHAG(XBM1,YBM1,ZBM1,DS)
            RR = np.sqrt(XBM1**2 + YBM1**2 + ZBM1**2)
            while (RR > RH):
#              print*,rr,rh
#            if (RR > RH):
                R1 = RR
                X1 = XBM1
                Y1 = YBM1
                Z1 = ZBM1
                #GOTO 55
                SHAG(XBM1,YBM1,ZBM1,DS)
                RR = np.sqrt(XBM1**2 + YBM1**2 + ZBM1**2)
            else:
                DR1 = np.abs(RH - R1)
                DR0 = np.abs(RH - RR)
                DR10 = DR1 + DR0
                if (DR10 != 0.):
                    DS = DS*(DR1/DR10)
                    RM = R1
                    SHAG(X1,Y1,Z1,DS)
                #ENDIF
#                    print*,rr,slac,sloc,x1,y1,z1
                SPHCAR(RR,SLAC,SLOC,X1,Y1,Z1,-1)
#                    print*,rr,slac,sloc,x1,y1,z1
                SLAC = 90. - SLAC*57.2957751
                SLOC = SLOC*57.2957751
            #ENDIF
        #endif
    #End of loop IHEM
    #77      continue
    #enddo

    if (n999 != 0): #goto 91
        if (NOBM == 1):
#  Interpolation of CGM latitudes if there is no B-min at this
#  magnetic field line
            rdel = jcs - jcn
            if (rdel == 0.):
                delat = 0.
            else:
                delat = (rslat - rnlat)/rdel
            #endif
            jdel = 0
            for jc in range(jcn+1,jcs-1):
                jdel = jdel + 1
                arlat[jc] = rnlat + delat*jdel
            #enddo
            RBM = 999.99
            SLAC = 999.99
            SLOC = 999.99
        else:          
# Geocentric latitude of the CGM equator
            rla = (rlan + rlas)/2.
#  Interpolation of the CGM latitudes in the Northern hemisphere
            rdel = SLAN - rla
            if (rdel == 0.):
                delat = 0.
            else:
                delat = rnlat/rdel
            #endif
            jdn = np.abs(rdel)
            jdel = 0
            for jc in range(jcn+1,jcn+jdn):
                jdel = jdel + 1
                arlat[jc] = rnlat - delat*jdel
            #enddo
#  Interpolation of the CGM latitudes in the Southern hemisphere

            rdel = SLAS - rla
            if (rdel == 0.):
                delat = 0.
            else:
                delat = rslat/rdel
            #endif
            jds = np.abs(rdel)
            jdel = 0
            for jc in range(jcs-1,jcs-jds,-1):
                jdel = jdel + 1
                arlat[jc] = rslat + delat*jdel
            #enddo
        # ENDIF
    #91 continue
#  Defining by interpolation the exact values of the CGM latitude
#  and longitude between two adjacent values

    L1 = 90. - SLAR + 1.
    if (SLAR < 0.):
        L2 = L1-1
    else:
        L2 = L1+1
    #ENDIF
    DSLA = np.abs(SLAR - np.int(SLAR))
    DELCLA = ARLAT[L2] - ARLAT[L1]
    DELCLO = ARLON[L2] - ARLON[L1]
    CLAR = ARLAT[L1] + DELCLA*DSLA
    CLOR = ARLON[L1] + DELCLO*DSLA
    return
#
#
def CORGEO(SLA,SLO,RH,DLA,DLO,CLA,CLO,PMI):
    '''*********************************************************************
        Calculates geocentric coordinates from corrected geomagnetic ones.
        The code is written by Vladimir Popov and Vladimir Papitashvili
        in mid-1980s; revised by V. Papitashvili in February 1999
        *********************************************************************'''
    global NM, IYR
#  This takes care if CLA is a dummy value (e.g., 999.99)

    jc = 0
    if (np.abs(cla) < 0.1):
#        write(7,2)
#   2    format(/
#     +'WARNING - No calculations within +/-0.1 degree near CGM equator')
        jc = 1
    #endif
    if (cla > 999.) or (jc == 1):
        SLA = 999.99
        SLO = 999.99
        DLA = 999.99
        DLO = 999.99
        PMI = 999.99
        return
    #endif
    
    NG = NM
    COL = 90. - CLA
    R = 10.
    R1 = R
    R0 = R
    COL = COL*0.017453293
    RLO = CLO*0.017453293
    SN = np.sin(COL)
    SN2 = SN*SN

#  The CGM latitude should be at least 0.01 deg. away of the CGM pole

    if (SN2 < 0.000000003): SN2 = 0.000000003
#      RFI = 1./SN2
    RFI = RH/SN2
    PMI = RFI
    if (PMI > 99.999): PMI = 999.99
    AA10 = R/RFI

#  RFI = R if COL = 90 deg.

    if (RFI > R): #GOTO 1
        SAA = AA10/(1.-AA10)
        SAQ = SQRT(SAA)
        SCLA = ATAN(SAQ)
        if (CLA < 0): SCLA = 3.14159265359 - SCLA
    else: #GOTO 3
        #1   SCLA = 1.57079632679
        SCLA = 1.57079632679
        R0 = RFI
    # 3 CALL SPHCAR(R0,SCLA,RLO,XM,YM,ZM,1)
    SPHCAR(R0,SCLA,RLO,XM,YM,ZM,1)
    GEOMAG(X,Y,Z,XM,YM,ZM,-1,IYR)
    RL = R0
    FRAC = -0.03/(1. + 3./(RL - 0.6))
    if (CLA < 0.): FRAC = -FRAC
    R = R0
    # 5    DS = R*FRAC
    while True:
        DS = R*FRAC
        NM = (1. + 9./R) + 0.5
        SHAG(X,Y,Z,DS)
        R = np.sqrt(X**2+Y**2+Z**2)
        if (R <= RH): break #GOTO 7
        R1 = R
        X1 = X
        Y1 = Y
        Z1 = Z
        #GOTO 5

#  Define intersection with the start surface

#    7   DR1 = ABS(RH - R1)
    DR1 = np.abs(RH - R1)
    DR0 = np.abs(RH - R)
    DR10 = DR1 + DR0
    if (DR10 != 0.):
        DS = DS*(DR1/DR10)
        SHAG(X1,Y1,Z1,DS)
    #ENDIF
    SPHCAR(R,GTET,GXLA,X1,Y1,Z1,-1)
    GTH = GTET*57.2957751
    SLO = GXLA*57.2957751
    SLA = 90. - GTH
    GEOMAG(X1,Y1,Z1,XM,YM,ZM,1,IYR)
    SPHCAR(RM,TH,PF,XM,YM,ZM,-1)
    DLO = PF*57.2957751
    DLA = 90. - TH*57.2957751

    NM = NG

#  Because CORGEO cannot check if the CGM --> GEO transformation is
#  performed correctly in the equatorial area (that is, where the IGRF
#  field line may never cross the dipole equatorial plane). Therefore,
#  the backward check is required for geocentric latitudes lower than
#  30 degrees (see the paper referenced in GEOLOW)

    if (np.abs(SLA) < 30.) or (np.abs(CLA) < 30.):
        GEOCOR(SLA,SLO,RH,DLS,DLS,CLAS,CLOS,PMS)
        
        if (CLAS > 999.): GEOLOW(SLA,SLO,RH,CLAS,CLOS,RBM,SLAC,SLOC)
        if (np.abs(np.abs(CLA)-np.abs(CLAS)) >= 1.):
#            write(7,22) CLA
#            22    format(/
#     +'WARNING - Selected CGM_Lat.=',f6.2,' is too close to geomagnetic'
#     +/'          equator where CGM coordinates are not defined')
            SLA = 999.99
            SLO = 999.99
            PMI = 999.99
        #ENDIF
    #ENDIF
    return
#
#
def GEOCOR(SLA,SLO,RH,DLA,DLO,CLA,CLO,PMI):
    '''*********************************************************************
        Calculates corrected geomagnetic coordinates from geocentric ones
        The code is written by Vladimir Popov and Vladimir Papitashvili
        in mid-1980s; revised by V. Papitashvili in February 1999
       *********************************************************************'''
       
    global NM,IYR

#  This takes care if SLA is a dummy value (e.g., 999.99)

    if (sla > 999.):
        CLA = 999.99
        CLO = 999.99
        DLA = 999.99
        DLO = 999.99
        PMI = 999.99
        return
    #endif

    NG = NM
    COL = 90. - SLA
    R = RH
    R1 = R
    COL = COL*0.017453293
    RLO = SLO*0.017453293
    SPHCAR(R,COL,RLO,X,Y,Z,1)
    GEOMAG(X,Y,Z,XM,YM,ZM,1,IYR)
    SPHCAR(RM,TH,PF,XM,YM,ZM,-1)
    SZM = ZM
    DLO = PF*57.2957751
    DCO = TH*57.2957751
    DLA = 90. - DCO
    RL = R/(np.sin(TH))**2
    FRAC = 0.03/(1. + 3./(RL - 0.6))
    
    if (SZM < 0.): FRAC = -FRAC

#  Error to determine the dipole equtorial plane: aprox. 0.5 arc min

    HHH = 0.0001571

#  Trace the IGRF magnetic field line to the dipole equatorial plane

#   1     DS = R*FRAC
#   3     NM = (1. + 9./R) + 0.5
    while True:
        NM = (1. + 9./R) + 0.5
        R1 = R
        X1 = X
        Y1 = Y
        Z1 = Z
        SHAG(X,Y,Z,DS)
        GEOMAG(X,Y,Z,XM,YM,ZM,1,IYR)
        SPHCAR(R,C,S,XM,YM,ZM,-1)

#  As tracing goes above (RH+10_Re), use the dipole field line

        if (R > 10.+RH): break #GOTO 9

#  If the field line returns to the start surface without crossing the
#  dipole equatorial plane, no CGM coordinates can be calculated

        if (R <= RH): break #GOTO 11

        DCL = C - 1.5707963268
        if (np.abs(DCL) <= HHH): break #GOTO 9
        RZM = ZM
        if (SZM > 0.) and (RZM > 0.):
            DS = R*FRAC
            continue #GOTO 1
        if (SZM < 0.) and (RZM < 0.):
            DS = R*FRAC
            continue #GOTO 1
        R = R1
        X = X1
        Y = Y1
        Z = Z1
        DS = DS/2.
        #GOTO 3
    
    if (R > 10.+RH):
#        9  CALL GEOMAG(X,Y,Z,XM,YM,ZM,1,IYR)
        GEOMAG(X,Y,Z,XM,YM,ZM,1,IYR)
        SPHCAR(R,GTET,GXLA,XM,YM,ZM,-1)
        ST = np.abs(np.sin(GTET))
        RRH = np.abs(RH/(R - RH*ST**2))
        CLA = 1.5707963 - np.atan(ST*SQRT(RRH))
        CLA = CLA*57.2957751
        CLO = GXLA*57.2957751
        if (SZM < 0.): CLA = -CLA
        SSLA = 90. - CLA
        SSLA = SSLA*0.017453293
        SN = np.sin(SSLA)
#       PMI = 1/(SN*SN)
        PMI = RH/(SN*SN)
#        GOTO 13
    elif (R <= RH):
#   11   CLA = 999.99
        CLA = 999.99
        CLO = 999.99
        PMI = 999.99
#   13    NM = NG
    NM = NG
    
    return
#
#
def SHAG(X,Y,Z,DS):
    '''*********************************************************************
        Similar to SUBR STEP from GEOPACK-1996 but SHAG takes into account
        only internal sources
        The code is re-written from Tsyganenko's subroutine STEP by
        Natalia and Vladimir Papitashvili in mid-1980s
       *********************************************************************'''
    global DS3
    
    DS3 = -DS/3.
    RIGHT(X,Y,Z,R11,R12,R13)
    RIGHT(X+R11,Y+R12,Z+R13,R21,R22,R23)
    RIGHT(X+.5*(R11+R21),Y+.5*(R12+R22),Z+.5*(R13+R23),R31,R32,R33)
    RIGHT(X+.375*(R11+3.*R31),Y+.375*(R12+3.*R32),Z+.375*(R13+3.*R33),R41,R42,R43)
    RIGHT(X+1.5*(R11-3.*R31+4.*R41),Y+1.5*(R12-3.*R32+4.*R42),Z+1.5*(R13-3.*R33+4.*R43),
          R51,R52,R53)
    X = X+.5*(R11+4.*R41+R51)
    Y = Y+.5*(R12+4.*R42+R52)
    Z = Z+.5*(R13+4.*R43+R53)
    
    return
#
#
def RIGHT(X,Y,Z,R1,R2,R3):
    '''*********************************************************************
        Similar to SUBR RHAND from GEOPACK-1996 but RIGHT takes into account
        only internal sources
        The code is re-written from Tsyganenko's subroutine RHAND
        by Natalia and Vladimir Papitashvili in mid-1980s
       *********************************************************************'''
    global DS3,NM,IYR
    
    SPHCAR(R,T,F,X,Y,Z,-1)
    IGRF(IYR,NM,R,T,F,BR,BT,BF)
    BSPCAR(T,F,BR,BT,BF,BX,BY,BZ)
    B = DS3/SQRT(BX**2+BY**2+BZ**2)
    R1 = BX*B
    R2 = BY*B
    R3 = BZ*B

    return
#
#
def IGRF(IY,NM,R,T,F,BR,BT,BF):
    '''*********************************************************************
        CALCULATES COMPONENTS OF THE MAIN (INTERNAL) GEOMAGNETIC FIELD IN SPHERICAL
        GEOGRAPHICAL COORDINATE SYSTEM, USING IAGA INTERNATIONAL GEOMAGNETIC REFERENCE
        MODEL COEFFICIENTS (e.g., http://www.ngdc.noaa.gov/IAGA/wg8/igrf2000.html)
        
        UPDATING THE COEFFICIENTS TO A GIVEN EPOCH IS MADE AUTOMATICALLY UPON THE FIRST
        CALL AND AFTER EVERY CHANGE OF THE PARAMETER IY.
        
        -----INPUT PARAMETERS:
            
            IY  -  YEAR NUMBER (FOUR-DIGIT; 1965 &LE IY &LE 2005)
            NM  -  HIGHEST ORDER OF SPHERICAL HARMONICS IN THE SCALAR POTENTIAL (NM &LE 10)
            R,T,F -  SPHERICAL COORDINATES (RADIUS R IN UNITS RE=6371.2 KM, GEOGRAPHIC
                COLATITUDE  T  AND LONGITUDE  F  IN RADIANS)
            
        -----OUTPUT PARAMETERS:
            
            BR,BT,BF - SPHERICAL COMPONENTS OF THE MAIN GEOMAGNETIC FIELD IN NANOTESLA
            
        LAST MODIFICATION:  JANUARY 5, 2001, BY: N. A. TSYGANENKO
        THE CODE WAS MODIFIED TO ACCEPT DATES THROUGH 2005.
        IT HAS ALSO BEEN SLIGHTLY SIMPLIFIED BY TAKING OUT SOME REDUNDANT STATEMENTS,
        AND A "SAVE" STATEMENT WAS ADDED, TO AVOID POTENTIAL PROBLEMS WITH SOME
        FORTRAN COMPILERS.
        
        MODIFIED TO DGRF TO ACCEPT YEARS FROM 1900 THROUGH 2005
        BY SCOTT BOARDSEN, NASA GSFC, OCTOBER 2004

        MODIFIED TO IGRF-10 WITH YEARS THROUGH 2010
        BY V. PAPITASHVILI, NOVEMBER 2005

        MODIFIED TO IGRF-11 WITH YEARS THROUGH 2015
        BY V. PAPITASHVILI, January 2011
        
        MODIFIED TO IGRF-12 WITH YEARS THROUGH 2020
        BY D. Bilitza, July 2017
       *********************************************************************'''

#      SAVE MA,IYR,G,H,REC

#      DIMENSION A(11),B(11),DG(45),DH(45),G(66),H(66),REC(66),
#     * G1900(66),G1905(66),G1910(66),G1915(66),G1920(66),G1925(66),
#     * G1930(66),G1935(66),G1940(66),G1945(66),G1950(66),G1955(66),
#     * G1960(66),G1965(66),G1970(66),G1975(66),G1980(66),G1985(66),
#     * G1990(66),G1995(66),G2000(66),G2005(66),G2010(66),G2015(66),
#     * H1900(66),H1905(66),H1910(66),H1915(66),H1920(66),H1925(66),
#     * H1930(66),H1935(66),H1940(66),H1945(66),H1950(66),H1955(66),
#     * H1960(66),H1965(66),H1970(66),H1975(66),H1980(66),H1985(66),
#     * H1990(66),H1995(66),H2000(66),H2005(66),H2010(66),H2015(66)
#
#      logical	mess
#    SAVE MA,IYR,G,H,REC       
    global konsol,mess
#    REC = np.zeros(66)
      
    G1900 = np.array([0., -31543.,  -2298.,   -677.,   2905.,    924.,   1022.,
                      -1469.,   1256.,    572.,    876.,    628.,    660.,   -361.,
                      134.,   -184.,    328.,    264.,      5.,    -86.,    -16.,
                      63.,     61.,    -11.,   -217.,    -58.,     59.,    -90.,
                      70.,    -55.,      0.,     34.,    -41.,    -21.,     18.,
                      6.,     11.,      8.,     -4.,     -9.,      1.,      2.,
                      -9.,      5.,      8.,      8.,     10.,      1.,    -11.,
                      12.,      1.,     -2.,      2.,     -1.,     -1.,     -3.,
                      -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,
                      2.,      2.,      0.])

    H1900 = np.array([0.,      0.,   5922.,      0.,  -1061.,   1121.,      0.,
                      -330.,      3.,    523.,      0.,    195.,    -69.,   -210.,
                      -75.,      0.,   -210.,     53.,    -33.,   -124.,      3.,
                      0.,     -9.,     83.,      2.,    -35.,     36.,    -69.,
                      0.,    -45.,    -13.,    -10.,     -1.,     28.,    -12.,
                      -22.,      0.,      8.,    -14.,      7.,    -13.,      5.,
                      16.,     -5.,    -18.,      0.,    -20.,     14.,      5.,
                      -3.,     -2.,      8.,     10.,     -2.,      2.,      0.,
                      2.,      1.,      2.,      6.,     -4.,      0.,     -2.,
                      4.,      0.,     -6.])
            
    G1905 = np.array([0., -31464.,  -2298.,   -728.,   2928.,   1041.,   1037.,
                      -1494.,   1239.,    635.,    880.,    643.,    653.,   -380.,
                      146.,   -192.,    328.,    259.,     -1.,    -93.,    -26.,
                      62.,     60.,    -11.,   -221.,    -57.,     57.,    -92.,
                      70.,    -54.,      0.,     33.,    -41.,    -20.,     18.,
                      6.,     11.,      8.,     -4.,     -9.,      1.,      2.,
                      -8.,      5.,      8.,      8.,     10.,      1.,    -11.,
                      12.,      1.,     -2.,      2.,      0.,     -1.,     -3.,
                      -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,
                      2.,      2.,      0.])
            
    H1905 = np.array([0.,      0.,   5909.,      0.,  -1086.,   1065.,      0.,
                      -357.,     34.,    480.,      0.,    203.,    -77.,   -201.,
                      -65.,      0.,   -193.,     56.,    -32.,   -125.,     11.,
                      0.,     -7.,     86.,      4.,    -32.,     32.,    -67.,
                      0.,    -46.,    -14.,    -11.,      0.,     28.,    -12.,
                      -22.,      0.,      8.,    -15.,      7.,    -13.,      5.,
                      16.,     -5.,    -18.,      0.,    -20.,     14.,      5.,
                      -3.,     -2.,      8.,     10.,     -2.,      2.,      0.,
                      2.,      1.,      2.,      6.,     -4.,      0.,     -2.,
                      4.,      0.,     -6.])
            
    G1910 = np.array([0., -31354.,  -2297.,   -769.,   2948.,   1176.,   1058.,
                      -1524.,   1223.,    705.,    884.,    660.,    644.,   -400.,
                      160.,   -201.,    327.,    253.,     -9.,   -102.,    -38.,
                      62.,     58.,    -11.,   -224.,    -54.,     54.,    -95.,
                      71.,    -54.,      1.,     32.,    -40.,    -19.,     18.,
                      6.,     11.,      8.,     -4.,     -9.,      1.,      2.,
                      -8.,      5.,      8.,      8.,     10.,      1.,    -11.,
                      12.,      1.,     -2.,      2.,      0.,     -1.,     -3.,
                      -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,
                      2.,      2.,      0.])
            
    H1910 = np.array([0.,      0.,   5898.,      0.,  -1128.,   1000.,      0.,
                      -389.,     62.,    425.,      0.,    211.,    -90.,   -189.,
                      -55.,      0.,   -172.,     57.,    -33.,   -126.,     21.,
                      0.,     -5.,     89.,      5.,    -29.,     28.,    -65.,
                      0.,    -47.,    -14.,    -12.,      1.,     28.,    -13.,
                      -22.,      0.,      8.,    -15.,      6.,    -13.,      5.,
                      16.,     -5.,    -18.,      0.,    -20.,     14.,      5.,
                      -3.,     -2.,      8.,     10.,     -2.,      2.,      0.,
                      2.,      1.,      2.,      6.,     -4.,      0.,     -2.,
                      4.,      0.,     -6.])
            
    G1915 = np.array([0., -31212.,  -2306.,   -802.,   2956.,   1309.,   1084.,
                      -1559.,   1212.,    778.,    887.,    678.,    631.,   -416.,
                      178.,   -211.,    327.,    245.,    -16.,   -111.,    -51.,
                      61.,     57.,    -10.,   -228.,    -51.,     49.,    -98.,
                      72.,    -54.,      2.,     31.,    -38.,    -18.,     19.,
                      6.,     11.,      8.,     -4.,     -9.,      2.,      3.,
                      -8.,      6.,      8.,      8.,     10.,      1.,    -11.,
                      12.,      1.,     -2.,      2.,      0.,     -1.,     -3.,
                      -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,
                      1.,      2.,      0.])
            
    H1915 = np.array([0.,      0.,   5875.,      0.,  -1191.,    917.,      0.,
                      -421.,     84.,    360.,      0.,    218.,   -109.,   -173.,
                      -51.,      0.,   -148.,     58.,    -34.,   -126.,     32.,
                      0.,     -2.,     93.,      8.,    -26.,     23.,    -62.,
                      0.,    -48.,    -14.,    -12.,      2.,     28.,    -15.,
                      -22.,      0.,      8.,    -15.,      6.,    -13.,      5.,
                      16.,     -5.,    -18.,      0.,    -20.,     14.,      5.,
                      -3.,     -2.,      8.,     10.,     -2.,      2.,      0.,
                      2.,      1.,      2.,      6.,     -4.,      0.,     -2.,
                      4.,      0.,     -6.])
            
    G1920 = np.array([0., -31060.,  -2317.,   -839.,   2959.,   1407.,   1111.,
                      -1600.,   1205.,    839.,    889.,    695.,    616.,   -424.,
                      199.,   -221.,    326.,    236.,    -23.,   -119.,    -62.,
                      61.,     55.,    -10.,   -233.,    -46.,     44.,   -101.,
                      73.,    -54.,      2.,     29.,    -37.,    -16.,     19.,
                      6.,     11.,      7.,     -3.,     -9.,      2.,      4.,
                      -7.,      6.,      8.,      8.,     10.,      1.,    -11.,
                      12.,      1.,     -2.,      2.,      0.,     -1.,     -3.,
                      -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,
                      1.,      3.,      0.])
            
    H1920 = np.array([0.,      0.,   5845.,      0.,  -1259.,    823.,      0.,
                      -445.,    103.,    293.,      0.,    220.,   -134.,   -153.,
                      -57.,      0.,   -122.,     58.,    -38.,   -125.,     43.,
                      0.,      0.,     96.,     11.,    -22.,     18.,    -57.,
                      0.,    -49.,    -14.,    -13.,      4.,     28.,    -16.,
                      -22.,      0.,      8.,    -15.,      6.,    -14.,      5.,
                      17.,     -5.,    -19.,      0.,    -20.,     14.,      5.,
                      -3.,     -2.,      9.,     10.,     -2.,      2.,      0.,
                      2.,      1.,      2.,      6.,     -4.,      0.,     -2.,
                      4.,      0.,     -6.])

    G1925 = np.array([0., -30926.,  -2318.,   -893.,   2969.,   1471.,   1140.,
                      -1645.,   1202.,    881.,    891.,    711.,    601.,   -426.,
                      217.,   -230.,    326.,    226.,    -28.,   -125.,    -69.,
                      61.,     54.,     -9.,   -238.,    -40.,     39.,   -103.,
                      73.,    -54.,      3.,     27.,    -35.,    -14.,     19.,
                      6.,     11.,      7.,     -3.,     -9.,      2.,      4.,
                      -7.,      7.,      8.,      8.,     10.,      1.,    -11.,
                      12.,      1.,     -2.,      2.,      0.,     -1.,     -3.,
                      -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,
                      1.,      3.,      0.])
            
    H1925 = np.array([0.,      0.,   5817.,      0.,  -1334.,    728.,      0.,
                      -462.,    119.,    229.,      0.,    216.,   -163.,   -130.,
                      -70.,      0.,    -96.,     58.,    -44.,   -122.,     51.,
                      0.,      3.,     99.,     14.,    -18.,     13.,    -52.,
                      0.,    -50.,    -14.,    -14.,      5.,     29.,    -17.,
                      -21.,      0.,      8.,    -15.,      6.,    -14.,      5.,
                      17.,     -5.,    -19.,      0.,    -20.,     14.,      5.,
                      -3.,     -2.,      9.,     10.,     -2.,      2.,      0.,
                      2.,      1.,      2.,      6.,     -4.,      0.,     -2.,
                      4.,      0.,     -6.])

    G1930 = np.array([0., -30805.,  -2316.,   -951.,   2980.,   1517.,   1172.,
                      -1692.,   1205.,    907.,    896.,    727.,    584.,   -422.,
                      234.,   -237.,    327.,    218.,    -32.,   -131.,    -74.,
                      60.,     53.,     -9.,   -242.,    -32.,     32.,   -104.,
                      74.,    -54.,      4.,     25.,    -34.,    -12.,     18.,
                      6.,     11.,      7.,     -3.,     -9.,      2.,      5.,
                      -6.,      8.,      8.,      8.,     10.,      1.,    -12.,
                      12.,      1.,     -2.,      3.,      0.,     -2.,     -3.,
                      -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,
                      1.,      3.,      0.])

    H1930 = np.array([0.,      0.,   5808.,      0.,  -1424.,    644.,      0.,
                      -480.,    133.,    166.,      0.,    205.,   -195.,   -109.,
                      -90.,      0.,    -72.,     60.,    -53.,   -118.,     58.,
                      0.,      4.,    102.,     19.,    -16.,      8.,    -46.,
                      0.,    -51.,    -15.,    -14.,      6.,     29.,    -18.,
                      -20.,      0.,      8.,    -15.,      5.,    -14.,      5.,
                      18.,     -5.,    -19.,      0.,    -20.,     14.,      5.,
                      -3.,     -2.,      9.,     10.,     -2.,      2.,      0.,
                      2.,      1.,      2.,      6.,     -4.,      0.,     -2.,
                      4.,      0.,     -6.])

    G1935 = np.array([0., -30715.,  -2306.,  -1018.,   2984.,   1550.,   1206.,
                      -1740.,   1215.,    918.,    903.,    744.,    565.,   -415.,
                      249.,   -241.,    329.,    211.,    -33.,   -136.,    -76.,
                      59.,     53.,     -8.,   -246.,    -25.,     25.,   -106.,
                      74.,    -53.,      4.,     23.,    -33.,    -11.,     18.,
                      6.,     11.,      7.,     -3.,     -9.,      1.,      6.,
                      -6.,      8.,      7.,      8.,     10.,      1.,    -12.,
                      11.,      1.,     -2.,      3.,      0.,     -2.,     -3.,
                      -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,
                      2.,      3.,      0.])

    H1935 = np.array([0.,      0.,   5812.,      0.,  -1520.,    586.,      0.,
                      -494.,    146.,    101.,      0.,    188.,   -226.,    -90.,
                      -114.,      0.,    -51.,     64.,    -64.,   -115.,     64.,
                      0.,      4.,    104.,     25.,    -15.,      4.,    -40.,
                      0.,    -52.,    -17.,    -14.,      7.,     29.,    -19.,
                      -19.,      0.,      8.,    -15.,      5.,    -15.,      5.,
                      18.,     -5.,    -19.,      0.,    -20.,     15.,      5.,
                      -3.,     -3.,      9.,     11.,     -2.,      2.,      0.,
                      2.,      1.,      2.,      6.,     -4.,      0.,     -1.,
                      4.,      0.,     -6.])

    G1940 = np.array([0., -30654.,  -2292.,  -1106.,   2981.,   1566.,   1240.,
                      -1790.,   1232.,    916.,    914.,    762.,    550.,   -405.,
                      265.,   -241.,    334.,    208.,    -33.,   -141.,    -76.,
                      57.,     54.,     -7.,   -249.,    -18.,     18.,   -107.,
                      74.,    -53.,      4.,     20.,    -31.,     -9.,     17.,
                      5.,     11.,      7.,     -3.,    -10.,      1.,      6.,
                      -5.,      9.,      7.,      8.,     10.,      1.,    -12.,
                      11.,      1.,     -2.,      3.,      1.,     -2.,     -3.,
                      -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,
                      2.,      3.,      0.])

    H1940 = np.array([0.,      0.,   5821.,      0.,  -1614.,    528.,      0.,
                      -499.,    163.,     43.,      0.,    169.,   -252.,    -72.,
                      -141.,      0.,    -33.,     71.,    -75.,   -113.,     69.,
                      0.,      4.,    105.,     33.,    -15.,      0.,    -33.,
                      0.,    -52.,    -18.,    -14.,      7.,     29.,    -20.,
                      -19.,      0.,      8.,    -14.,      5.,    -15.,      5.,
                      19.,     -5.,    -19.,      0.,    -21.,     15.,      5.,
                      -3.,     -3.,      9.,     11.,     -2.,      2.,      0.,
                      2.,      1.,      2.,      6.,     -4.,      0.,     -1.,
                      4.,      0.,     -6.])

    G1945 = np.array([0., -30594.,  -2285.,  -1244.,   2990.,   1578.,   1282.,
                      -1834.,   1255.,    913.,    944.,    776.,    544.,   -421.,
                      304.,   -253.,    346.,    194.,    -20.,   -142.,    -82.,
                      59.,     57.,      6.,   -246.,    -25.,     21.,   -104.,
                      70.,    -40.,      0.,      0.,    -29.,    -10.,     15.,
                      29.,     13.,      7.,     -8.,     -5.,      9.,      7.,
                      -10.,      7.,      2.,      5.,    -21.,      1.,    -11.,
                      3.,     16.,     -3.,     -4.,     -3.,     -4.,     -3.,
                      11.,      1.,      2.,     -5.,     -1.,      8.,     -1.,
                      -3.,      5.,     -2.])

    H1945 = np.array([0.,      0.,   5810.,      0.,  -1702.,    477.,      0.,
                      -499.,    186.,    -11.,      0.,    144.,   -276.,    -55.,
                      -178.,      0.,    -12.,     95.,    -67.,   -119.,     82.,
                      0.,      6.,    100.,     16.,     -9.,    -16.,    -39.,
                      0.,    -45.,    -18.,      2.,      6.,     28.,    -17.,
                      -22.,      0.,     12.,    -21.,    -12.,     -7.,      2.,
                      18.,      3.,    -11.,      0.,    -27.,     17.,     29.,
                      -9.,      4.,      9.,      6.,      1.,      8.,      0.,
                      5.,      1.,    -20.,     -1.,     -6.,      6.,     -4.,
                      -2.,      0.,     -2.])

    G1950 = np.array([0., -30554.,  -2250.,  -1341.,   2998.,   1576.,   1297.,
                      -1889.,   1274.,    896.,    954.,    792.,    528.,   -408.,
                      303.,   -240.,    349.,    211.,    -20.,   -147.,    -76.,
                      54.,     57.,      4.,   -247.,    -16.,     12.,   -105.,
                      65.,    -55.,      2.,      1.,    -40.,     -7.,      5.,
                      19.,     22.,     15.,     -4.,     -1.,     11.,     15.,
                      -13.,      5.,     -1.,      3.,     -7.,     -1.,    -25.,
                      10.,      5.,     -5.,     -2.,      3.,      8.,     -8.,
                      4.,     -1.,     13.,     -4.,      4.,     12.,      3.,
                      2.,     10.,      3.])

    H1950 = np.array([0.,      0.,   5815.,      0.,  -1810.,    381.,      0.,
                      -476.,    206.,    -46.,      0.,    136.,   -278.,    -37.,
                      -210.,      0.,      3.,    103.,    -87.,   -122.,     80.,
                      0.,     -1.,     99.,     33.,    -12.,    -12.,    -30.,
                      0.,    -35.,    -17.,      0.,     10.,     36.,    -18.,
                      -16.,      0.,      5.,    -22.,      0.,    -21.,     -8.,
                      17.,     -4.,    -17.,      0.,    -24.,     19.,     12.,
                      2.,      2.,      8.,      8.,    -11.,     -7.,      0.,
                      13.,     -2.,    -10.,      2.,     -3.,      6.,     -3.,
                      6.,     11.,      8.])

    G1955 = np.array([0., -30500.,  -2215.,  -1440.,   3003.,   1581.,   1302.,
                      -1944.,   1288.,    882.,    958.,    796.,    510.,   -397.,
                      290.,   -229.,    360.,    230.,    -23.,   -152.,    -69.,
                      47.,     57.,      3.,   -247.,     -8.,      7.,   -107.,
                      65.,    -56.,      2.,     10.,    -32.,    -11.,      9.,
                      18.,     11.,      9.,     -6.,    -14.,      6.,     10.,
                      -7.,      6.,      9.,      4.,      9.,     -4.,     -5.,
                      2.,      4.,      1.,      2.,      2.,      5.,     -3.,
                      -5.,     -1.,      2.,     -3.,      7.,      4.,     -2.,
                      6.,     -2.,      0.])

    H1955 = np.array([0.,      0.,   5820.,      0.,  -1898.,    291.,      0.,
                      -462.,    216.,    -83.,      0.,    133.,   -274.,    -23.,
                      -230.,      0.,     15.,    110.,    -98.,   -121.,     78.,
                      0.,     -9.,     96.,     48.,    -16.,    -12.,    -24.,
                      0.,    -50.,    -24.,     -4.,      8.,     28.,    -20.,
                      -18.,      0.,     10.,    -15.,      5.,    -23.,      3.,
                      23.,     -4.,    -13.,      0.,    -11.,     12.,      7.,
                      6.,     -2.,     10.,      7.,     -6.,      5.,      0.,
                      -4.,      0.,     -8.,     -2.,     -4.,      1.,     -3.,
                      7.,     -1.,     -3.])

    G1960 = np.array([0., -30421.,  -2169.,  -1555.,   3002.,   1590.,   1302.,
                      -1992.,   1289.,    878.,    957.,    800.,    504.,   -394.,
                      269.,   -222.,    362.,    242.,    -26.,   -156.,    -63.,
                      46.,     58.,      1.,   -237.,     -1.,     -2.,   -113.,
                      67.,    -56.,      5.,     15.,    -32.,     -7.,     17.,
                      8.,     15.,      6.,     -4.,    -11.,      2.,     10.,
                      -5.,     10.,      8.,      4.,      6.,      0.,     -9.,
                      1.,      4.,     -1.,     -2.,      3.,     -1.,      1.,
                      -3.,      4.,      0.,     -1.,      4.,      6.,      1.,
                      -1.,      2.,      0.])

    H1960 = np.array([0.,      0.,   5791.,      0.,  -1967.,    206.,      0.,
                      -414.,    224.,   -130.,      0.,    135.,   -278.,      3.,
                      -255.,      0.,     16.,    125.,   -117.,   -114.,     81.,
                      0.,    -10.,     99.,     60.,    -20.,    -11.,    -17.,
                      0.,    -55.,    -28.,     -6.,      7.,     23.,    -18.,
                      -17.,      0.,     11.,    -14.,      7.,    -18.,      4.,
                      23.,      1.,    -20.,      0.,    -18.,     12.,      2.,
                      0.,     -3.,      9.,      8.,      0.,      5.,      0.,
                      4.,      1.,      0.,      2.,     -5.,      1.,     -1.,
                      6.,      0.,     -7.])

    G1965 = np.array([0., -30334.,  -2119.,  -1662.,   2997.,   1594.,   1297.,
                      -2038.,   1292.,    856.,    957.,    804.,    479.,   -390.,
                      252.,   -219.,    358.,    254.,    -31.,   -157.,    -62.,
                      45.,     61.,      8.,   -228.,      4.,      1.,   -111.,
                      75.,    -57.,      4.,     13.,    -26.,     -6.,     13.,
                      1.,     13.,      5.,     -4.,    -14.,      0.,      8.,
                      -1.,     11.,      4.,      8.,     10.,      2.,    -13.,
                      10.,     -1.,     -1.,      5.,      1.,     -2.,     -2.,
                      -3.,      2.,     -5.,     -2.,      4.,      4.,      0.,
                      2.,      2.,      0.])

    H1965 = np.array([0.,      0.,   5776.,      0.,  -2016.,    114.,      0.,
                      -404.,    240.,   -165.,      0.,    148.,   -269.,     13.,
                      -269.,      0.,     19.,    128.,   -126.,    -97.,     81.,
                      0.,    -11.,    100.,     68.,    -32.,     -8.,     -7.,
                      0.,    -61.,    -27.,     -2.,      6.,     26.,    -23.,
                      -12.,      0.,      7.,    -12.,      9.,    -16.,      4.,
                      24.,     -3.,    -17.,      0.,    -22.,     15.,      7.,
                      -4.,     -5.,     10.,     10.,     -4.,      1.,      0.,
                      2.,      1.,      2.,      6.,     -4.,      0.,     -2.,
                      3.,      0.,     -6.])

    G1970 = np.array([0., -30220.,  -2068.,  -1781.,   3000.,   1611.,   1287.,
                      -2091.,   1278.,    838.,    952.,    800.,    461.,   -395.,
                      234.,   -216.,    359.,    262.,    -42.,   -160.,    -56.,
                      43.,     64.,     15.,   -212.,      2.,      3.,   -112.,
                      72.,    -57.,      1.,     14.,    -22.,     -2.,     13.,
                      -2.,     14.,      6.,     -2.,    -13.,     -3.,      5.,
                      0.,     11.,      3.,      8.,     10.,      2.,    -12.,
                      10.,     -1.,      0.,      3.,      1.,     -1.,     -3.,
                      -3.,      2.,     -5.,     -1.,      6.,      4.,      1.,
                      0.,      3.,     -1.])

    H1970 = np.array([0.,      0.,   5737.,      0.,  -2047.,     25.,      0.,
                      -366.,    251.,   -196.,      0.,    167.,   -266.,     26.,
                      -279.,      0.,     26.,    139.,   -139.,    -91.,     83.,
                      0.,    -12.,    100.,     72.,    -37.,     -6.,      1.,
                      0.,    -70.,    -27.,     -4.,      8.,     23.,    -23.,
                      -11.,      0.,      7.,    -15.,      6.,    -17.,      6.,
                      21.,     -6.,    -16.,      0.,    -21.,     16.,      6.,
                      -4.,     -5.,     10.,     11.,     -2.,      1.,      0.,
                      1.,      1.,      3.,      4.,     -4.,      0.,     -1.,
                      3.,      1.,     -4.])

    G1975 = np.array([0., -30100.,  -2013.,  -1902.,   3010.,   1632.,   1276.,
                      -2144.,   1260.,    830.,    946.,    791.,    438.,   -405.,
                      216.,   -218.,    356.,    264.,    -59.,   -159.,    -49.,
                      45.,     66.,     28.,   -198.,      1.,      6.,   -111.,
                      71.,    -56.,      1.,     16.,    -14.,      0.,     12.,
                      -5.,     14.,      6.,     -1.,    -12.,     -8.,      4.,
                      0.,     10.,      1.,      7.,     10.,      2.,    -12.,
                      10.,     -1.,     -1.,      4.,      1.,     -2.,     -3.,
                      -3.,      2.,     -5.,     -2.,      5.,      4.,      1.,
                      0.,      3.,     -1.])

    H1975 = np.array([0.,      0.,   5675.,      0.,  -2067.,    -68.,      0.,
                      -333.,    262.,   -223.,      0.,    191.,   -265.,     39.,
                      -288.,      0.,     31.,    148.,   -152.,    -83.,     88.,
                      0.,    -13.,     99.,     75.,    -41.,     -4.,     11.,
                      0.,    -77.,    -26.,     -5.,     10.,     22.,    -23.,
                      -12.,      0.,      6.,    -16.,      4.,    -19.,      6.,
                      18.,    -10.,    -17.,      0.,    -21.,     16.,      7.,
                      -4.,     -5.,     10.,     11.,     -3.,      1.,      0.,
                      1.,      1.,      3.,      4.,     -4.,     -1.,     -1.,
                      3.,      1.,     -5.])

    G1980 = np.array([0., -29992.,  -1956.,  -1997.,   3027.,   1663.,   1281.,
                      -2180.,   1251.,    833.,    938.,    782.,    398.,   -419.,
                      199.,   -218.,    357.,    261.,    -74.,   -162.,    -48.,
                      48.,     66.,     42.,   -192.,      4.,     14.,   -108.,
                      72.,    -59.,      2.,     21.,    -12.,      1.,     11.,
                      -2.,     18.,      6.,      0.,    -11.,     -7.,      4.,
                      3.,      6.,     -1.,      5.,     10.,      1.,    -12.,
                      9.,     -3.,     -1.,      7.,      2.,     -5.,     -4.,
                      -4.,      2.,     -5.,     -2.,      5.,      3.,      1.,
                      2.,      3.,      0.])

    H1980 = np.array([0.,      0.,   5604.,      0.,  -2129.,   -200.,      0.,
                      -336.,    271.,   -252.,      0.,    212.,   -257.,     53.,
                      -297.,      0.,     46.,    150.,   -151.,    -78.,     92.,
                      0.,    -15.,     93.,     71.,    -43.,     -2.,     17.,
                      0.,    -82.,    -27.,     -5.,     16.,     18.,    -23.,
                      -10.,      0.,      7.,    -18.,      4.,    -22.,      9.,
                      16.,    -13.,    -15.,      0.,    -21.,     16.,      9.,
                      -5.,     -6.,      9.,     10.,     -6.,      2.,      0.,
                      1.,      0.,      3.,      6.,     -4.,      0.,     -1.,
                      4.,      0.,     -6.])

    G1985 = np.array([0., -29873.,  -1905.,  -2072.,   3044.,   1687.,   1296.,
                      -2208.,   1247.,    829.,    936.,    780.,    361.,   -424.,
                      170.,   -214.,    355.,    253.,    -93.,   -164.,    -46.,
                      53.,     65.,     51.,   -185.,      4.,     16.,   -102.,
                      74.,    -62.,      3.,     24.,     -6.,      4.,     10.,
                      0.,     21.,      6.,      0.,    -11.,     -9.,      4.,
                      4.,      4.,     -4.,      5.,     10.,      1.,    -12.,
                      9.,     -3.,     -1.,      7.,      1.,     -5.,     -4.,
                      -4.,      3.,     -5.,     -2.,      5.,      3.,      1.,
                      2.,      3.,      0.])


    H1985 = np.array([0.,      0.,   5500.,      0.,  -2197.,   -306.,      0.,
                      -310.,    284.,   -297.,      0.,    232.,   -249.,     69.,
                      -297.,      0.,     47.,    150.,   -154.,    -75.,     95.,
                      0.,    -16.,     88.,     69.,    -48.,     -1.,     21.,
                      0.,    -83.,    -27.,     -2.,     20.,     17.,    -23.,
                      -7.,      0.,      8.,    -19.,      5.,    -23.,     11.,
                      14.,    -15.,    -11.,      0.,    -21.,     15.,      9.,
                      -6.,     -6.,      9.,      9.,     -7.,      2.,      0.,
                      1.,      0.,      3.,      6.,     -4.,      0.,     -1.,
                      4.,      0.,     -6.])

    G1990 = np.array([0., -29775.,  -1848.,  -2131.,   3059.,   1686.,   1314.,
                      -2239.,   1248.,    802.,    939.,    780.,    325.,   -423.,
                      141.,   -214.,    353.,    245.,   -109.,   -165.,    -36.,
                      61.,     65.,     59.,   -178.,      3.,     18.,    -96.,
                      77.,    -64.,      2.,     26.,     -1.,      5.,      9.,
                      0.,     23.,      5.,     -1.,    -10.,    -12.,      3.,
                      4.,      2.,     -6.,      4.,      9.,      1.,    -12.,
                      9.,     -4.,     -2.,      7.,      1.,     -6.,     -3.,
                      -4.,      2.,     -5.,     -2.,      4.,      3.,      1.,
                      3.,      3.,      0.])

    H1990 = np.array([0.,      0.,   5406.,      0.,  -2279.,   -373.,      0.,
                      -284.,    293.,   -352.,      0.,    247.,   -240.,     84.,
                      -299.,      0.,     46.,    154.,   -153.,    -69.,     97.,
                      0.,    -16.,     82.,     69.,    -52.,      1.,     24.,
                      0.,    -80.,    -26.,      0.,     21.,     17.,    -23.,
                      -4.,      0.,     10.,    -19.,      6.,    -22.,     12.,
                      12.,    -16.,    -10.,      0.,    -20.,     15.,     11.,
                      -7.,     -7.,      9.,      8.,     -7.,      2.,      0.,
                      2.,      1.,      3.,      6.,     -4.,      0.,     -2.,
                      3.,     -1.,     -6.])

    G1995 = np.array([0., -29692.,  -1784.,  -2200.,   3070.,   1681.,   1335.,
                      -2267.,   1249.,    759.,    940.,    780.,    290.,   -418.,
                      122.,   -214.,    352.,    235.,   -118.,   -166.,    -17.,
                      68.,     67.,     68.,   -170.,     -1.,     19.,    -93.,
                      77.,    -72.,      1.,     28.,      5.,      4.,      8.,
                      -2.,     25.,      6.,     -6.,     -9.,    -14.,      9.,
                      6.,     -5.,     -7.,      4.,      9.,      3.,    -10.,
                      8.,     -8.,     -1.,     10.,     -2.,     -8.,     -3.,
                      -6.,      2.,     -4.,     -1.,      4.,      2.,      2.,
                      5.,      1.,      0.])

    H1995 = np.array([0.,      0.,   5306.,      0.,  -2366.,   -413.,      0.,
                      -262.,    302.,   -427.,      0.,    262.,   -236.,     97.,
                      -306.,      0.,     46.,    165.,   -143.,    -55.,    107.,
                      0.,    -17.,     72.,     67.,    -58.,      1.,     36.,
                      0.,    -69.,    -25.,      4.,     24.,     17.,    -24.,
                      -6.,      0.,     11.,    -21.,      8.,    -23.,     15.,
                      11.,    -16.,     -4.,      0.,    -20.,     15.,     12.,
                      -6.,     -8.,      8.,      5.,     -8.,      3.,      0.,
                      1.,      0.,      4.,      5.,     -5.,     -1.,     -2.,
                      1.,     -2.,     -7.])

    G2000 = np.array([0.0,-29619.4, -1728.2, -2267.7,  3068.4,  1670.9,  1339.6,
                      -2288.0,  1252.1,   714.5,   932.3,   786.8,   250.0,  -403.0,
                      111.3,  -218.8,   351.4,   222.3,  -130.4,  -168.6,   -12.9,
                      72.3,    68.2,    74.2,  -160.9,    -5.9,    16.9,   -90.4,
                      79.0,   -74.0,     0.0,    33.3,     9.1,     6.9,     7.3,
                      -1.2,    24.4,     6.6,    -9.2,    -7.9,   -16.6,     9.1,
                      7.0,    -7.9,    -7.0,     5.0,     9.4,     3.0,    -8.4,
                      6.3,    -8.9,    -1.5,     9.3,    -4.3,    -8.2,    -2.6,
                      -6.0,     1.7,    -3.1,    -0.5,     3.7,     1.0,     2.0,
                      4.2,     0.3,    -1.1])

    H2000 = np.array([0.0,     0.0,  5186.1,     0.0, -2481.6,  -458.0,     0.0,
                      -227.6,   293.4,  -491.1,     0.0,   272.6,  -231.9,   119.8,
                      -303.8,     0.0,    43.8,   171.9,  -133.1,   -39.3,   106.3,
                      0.0,   -17.4,    63.7,    65.1,   -61.2,     0.7,    43.8,
                      0.0,   -64.6,   -24.2,     6.2,    24.0,    14.8,   -25.4,
                      -5.8,     0.0,    11.9,   -21.5,     8.5,   -21.5,    15.5,
                      8.9,   -14.9,    -2.1,     0.0,   -19.7,    13.4,    12.5,
                      -6.2,    -8.4,     8.4,     3.8,    -8.2,     4.8,     0.0,
                      1.7,     0.0,     4.0,     4.9,    -5.9,    -1.2,    -2.9,
                      0.0,    -2.2,    -7.4])

    G2005 = np.array([0.00,-29554.63,-1669.05,-2337.24, 3047.69, 1657.76, 1336.30,
                      -2305.83,  1246.39,  672.51,  920.55,  797.96,  210.65, -379.86,
                      100.00,  -227.00,  354.41,  208.95, -136.54, -168.05,  -13.55,
                      73.60,    69.56,   76.74, -151.34,  -14.58,   14.58,  -86.36,
                      79.88,   -74.46,   -1.65,   38.73,   12.30,    9.37,    5.42,
                      1.94,    24.80,    7.62,  -11.73,   -6.88,  -18.11,   10.17,
                      9.36,   -11.25,   -4.87,    5.58,    9.76,    3.58,   -6.94,
                      5.01,   -10.76,   -1.25,    8.76,   -6.66,   -9.22,   -2.17,
                      -6.12,     1.42,   -2.35,   -0.15,    3.06,    0.29,    2.06,
                      3.77,    -0.21,   -2.09])

    H2005 = np.array([0.00,    0.00, 5077.99,     0.00,-2594.50, -515.43,    0.00,
                      -198.86,  269.72, -524.72,     0.00,  282.07, -225.23,  145.15,
                      -305.36,    0.00,   42.72,   180.25, -123.45,  -19.57,  103.85,
                      0.00,  -20.33,   54.75,    63.63,  -63.53,    0.24,   50.94,
                      0.00,  -61.14,  -22.57,     6.82,   25.35,   10.93,  -26.32,
                      -4.64,    0.00,   11.20,   -20.88,    9.83,  -19.71,   16.22,
                      7.61,  -12.76,   -0.06,     0.00,  -20.11,   12.69,   12.67,
                      -6.72,   -8.16,    8.10,     2.92,   -7.73,    6.01,    0.00,
                      2.19,    0.10,    4.46,     4.76,   -6.58,   -1.01,   -3.47,
                      -0.86,   -2.31,   -7.93])

    G2010 = np.array([0.0,-29496.6, -1586.4,  -2396.1,  3026.3,  1668.2,  1339.8,
                      -2326.5,  1232.1,   633.7,    912.7,   809.0,   166.6,  -356.8,
                      89.4,  -230.9,   357.3,    200.3,  -141.1,  -163.2,    -8.0,
                      72.8,    68.7,    75.9,   -141.4,   -22.8,    13.1,   -78.1,
                      80.4,   -75.0,    -4.6,     45.2,    14.0,    10.5,     1.6,
                      4.9,    24.4,     8.2,    -14.5,    -5.6,   -19.3,    11.6,
                      10.9,   -14.1,    -3.5,      5.5,     9.4,     3.5,    -5.3,
                      3.1,   -12.4,    -0.8,      8.4,    -8.4,   -10.1,    -1.9,
                      -6.2,     0.9,    -1.1,     -0.2,     2.5,    -0.3,     2.1,
                      3.1,    -1.0,    -2.8])

    H2010 = np.array([0.0,     0.0,  4944.3,      0.0, -2708.5,  -575.7,     0.0,
                      -160.4,   251.8,  -537.0,      0.0,   286.5,  -211.0,   164.5,
                      -309.7,     0.0,    44.6,    189.0,  -118.1,     0.0,   101.0,
                      0.0,   -20.9,    44.2,     61.5,   -66.3,     3.0,    55.4,
                      0.0,   -57.8,   -21.2,      6.5,    25.0,     7.0,   -27.6,
                      -3.3,     0.0,    10.8,    -20.0,    11.8,   -17.4,    16.7,
                      7.0,   -10.7,     1.6,      0.0,   -20.5,    11.5,    12.8,
                      -7.1,    -7.4,     8.0,      2.1,    -6.1,     7.0,     0.0,
                      2.7,    -0.1,     4.7,      4.4,    -7.2,    -1.0,    -4.0,
                      -2.0,    -2.0,    -8.3])

    G2015 = np.array([0.0,-29442.0, -1501.0,  -2445.1,  3012.9,  1676.7,  1350.7,
                      -2352.3,  1225.6,   582.0,    907.6,   813.7,   120.4,  -334.9,
                      70.4,  -232.6,   360.1,    192.4,  -140.9,  -157.5,     4.1,
                      70.0,    67.7,    72.7,   -129.9,   -28.9,    13.2,   -70.9,
                      81.6,   -76.1,    -6.8,     51.8,    15.0,     9.4,    -2.8,
                      6.8,    24.2,     8.8,    -16.9,    -3.2,   -20.6,    13.4,
                      11.7,   -15.9,    -2.0,      5.4,     8.8,     3.1,    -3.3,
                      0.7,   -13.3,    -0.1,      8.7,    -9.1,   -10.5,    -1.9,
                      -6.3,     0.1,     0.5,     -0.5,     1.8,    -0.7,     2.1,
                      2.4,    -1.8,    -3.6])


    H2015 = np.array([0.0,     0.0,  4797.1,      0.0, -2845.6,  -641.9,     0.0,
                      -115.3,   244.9,  -538.4,      0.0,   283.3,  -188.7,   180.9,
                      -329.5,     0.0,    47.3,    197.0,  -119.3,    16.0,   100.2,
                      0.0,   -20.8,    33.2,     58.9,   -66.7,     7.3,    62.6,
                      0.0,   -54.1,   -19.5,      5.7,    24.4,     3.4,   -27.4,
                      -2.2,     0.0,    10.1,    -18.3,    13.3,   -14.6,    16.2,
                      5.7,    -9.1,     2.1,      0.0,   -21.6,    10.8,    11.8,
                      -6.8,    -6.9,     7.8,      1.0,    -4.0,     8.4,     0.0,
                      3.2,    -0.4,     4.6,      4.4,    -7.9,    -0.6,    -4.2,
                      -2.8,    -1.2,    -8.7])

    DG = np.array([0.0,    10.3,    18.1,     -8.7,    -3.3,     2.1,     3.4,
                   -5.5,    -0.7,   -10.1,     -0.7,     0.2,    -9.1,     4.1,
                   -4.3,    -0.2,     0.5,     -1.3,    -0.1,     1.4,     3.9,
                   -0.3,    -0.1,    -0.7,      2.1,    -1.2,     0.3,     1.6,
                   0.3,    -0.2,    -0.5,      1.3,     0.1,    -0.6,    -0.8,
                   0.2,     0.2,     0.0,     -0.6,     0.5,    -0.2,     0.4,
                   0.1,    -0.4,     0.3])

    DH = np.array([0.0,     0.0,   -26.6,      0.0,   -27.4,   -14.1,     0.0,
                   8.2,    -0.4,     1.8,      0.0,    -1.3,     5.3,     2.9,
                   -5.2,     0.0,     0.6,      1.7,    -1.2,     3.4,     0.0,
                   0.0,     0.0,    -2.1,     -0.7,     0.2,     0.9,     1.0,
                   0.0,     0.8,     0.4,     -0.2,    -0.3,    -0.6,     0.1,
                   -0.2,     0.0,    -0.3,      0.3,     0.1,     0.5,    -0.2,
                   -0.3,     0.3,     0.0])

#
#
    MA = 0
    IYR = 0
    
    if (MA != 1): #GOTO 10
    #10    MA=1
        MA=1 
#
#    DO 20 N=1,11
        for N in range(1,11):
            N2=2*N-1
            N2=N2*(N2-2)
#        DO 20 M=1,N
            for M in range(1,N):
                MN = N*(N-1)/2+M
                REC[MN] = float((N-M)*(N+M-2))/float(N2)
#20    REC(MN)=FLOAT((N-M)*(N+M-2))/FLOAT(N2)
    elif (IY != IYR):
        pass #GOTO 30
    else: #GOTO 130
    
#
#30    IYR=IY
        ps = '   IGRF: GIVEN YEAR %5d IS OUT OF INTERVAL 1900-2015\n   *** CALCULATIONS WILL BE DONE FOR YEAR = %5d ***'

        IYR=IY
        if (IYR < 1900): IYR=1900
        if (IYR > 2020): IYR=2020
        if ((IY != IYR) and mess): print(ps,IY,IYR)  #WRITE (konsol,999)IY,IYR

#	include 'igrf_goto.h'

#       INTERPOLATE BETWEEN YEARS 
        if (IYR < 1905):# GOTO 1900      !INTERPOLATE BETWEEN 1900 - 1905
            #1900  F2=(IYR-1900)/5.
            F2=(IYR-1900)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1900[N]*F1+G1905[N]*F2
                H[N]=H1900[N]*F1+H1905[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1910): #GOTO 1905      !INTERPOLATE BETWEEN 1905 - 1910
            #1905  F2=(IYR-1905)/5.
            F2=(IYR-1905)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1905[N]*F1+G1910[N]*F2
                H[N]=H1905[N]*F1+H1910[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1915): #GOTO 1910      !INTERPOLATE BETWEEN 1910 - 1915
            #1910  F2=(IYR-1910)/5.
            F2=(IYR-1910)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1910[N]*F1+G1915[N]*F2
                H[N]=H1910[N]*F1+H1915[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1920): #GOTO 1915      !INTERPOLATE BETWEEN 1915 - 1920
            #1915  F2=(IYR-1915)/5.
            F2=(IYR-1915)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1915[N]*F1+G1920[N]*F2
                H[N]=H1915[N]*F1+H1920[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1925): #GOTO 1920      !INTERPOLATE BETWEEN 1920 - 1925
            #1920  F2=(IYR-1920)/5.
            F2=(IYR-1920)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1920[N]*F1+G1925[N]*F2
                H[N]=H1920[N]*F1+H1925[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1930): #GOTO 1925      !INTERPOLATE BETWEEN 1925 - 1930
            #1925  F2=(IYR-1925)/5.
            F2=(IYR-1925)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1925[N]*F1+G1930[N]*F2
                H[N]=H1925[N]*F1+H1930[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1935): #GOTO 1930      !INTERPOLATE BETWEEN 1930 - 1935
            #1930  F2=(IYR-1930)/5.
            F2=(IYR-1930)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1930[N]*F1+G1935[N]*F2
                H[N]=H1930[N]*F1+H1935[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1940): #GOTO 1935      !INTERPOLATE BETWEEN 1935 - 1940
            #1935  F2=(IYR-1935)/5.
            F2=(IYR-1935)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1935[N]*F1+G1940[N]*F2
                H[N]=H1935[N]*F1+H1940[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1945): #GOTO 1940      !INTERPOLATE BETWEEN 1940 - 1945
            #1940  F2=(IYR-1940)/5.
            F2=(IYR-1940)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1940[N]*F1+G1945[N]*F2
                H[N]=H1940[N]*F1+H1945[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1950): #GOTO 1945      !INTERPOLATE BETWEEN 1945 - 1950
            #1945  F2=(IYR-1945)/5.
            F2=(IYR-1945)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1945[N]*F1+G1950[N]*F2
                H[N]=H1945[N]*F1+H1950[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1955): #GOTO 1950      !INTERPOLATE BETWEEN 1950 - 1955
            #1950  F2=(IYR-1950)/5.
            F2=(IYR-1950)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1950[N]*F1+G1955[N]*F2
                H[N]=H1950[N]*F1+H1955[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1960): #GOTO 1955      !INTERPOLATE BETWEEN 1955 - 1960
            #1955  F2=(IYR-1955)/5.
            F2=(IYR-1955)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1955[N]*F1+G1960[N]*F2
                H[N]=H1955[N]*F1+H1960[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1965): #GOTO 1960      !INTERPOLATE BETWEEN 1960 - 1965
            #1960  F2=(IYR-1960)/5.
            F2=(IYR-1960)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1960[N]*F1+G1965[N]*F2
                H[N]=H1960[N]*F1+H1965[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1970): #GOTO 1965      !INTERPOLATE BETWEEN 1965 - 1970
            #1965  F2=(IYR-1965)/5.
            F2=(IYR-1965)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1965[N]*F1+G1970[N]*F2
                H[N]=H1965[N]*F1+H1970[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1975): #GOTO 1970      !INTERPOLATE BETWEEN 1970 - 1975
            #1970  F2=(IYR-1970)/5.
            F2=(IYR-1970)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1970[N]*F1+G1975[N]*F2
                H[N]=H1970[N]*F1+H1975[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1980): #GOTO 1975      !INTERPOLATE BETWEEN 1975 - 1980
            #1975  F2=(IYR-1975)/5.
            F2=(IYR-1975)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1975[N]*F1+G1980[N]*F2
                H[N]=H1975[N]*F1+H1980[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1985): #GOTO 1980      !INTERPOLATE BETWEEN 1980 - 1985:
            #1980  F2=(IYR-1980)/5.
            F2=(IYR-1980)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1980[N]*F1+G1985[N]*F2
                H[N]=H1980[N]*F1+H1985[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1990): #GOTO 1985      !INTERPOLATE BETWEEN 1985 - 1990:
            #1985  F2=(IYR-1985)/5.
            F2=(IYR-1985)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1985[N]*F1+G1990[N]*F2
                H[N]=H1985[N]*F1+H1990[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 1995): #GOTO 1990      !INTERPOLATE BETWEEN 1990 - 1995
            #1990  F2=(IYR-1990)/5.
            F2=(IYR-1990)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1990[N]*F1+G1995[N]*F2
                H[N]=H1990[N]*F1+H1995[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 2000): #GOTO 1995      !INTERPOLATE BETWEEN 1995 - 2000
            #1995  F2=(IYR-1995)/5.
            F2=(IYR-1995)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G1995[N]*F1+G2000[N]*F2
                H[N]=H1995[N]*F1+H2000[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 2005): #GOTO 2000      !INTERPOLATE BETWEEN 2000 - 2005
            #2000  F2=(IYR-2000)/5.
            F2=(IYR-2000)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G2000[N]*F1+G2005[N]*F2
                H[N]=H2000[N]*F1+H2005[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 2010): #GOTO 2005      !INTERPOLATE BETWEEN 2005 - 2010
            #2005  F2=(IYR-2005)/5.
            F2=(IYR-2005)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G2005[N]*F1+G2010[N]*F2
                H[N]=H2005[N]*F1+H2010[N]*F2
            #ENDDO
        #GOTO 300
        elif (IYR < 2015): #GOTO 2010      !INTERPOLATE BETWEEN 2010 - 2015:
            #2010  F2=(IYR-2010)/5.
            F2=(IYR-2010)/5.
            F1=1.-F2
            for N in range(1,66):
                G[N]=G2010[N]*F1+G2015[N]*F2
                H[N]=H2010[N]*F1+H2015[N]*F2
            #ENDDO
        #GOTO 300
        else:
#
#       EXTRAPOLATE BEYOND 2015:
#
            DT=float(IYR)-2015.
            for N in range(1,66):
                G[N]=G2015[N]
                H[N]=H2015[N]
                if (N <= 45):
                    G[N]=G[N]+DG[N]*DT
                    H[N]=H[N]+DH[N]*DT
                #40    CONTINUE
        #GOTO 300

#   COEFFICIENTS FOR A GIVEN YEAR HAVE BEEN CALCULATED; NOW MULTIPLY
#   THEM BY SCHMIDT NORMALIZATION FACTORS:

#300   S=1.
        S=1.
        for N in range(2,11):
            MN=N*(N-1)/2+1
            S=S*float(2*N-3)/float(N-1)
            G[MN]=G[MN]*S
            H[MN]=H[MN]*S
            P=S
            for M in range(2,N):
                AA=1.
                if (M == 2): AA=2.
                P=P*np.sqrt(AA*float(N-M+1)/float(N+M-2))
                MNN=MN+M-1
                G[MNN]=G[MNN]*P
                H[MNN]=H[MNN]*P   #120

#     NOW CALCULATE THE FIELD COMPONENTS
#     (IN CASE OF MULTIPLE INVOCATIONS WITH THE SAME VALUES OF IY AND NM,
#      CALCULATIONS START RIGHT HERE):

#130   PP=1./R
    PP=1./R
    P=PP
        
    K=NM+1
    for N in range(1,K):
        P=P*PP
        A[N]=P
        B[N]=P*N
#150      B(N)=P*N
    P=1.
    D=0.
    BBR=0.
    BBT=0.
    BBF=0.
    U=T
    CF=np.cos(F)
    SF=np.sin(F)
    C=np.cos(U)
    S=np.sin(U)
    for M in range(1,K):  #200
            #if (M.EQ.1) GOTO 160
        if (M != 1):
            MM=M-1
            W=X
            X=W*CF+Y*SF
            Y=Y*CF-W*SF
                #GOTO 170
        else:
                #160      X=0.
            X=0.
            Y=1.
            #170      Q=P
        Q=P
        Z=D
        BI=0.
        P2=0.
        D2=0.
        for N in range(M,K):
            AN=A[N]
            MN=N*(N-1)/2+M
            E=G[MN]
            HH=H[MN]
            W=E*Y+HH*X
            BBR=BBR+B[N]*W*Q
            BBT=BBT-AN*W*Z
                #IF(M.EQ.1) GOTO 180
            if (M != 1):
                QQ=Q
                if (S < 1.0E-5): QQ=Z
                BI=BI+AN*(E*X-HH*Y)*QQ
            XK=REC[MN]
                #180         XK=REC(MN)
            DP=C*Z-S*Q-XK*D2
            PM=C*Q-XK*P2
            D2=Z
            P2=Q
            Z=DP
            #190        Q=PM
            Q=PM
        D=S*D+C*P
        P=S*P
            #IF(M.EQ.1) GOTO 200
        if (M != 1):
            BI=BI*MM
            BBF=BBF+BI
        #200   CONTINUE
#
    BR=BBR
    BT=BBT
        #IF(S.LT.1.E-5) GOTO 210
    if (S >= 1.E-5): 
        BF=BBF/S
        return
#210   IF(C.LT.0.) BBF=-BBF
    if (C < 0.): BBF=-BBF
    BF=BBF
    return
#
#999   FORMAT(/
#     * '   IGRF: GIVEN YEAR',I5,' IS OUT OF INTERVAL 1900-2015'/,
#     * '   *** CALCULATIONS WILL BE DONE FOR YEAR =',I5,' ***'/)
#      END
#
#
def RECALC(IYR,IDAY,IHOUR,MIN,ISEC):
    '''*********************************************************************
        If only IYR is given then CALL RECALC(IYR,0,25,0,0)
        THIS IS A MODIFIED VERSION OF THE SUBROUTINE RECOMP WRITTEN BY
        N. A. TSYGANENKO. SINCE I WANT TO USE IT IN PLACE OF SUBROUTINE
        RECALC, I HAVE RENAMED THIS ROUTINE RECALC AND ELIMINATED THE
        ORIGINAL RECALC FROM THIS VERSION OF THE <GEOPACK.FOR> PACKAGE
        THIS WAY ALL ORIGINAL CALLS TO RECALC WILL CONTINUE TO WORK WITHOUT
        HAVING TO CHANGE THEM TO CALLS TO RECOMP.
        
        AN ALTERNATIVE VERSION OF THE SUBROUTINE RECALC FROM THE GEOPACK
        PACKAGE BASED ON A DIFFERENT APPROACH TO DERIVATION OF ROTATION
        MATRIX ELEMENTS
        
        THIS SUBROUTINE WORKS BY 20% FASTER THAN RECALC AND IS EASIER TO
        UNDERSTAND
        #####################################################
        #  WRITTEN BY  N.A. TSYGANENKO ON DECEMBER 1, 1991  #
        #####################################################
        Modified by Mauricio Peredo, Hughes STX at NASA/GSFC Code 695,
        September 1992
        
        Modified to accept years up to year 2000 and updated IGRF coeficients
        from 1945 (updated by V. Papitashvili, February 1995)
        
        Modified to accept years up to 2005 (V. Papitashvili, January 2001)
        
        Modified to accept years from 1900 through 2010 using the DGRF & 
        IGRF-10 coeficients (updated by V. Papitashvili, November 2005)
        
        Modified to accept years up to 2015 (V. Papitashvili, January 2011)
        
        Modified to accept years up to 2020 (D. Bilitza, October 2015)
        
        OTHER SUBROUTINES CALLED BY THIS ONE: SUN
        
        IYR = YEAR NUMBER (FOUR DIGITS)
        IDAY = DAY OF YEAR (DAY 1 = JAN 1)
        IHOUR = HOUR OF DAY (00 TO 23)
        MIN = MINUTE OF HOUR (00 TO 59)
        ISEC = SECONDS OF DAY(00 TO 59)
       *********************************************************************'''
#    IMPLICIT NONE
#
#        REAL ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,CPS,
#     1       SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,
#     2       A33,DS3,F2,F1,G10,G11,H11,DT,SQ,SQQ,SQR,S1,S2,
#     3       S3,CGST,SGST,DIP1,DIP2,DIP3,Y1,Y2,Y3,Y,Z1,Z2,Z3,DJ,
#     4       T,OBLIQ,DZ1,DZ2,DZ3,DY1,DY2,DY3,EXMAGX,EXMAGY,EXMAGZ,
#     5       EYMAGX,EYMAGY,GST,SLONG,SRASN,SDEC,BA(8),DECARG
#
#        INTEGER IYR,IDAY,IHOUR,MIN,ISEC,K,IY,IDE,IYE,konsol
#        logical	mess

    global ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,CPS,SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3,K,IY,BA
    global konsol,mess  

#      DATA IYE,IDE/2*0/
    IYE = 2
    IDE = 0
    if (IYR == IYE) and (IDAY == IDE): #GOTO 5
        pass
    else:

#  IYE AND IDE ARE THE CURRENT VALUES OF YEAR AND DAY NUMBER

        IY=IYR
        IDE=IDAY
        if (IY < 1900): IY=1900
#      IF(IY.GT.2015) IY=2015
        if (IY > 2020): IY=2020

#  WE ARE RESTRICTED BY THE INTERVAL 1900-2015, FOR WHICH THE DGRF & IGRF-11
#  COEFFICIENTS ARE KNOWN; IF IYR IS OUTSIDE THIS INTERVAL, THE
#  SUBROUTINE GIVES A WARNING (BUT DOES NOT REPEAT IT AT THE NEXT CALLS)

        #if ((IY != IYR) and mess) write(konsol,10) IYR,IY
        IYE=IY

#  LINEAR INTERPOLATION OF THE GEODIPOLE MOMENT COMPONENTS BETWEEN THE
#  VALUES FOR THE NEAREST EPOCHS:

        if (IY < 1905): #THEN                             !1900-1905
            F2=(float(IY)+float(IDAY)/365.-1900.)/5.
            F1=1.0-F2
            G10=31543.*F1+31464.*F2
            G11=-2298.*F1-2298.*F2
            H11= 5922.*F1+5909.*F2
        elif (IY < 1910): #THEN                         !1905-1910
            F2=(float(IY)+float(IDAY)/365.-1905.)/5.
            F1=1.0-F2
            G10=31464.*F1+31354.*F2
            G11=-2298.*F1-2297.*F2
            H11= 5909.*F1+5898.*F2
        elif (IY < 1915): #THEN                         !1910-1915
            F2=(float(IY)+float(IDAY)/365.-1910.)/5.
            F1=1.0-F2
            G10=31354.*F1+31212.*F2
            G11=-2297.*F1-2306.*F2
            H11= 5898.*F1+5875.*F2
        elif (IY < 1920): #THEN                         !1915-1920
            F2=(float(IY)+float(IDAY)/365.-1915.)/5.
            F1=1.0-F2
            G10=31212.*F1+31060.*F2
            G11=-2306.*F1-2317.*F2
            H11= 5875.*F1+5845.*F2
        elif (IY < 1925): #THEN                         !1920-1925
            F2=(float(IY)+float(IDAY)/365.-1920.)/5.
            F1=1.0-F2
            G10=31060.*F1+30926.*F2
            G11=-2317.*F1-2318.*F2
            H11= 5845.*F1+5817.*F2
        elif (IY < 1930): #THEN                         !1925-1930
            F2=(float(IY)+float(IDAY)/365.-1925.)/5.
            F1=1.0-F2
            G10=30926.*F1+30805.*F2
            G11=-2318.*F1-2316.*F2
            H11= 5817.*F1+5808.*F2
        elif (IY < 1935): #THEN                        !1930-1935
            F2=(float(IY)+float(IDAY)/365.-1930.)/5.
            F1=1.0-F2
            G10=30805.*F1+30715.*F2
            G11=-2316.*F1-2306.*F2
            H11= 5808.*F1+5812.*F2
        elif (IY < 1940): #THEN                        !1935-1940
            F2=(float(IY)+float(IDAY)/365.-1935.)/5.
            F1=1.0-F2
            G10=30715.*F1+30654.*F2
            G11=-2306.*F1-2292.*F2
            H11= 5812.*F1+5821.*F2
        elif (IY < 1945): #THEN                        !1940-1945
            F2=(float(IY)+float(IDAY)/365.-1940.)/5.
            F1=1.0-F2
            G10=30654.*F1+30594.*F2
            G11=-2292.*F1-2285.*F2
            H11= 5821.*F1+5810.*F2
        elif (IY < 1950): #THEN                        !1945-1950
            F2=(float(IY)+float(IDAY)/365.-1945.)/5.
            F1=1.0-F2
            G10=30594.*F1+30554.*F2
            G11=-2285.*F1-2250.*F2
            H11= 5810.*F1+5815.*F2
        elif (IY < 1955): #THEN                        !1950-1955
            F2=(float(IY)+float(IDAY)/365.-1950.)/5.
            F1=1.0-F2
            G10=30554.*F1+30500.*F2
            G11=-2250.*F1-2215.*F2
            H11= 5815.*F1+5820.*F2
        elif (IY < 1960): #THEN                        !1955-1960
            F2=(float(IY)+float(IDAY)/365.-1955.)/5.
            F1=1.0-F2
            G10=30500.*F1+30421.*F2
            G11=-2215.*F1-2169.*F2
            H11= 5820.*F1+5791.*F2
        elif (IY < 1965): #THEN                        !1960-1965
            F2=(float(IY)+float(IDAY)/365.-1960.)/5.
            F1=1.0-F2
            G10=30421.*F1+30334.*F2
            G11=-2169.*F1-2119.*F2
            H11= 5791.*F1+5776.*F2
        elif (IY < 1970): #THEN                        !1965-1970
            F2=(float(IY)+float(IDAY)/365.-1965.)/5.
            F1=1.0-F2
            G10=30334.*F1+30220.*F2
            G11=-2119.*F1-2068.*F2
            H11= 5776.*F1+5737.*F2
        elif (IY < 1975): #THEN                        !1970-1975
            F2=(float(IY)+float(IDAY)/365.-1970.)/5.
            F1=1.0-F2
            G10=30220.*F1+30100.*F2
            G11=-2068.*F1-2013.*F2
            H11= 5737.*F1+5675.*F2
        elif (IY < 1980): #THEN                        !1975-1980
            F2=(float(IY)+float(IDAY)/365.-1975.)/5.
            F1=1.0-F2
            G10=30100.*F1+29992.*F2
            G11=-2013.*F1-1956.*F2
            H11= 5675.*F1+5604.*F2
        elif (IY < 1985): #THEN                        !1980-1985
            F2=(float(IY)+float(IDAY)/365.-1980.)/5.
            F1=1.0-F2
            G10=29992.*F1+29873.*F2
            G11=-1956.*F1-1905.*F2
            H11= 5604.*F1+5500.*F2
        elif (IY < 1990): #THEN                        !1985-1990
            F2=(float(IY)+float(IDAY)/365.-1985.)/5.
            F1=1.0-F2
            G10=29873.*F1+29775.*F2
            G11=-1905.*F1-1848.*F2
            H11= 5500.*F1+5406.*F2
        elif (IY < 1995): #THEN                        !1990-1995
            F2=(float(IY)+float(IDAY)/365.-1990.)/5.
            F1=1.0-F2
            G10=29775.*F1+29692.*F2
            G11=-1848.*F1-1784.*F2
            H11= 5406.*F1+5306.*F2
        elif (IY < 2000): #THEN                        !1995-2000
            F2=(float(IY)+float(IDAY)/365.-1995.)/5.
            F1=1.0-F2
            G10=29692.*F1+29619.4*F2
            G11=-1784.*F1-1728.2*F2
            H11= 5306.*F1+5186.1*F2
        elif (IY < 2005): #THEN                        !2000-2005
            F2=(float(IY)+float(IDAY)/365.-2000.)/5.
            F1=1.0-F2
            G10=29619.4*F1+29554.63*F2
            G11=-1728.2*F1-1669.05*F2
            H11= 5186.1*F1+5077.99*F2
        elif (IY < 2010): #THEN                        !2005-2010
            F2=(float(IY)+float(IDAY)/365.-2005.)/5.
            F1=1.0-F2
            G10=29554.63*F1+29496.57*F2
            G11=-1669.05*F1-1586.42*F2
            H11= 5077.99*F1+4944.26*F2
        elif (IY < 2015): #THEN                        !2010-2015
            F2=(float(IY)+float(IDAY)/365.-2010.)/5.
            F1=1.0-F2
            G10=29496.57*F1+29442.0*F2
            G11=-1586.42*F1-1501.0*F2
            H11= 4944.26*F1+4797.1*F2
        else:                                            #!2015-2020
            DT=float(IY)+float(IDAY)/365.-2015.
            G10=29442.0-10.3*DT
            G11=-1501.0+18.1*DT
            H11= 4797.1-26.6*DT
        #ENDIF

#  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EzMAG IN GEO COORD
#  SYSTEM:
#  SIN(TETA0)*COS(LAMBDA0), SIN(TETA0)*SIN(LAMBDA0), AND COS(TETA0)
#         ST0 * CL0                ST0 * SL0                CT0

        SQ=G11**2+H11**2
        SQQ=np.sqrt(SQ)
        SQR=np.sqrt(G10**2+SQ)
        SL0=-H11/SQQ
        CL0=-G11/SQQ
        ST0=SQQ/SQR
        CT0=G10/SQR
        STCL=ST0*CL0
        STSL=ST0*SL0
        CTSL=CT0*SL0
        CTCL=CT0*CL0

#  THE CALCULATIONS ARE TERMINATED IF ONLY GEO-MAG TRANSFORMATION
#  IS TO BE DONE  (IHOUR>24 IS THE AGREED CONDITION FOR THIS CASE):

#   5   IF (IHOUR.GT.24) RETURN
    if (IHOUR > 24): return
    
    SUN(IY,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)

#  S1,S2, AND S3 ARE THE COMPONENTS OF THE UNIT VECTOR EXGSM=EXGSE
#  IN THE SYSTEM GEI POINTING FROM THE EARTH'S CENTER TO THE SUN:

    S1=np.cos(SRASN)*np.cos(SDEC)
    S2=np.sin(SRASN)*np.cos(SDEC)
    S3=np.sin(SDEC)
    CGST=np.cos(GST)
    SGST=np.sin(GST)

#  DIP1, DIP2, AND DIP3 ARE THE COMPONENTS OF THE UNIT VECTOR
#  EZSM=EZMAG IN THE SYSTEM GEI:

    DIP1=STCL*CGST-STSL*SGST
    DIP2=STCL*SGST+STSL*CGST
    DIP3=CT0

#  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EYGSM IN THE SYSTEM
#  GEI BY TAKING THE VECTOR PRODUCT D x S AND NORMALIZING IT TO UNIT
#  LENGTH:

    Y1=DIP2*S3-DIP3*S2
    Y2=DIP3*S1-DIP1*S3
    Y3=DIP1*S2-DIP2*S1
    Y=np.sqrt(Y1*Y1+Y2*Y2+Y3*Y3)
    Y1=Y1/Y
    Y2=Y2/Y
    Y3=Y3/Y

#  THEN IN THE GEI SYSTEM THE UNIT VECTOR Z=EZGSM=EXGSM x EYGSM=S x Y
#  HAS THE COMPONENTS:

    Z1=S2*Y3-S3*Y2
    Z2=S3*Y1-S1*Y3
    Z3=S1*Y2-S2*Y1

#  THE VECTOR EZGSE (HERE DZ) IN GEI HAS THE COMPONENTS (0,-SIN(DELTA),
#  COS(DELTA)) = (0.,-0.397823,0.917462); HERE DELTA = 23.44214 DEG FOR
#  THE EPOCH 1978 (SEE THE BOOK BY GUREVICH OR OTHER ASTRONOMICAL
#  HANDBOOKS). HERE THE MOST ACCURATE TIME-DEPENDENT FORMULA IS USED:

    DJ=float(365*(IY-1900)+(IY-1901)/4 +IDAY)-0.5+float(ISEC)/86400.
    T=DJ/36525.
    OBLIQ=(23.45229-0.0130125*T)/57.2957795
    DZ1=0.
    DZ2=-np.sin(OBLIQ)
    DZ3=np.cos(OBLIQ)

#  THEN THE UNIT VECTOR EYGSE IN GEI SYSTEM IS THE VECTOR PRODUCT DZ x S

    DY1=DZ2*S3-DZ3*S2
    DY2=DZ3*S1-DZ1*S3
    DY3=DZ1*S2-DZ2*S1

#  THE ELEMENTS OF THE MATRIX GSE TO GSM ARE THE SCALAR PRODUCTS:
#  CHI=EM22=(EYGSM,EYGSE), SHI=EM23=(EYGSM,EZGSE),
#  EM32=(EZGSM,EYGSE)=-EM23, AND EM33=(EZGSM,EZGSE)=EM22

    CHI=Y1*DY1+Y2*DY2+Y3*DY3
    SHI=Y1*DZ1+Y2*DZ2+Y3*DZ3
    DECARG=SHI
    if (np.abs(DECARG) > 1.): DECARG = 1.0 * DECARG/np.abs(DECARG) #SIGN(1.,DECARG)
    HI=math.asin(DECARG)

#  TILT ANGLE: PSI=ARCSIN(DIP,EXGSM)

    SPS=DIP1*S1+DIP2*S2+DIP3*S3
    CPS=np.sqrt(1.-SPS**2)
    DECARG=SPS
    if (np.abs(DECARG) > 1.): DECARG=1.0 * DECARG/np.abs(DECARG)   #SIGN(1.,DECARG)
    PSI=math.asin(DECARG)

#  THE ELEMENTS OF THE MATRIX MAG TO SM ARE THE SCALAR PRODUCTS:
#  CFI=GM22=(EYSM,EYMAG), SFI=GM23=(EYSM,EXMAG); THEY CAN BE DERIVED
#  AS FOLLOWS:

#  IN GEO THE VECTORS EXMAG AND EYMAG HAVE THE COMPONENTS
#  (CT0*CL0,CT0*SL0,-ST0) AND (-SL0,CL0,0), RESPECTIVELY. HENCE, IN
#  GEI SYSTEM THE COMPONENTS ARE:
#  EXMAG:    CT0*CL0*COS(GST)-CT0*SL0*SIN(GST)
#            CT0*CL0*SIN(GST)+CT0*SL0*COS(GST)
#            -ST0
#  EYMAG:    -SL0*COS(GST)-CL0*SIN(GST)
#            -SL0*SIN(GST)+CL0*COS(GST)
#             0
#  THE COMPONENTS OF EYSM IN GEI WERE FOUND ABOVE AS Y1, Y2, AND Y3;
#  NOW WE ONLY HAVE TO COMBINE THE QUANTITIES INTO SCALAR PRODUCTS:

    EXMAGX=CT0*(CL0*CGST-SL0*SGST)
    EXMAGY=CT0*(CL0*SGST+SL0*CGST)
    EXMAGZ=-ST0
    EYMAGX=-(SL0*CGST+CL0*SGST)
    EYMAGY=-(SL0*SGST-CL0*CGST)
    CFI=Y1*EYMAGX+Y2*EYMAGY
    SFI=Y1*EXMAGX+Y2*EXMAGY+Y3*EXMAGZ

    XMUT=(math.atan2(SFI,CFI)+3.1415926536)*3.8197186342

#  THE ELEMENTS OF THE MATRIX GEO TO GSM ARE THE SCALAR PRODUCTS:

#  A11=(EXGEO,EXGSM), A12=(EYGEO,EXGSM), A13=(EZGEO,EXGSM),
#  A21=(EXGEO,EYGSM), A22=(EYGEO,EYGSM), A23=(EZGEO,EYGSM),
#  A31=(EXGEO,EZGSM), A32=(EYGEO,EZGSM), A33=(EZGEO,EZGSM),

#  ALL THE UNIT VECTORS IN BRACKETS ARE ALREADY DEFINED IN GEI:

#  EXGEO=(CGST,SGST,0), EYGEO=(-SGST,CGST,0), EZGEO=(0,0,1)
#  EXGSM=(S1,S2,S3),  EYGSM=(Y1,Y2,Y3),   EZGSM=(Z1,Z2,Z3)
#  AND  THEREFORE:

    A11=S1*CGST+S2*SGST
    A12=-S1*SGST+S2*CGST
    A13=S3
    A21=Y1*CGST+Y2*SGST
    A22=-Y1*SGST+Y2*CGST
    A23=Y3
    A31=Z1*CGST+Z2*SGST
    A32=-Z1*SGST+Z2*CGST
    A33=Z3

# 10   FORMAT(/
#     * ' RECALC: GIVEN YEAR',I5,' IS OUT OF INTERVAL 1900-2010'/,
#     * '   *** CALCULATIONS WILL BE DONE FOR YEAR =',I5,' ***'/)

    return
#
#
def SPHCAR(R,TETA,PHI,X,Y,Z,J):
    '''*********************************************************************
        CONVERTS SPHERICAL COORDS INTO CARTESIAN ONES AND VICA VERSA
        (TETA AND PHI IN RADIANS)
                         J>0            J<0
       -----INPUT:   J,R,TETA,PHI     J,X,Y,Z
       ----OUTPUT:      X,Y,Z        R,TETA,PHI
       AUTHOR: NIKOLAI A. TSYGANENKO, INSTITUTE OF PHYSICS, ST.-PETERSBURG
       STATE UNIVERSITY, STARY PETERGOF 198904, ST.-PETERSBURG, RUSSIA
       (now the NASA Goddard Space Fligth Center, Greenbelt, Maryland)
       *********************************************************************'''
    if (J > 0):
        SQ=R*np.sin(TETA)
        X=SQ*np.cos(PHI)
        Y=SQ*np.sin(PHI)
        Z=R*np.cos(TETA)
        return
    SQ=X**2+Y**2
    R=np.sqrt(SQ+Z**2)
    if (SQ != 0.):
        SQ=np.sqrt(SQ)
        PHI=math.atan2(Y,X)
        TETA=math.atan2(SQ,Z)
        if (PHI < 0.): PHI=PHI+6.28318531
        return
    PHI=0.
    if (Z < 0.):
        TETA=3.141592654
        return
    TETA=0.
    return
#      IF(J.GT.0) GOTO 3
#      SQ=X**2+Y**2
#      R=SQRT(SQ+Z**2)
#      IF (SQ.NE.0.) GOTO 2
#      PHI=0.
#      IF (Z.LT.0.) GOTO 1
#      TETA=0.
#      RETURN
#  1   TETA=3.141592654
#      RETURN
#  2   SQ=SQRT(SQ)
#      PHI=ATAN2(Y,X)
#      TETA=ATAN2(SQ,Z)
#      IF (PHI.LT.0.) PHI=PHI+6.28318531
#      RETURN
#  3   SQ=R*SIN(TETA)
#      X=SQ*COS(PHI)
#      Y=SQ*SIN(PHI)
#      Z=R*COS(TETA)

def BSPCAR(TETA,PHI,BR,BTET,BPHI,BX,BY,BZ):
    '''*********************************************************************
        CALCULATES CARTESIAN FIELD COMPONENTS FROM SPHERICAL ONES
        -----INPUT:   TETA,PHI - SPHERICAL ANGLES OF THE POINT IN RADIANS
              BR,BTET,BPHI -  SPHERICAL COMPONENTS OF THE FIELD
        -----OUTPUT:  BX,BY,BZ - CARTESIAN COMPONENTS OF THE FIELD
        AUTHOR: NIKOLAI A. TSYGANENKO, INSTITUTE OF PHYSICS, ST.-PETERSBURG
        STATE UNIVERSITY, STARY PETERGOF 198904, ST.-PETERSBURG, RUSSIA
        (now the NASA Goddard Space Fligth Center, Greenbelt, Maryland)
       *********************************************************************'''
    S=np.sin(TETA)
    C=np.cos(TETA)
    SF=np.sin(PHI)
    CF=np.cos(PHI)
    BE=BR*S+BTET*C
    BX=BE*CF-BPHI*SF
    BY=BE*SF+BPHI*CF
    BZ=BR*C-BTET*S
    return
#
#
def GEOMAG(XGEO,YGEO,ZGEO,XMAG,YMAG,ZMAG,J,IYR):
    '''*********************************************************************
        CONVERTS GEOCENTRIC (GEO) TO DIPOLE (MAG) COORDINATES OR VICA VERSA.
        IYR IS YEAR NUMBER (FOUR DIGITS).
                           J>0                J<0
        -----INPUT:  J,XGEO,YGEO,ZGEO,IYR   J,XMAG,YMAG,ZMAG,IYR
        -----OUTPUT:    XMAG,YMAG,ZMAG        XGEO,YGEO,ZGEO
        AUTHOR: NIKOLAI A. TSYGANENKO, INSTITUTE OF PHYSICS, ST.-PETERSBURG
        STATE UNIVERSITY, STARY PETERGOF 198904, ST.-PETERSBURG, RUSSIA
        (now the NASA Goddard Space Fligth Center, Greenbelt, Maryland)
       *********************************************************************'''

    global ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,AB,K,IY,BB
    II = 1
    if (IYR != II):
        II = IYR
        RECALC(II,0,25,0,0)
        if (J < 0):
            XGEO=XMAG*CTCL-YMAG*SL0+ZMAG*STCL
            YGEO=XMAG*CTSL+YMAG*CL0+ZMAG*STSL
            ZGEO=ZMAG*CT0-XMAG*ST0
            return
        XMAG=XGEO*CTCL+YGEO*CTSL-ZGEO*ST0
        YMAG=YGEO*CL0-XGEO*SL0
        ZMAG=XGEO*STCL+YGEO*STSL+ZGEO*CT0
        return
    return        
#      IF(IYR.EQ.II) GOTO 1
#      II=IYR
#      CALL RECALC(II,0,25,0,0)
#  1   CONTINUE
#      IF(J.LT.0) GOTO 2
#      XMAG=XGEO*CTCL+YGEO*CTSL-ZGEO*ST0
#      YMAG=YGEO*CL0-XGEO*SL0
#      ZMAG=XGEO*STCL+YGEO*STSL+ZGEO*CT0
#      RETURN
#  2   XGEO=XMAG*CTCL-YMAG*SL0+ZMAG*STCL
#      YGEO=XMAG*CTSL+YMAG*CL0+ZMAG*STSL
#      ZGEO=ZMAG*CT0-XMAG*ST0
#      RETURN
#      END
#
#
def MAGSM(XMAG,YMAG,ZMAG,XSM,YSM,ZSM,J):
    '''*********************************************************************
        CONVERTS DIPOLE (MAG) TO SOLAR MAGNETIC (SM) COORDINATES OR VICA VERSA
                        J>0                J<0
       -----INPUT: J,XMAG,YMAG,ZMAG     J,XSM,YSM,ZSM
       ----OUTPUT:    XSM,YSM,ZSM       XMAG,YMAG,ZMAG
       ATTENTION: SUBROUTINE RECALC MUST BE CALLED BEFORE MAGSM IN TWO CASES
       /A/  BEFORE THE FIRST USE OF MAGSM
       /B/  IF THE CURRENT VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC ARE
            DIFFERENT FROM THOSE IN THE PRECEDING CALL OF  MAGSM
        AUTHOR: NIKOLAI A. TSYGANENKO, INSTITUTE OF PHYSICS, ST.-PETERSBURG
        STATE UNIVERSITY, STARY PETERGOF 198904, ST.-PETERSBURG, RUSSIA
        (now the NASA Goddard Space Fligth Center, Greenbelt, Maryland)
       *********************************************************************'''

    global A,SFI,CFI,B,AB,K,IY,BA
    
    if (J < 0):
        XMAG=XSM*CFI+YSM*SFI
        YMAG=YSM*CFI-XSM*SFI
        ZMAG=ZSM
        return
    XSM=XMAG*CFI-YMAG*SFI
    YSM=XMAG*SFI+YMAG*CFI
    ZSM=ZMAG
    return
#      IF (J.LT.0) GOTO 1
#      XSM=XMAG*CFI-YMAG*SFI
#      YSM=XMAG*SFI+YMAG*CFI
#      ZSM=ZMAG
#      RETURN
#  1   XMAG=XSM*CFI+YSM*SFI
#      YMAG=YSM*CFI-XSM*SFI
#      ZMAG=ZSM
#
#      RETURN
#      END
#
#
def SMGSM(XSM,YSM,ZSM,XGSM,YGSM,ZGSM,J):
    '''*********************************************************************
        CONVERTS SOLAR MAGNETIC (SM) TO SOLAR MAGNETOSPHERIC (GSM) COORDINATES
        OR VICA VERSA.
                        J>0                   J<0
        -----INPUT: J,XSM,YSM,ZSM        J,XGSM,YGSM,ZGSM
        ----OUTPUT:  XGSM,YGSM,ZGSM       XSM,YSM,ZSM
        ATTENTION: SUBROUTINE RECALC MUST BE CALLED BEFORE SMGSM IN TWO CASES
        /A/  BEFORE THE FIRST USE OF SMGSM
        /B/  IF THE CURRENT VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC ARE
          DIFFERENT FROM THOSE IN THE PRECEDING CALL OF SMGSM
        AUTHOR: NIKOLAI A. TSYGANENKO, INSTITUTE OF PHYSICS, ST.-PETERSBURG
        STATE UNIVERSITY, STARY PETERGOF 198904, ST.-PETERSBURG, RUSSIA
        (now the NASA Goddard Space Fligth Center, Greenbelt, Maryland)
       *********************************************************************'''

    global A,SPS,CPS,B,K,IY,AB
    if (J < 0):
        XSM=XGSM*CPS-ZGSM*SPS
        YSM=YGSM
        ZSM=XGSM*SPS+ZGSM*CPS
        return
    XGSM=XSM*CPS+ZSM*SPS
    YGSM=YSM
    ZGSM=ZSM*CPS-XSM*SPS
    return
#      IF (J.LT.0) GOTO 1
#      XGSM=XSM*CPS+ZSM*SPS
#      YGSM=YSM
#      ZGSM=ZSM*CPS-XSM*SPS
#      RETURN
#  1   XSM=XGSM*CPS-ZGSM*SPS
#      YSM=YGSM
#      ZSM=XGSM*SPS+ZGSM*CPS
#
#      RETURN
#      END
#
#
def CLCMLT(IYYYY,DDD,UTHR,GLAT,GLON,MLT):
    '''--------------------------------------------------------------------
        calculates magnetic local time
        Inputs:
             IYYYY..Year as YYYY, e.g. 1998
             DDD..day of year (1.1. = 0)
             UTHR..universal time in decimal hours
             GLAT,GLON..latitude north and longitude east in degrees
        Output:
             MLT..magnetic local time in decimal hours
        Required subroutines: DPMTRX
       --------------------------------------------------------------------'''

    global DTOR,PI

    XXM = np.zeros(3)
    YYM = np.zeros(3)
    ZZM = np.zeros(3)
    SA = np.zeros(3)
    SG = np.zeros(3)
    
    XG=np.cos(GLAT*DTOR)*np.cos(GLON*DTOR)
    YG=np.cos(GLAT*DTOR)*np.sin(GLON*DTOR)
    ZG=np.sin(GLAT*DTOR)
    DPMTRX(IYYYY,DDD,XXM,YYM,ZZM)
       
#       transform
    XM=XXM[1]*XG+XXM[2]*YG+XXM[3]*ZG
    YM=YYM[1]*XG+YYM[2]*YG+YYM[3]*ZG
    ZM=ZZM[1]*XG+ZZM[2]*YG+ZZM[3]*ZG
#       
    IHOUR=int(UTHR)
    MIN=int((UTHR-IHOUR)*60)
    ISEC=int((UTHR-IHOUR-MIN/60.0)*3600)
    SUN (IYYYY,DDD+1,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
    BE=GST
    CAL=np.cos(SRASN)
    SA[3]=np.sin(SDEC)
    SA[1]=np.cos(SDEC)
    SA[2]=SA[1]*np.sin(SRASN)
    SA[1]=SA[1]*CAL
    S=np.sin(BE)
    C=np.cos(BE)
    SG[1]=C*SA[1]+S*SA[2]
    SG[2]=C*SA[2]-S*SA[1]
    SG[3]=SA[3]       
#       transform
    SM[1]=XXM[1]*SG[1]+XXM[2]*SG[2]+XXM[3]*SG[3]
    SM[2]=YYM[1]*SG[1]+YYM[2]*SG[2]+YYM[3]*SG[3]
    SM[3]=ZZM[1]*SG[1]+ZZM[2]*SG[2]+ZZM[3]*SG[3]
#      
    LAM=math.atan2(YM,XM)
    LAMS=math.atan2(SM[2],SM[1])
    DELLAM=LAM-LAMS
    if (DELLAM < 0.): DELLAM=DELLAM+2*PI
    MLT=np.mod(DELLAM/PI*12.+12.,24.)
    return
#
#
def DPMTRX(IYYYY,DDD,XM,YM,ZM):
    '''------------------------------------------------------------------------
        calculates othonormal matrix (columns XM,YM,ZM) for transformation
        from geographic to magnetic coordinates
        Inputs:
                IYYYY..year
                DDD..day of year (1.1 = 0)
        Outputs:
                XM,YM,ZM..colums of the matrix
        Notes:
            MX(N),MY(N),MZ(N)..coordinates of the B vector in geographic system
                               for years stored in YR(N)
                               N..number of elements of arrays MX,MY,MZ and YR
       ------------------------------------------------------------------------'''
    
    XM = np.zeros(3)
    YM = np.zeros(3)
    ZM = np.zeros(3)
    YR = np.zeros(10)
    MX = np.zeros(10)
    MY = np.zeros(10)
    MZ = np.zeros(10)
    
    global GHI1,GHI2,GHI3
    
    N = 10
# IGRF coefficients (dipole) calculated in FELDCOF in IGRF.FOR
    MXI = -GHI2
    MYI = -GHI3
    MZI = -GHI1

# normalization of the vector of the dipole exis of the magnetic field
    M=np.sqrt(MXI*MXI+MYI*MYI+MZI*MZI)
    MYZ=np.sqrt(MYI*MYI+MZI*MZI)
    ZM[1]=MXI/M
    ZM[2]=MYI/M
    ZM[3]=MZI/M
    ZM12=np.sqrt(ZM[1]*ZM[1]+ZM[2]*ZM[2])
    YM[1]=-ZM[2]/ZM12
    YM[2]=ZM[1]/ZM12
    YM[3]=0.
    XM[1]=YM[2]*ZM[3]-YM[3]*ZM[2]
    XM[2]=YM[3]*ZM[1]-YM[1]*ZM[3]
    XM[3]=YM[1]*ZM[2]-YM[2]*ZM[1]
    return
#
# --------------------- end IGRF.FOR ----------------------------------
