C====================================================================
C
C UMAT FOR ANISOTROPIC ELASTIC VISCOPLASTIC DAMAGE CONSTITUTIVE LAW
C
C     BY JAKOB SCHWIEDRZIK - ISTB - University of Bern - 2012
C
C     ADAPTED BY DENIS SCHENK - ISTB - UNIVERSITY OF BERN - 2018
C
C     ADAPTED BY DENIS SCHENK - ARTORG - UNIVERSITY OF BERN - 2021
C
C====================================================================
C
C     Modifications 
C   
C     CHANGED TO QUADRIC YIELD SURFACE - JJS, ISTB, 07/2012
C     
C     ADDED PRIMAL CPPA (PEREZ-FOGUET 2002) BASED ON VARIATIONAL 
C     STRUCTURE AS FAIL SAFE FOR NEWTON CPPA TO ENSURE GLOBAL 
C     CONVERGENCE - JJS, ISTB, 08/2012
C
C     ADDED NEW MATERIAL CONSTANTS: FABRIC-BASED ORTHOTROPIC
C     FOR TRABECULAR BONE (PHZ 2013) AND ISOTROPIC FOR CORTICAL
C     BONE (JJS 2013) - JJS, ISTB, 02/2013                           
C
C     ADAPTED THE UMAT FOR MIXTURE BETWEEN TRABECULAR AND
C     CORTICAL BONE VIA PARTIAL BONE VOLUME (PBV)- GM, ISTB, 06/2017
C
C     IMPLEMENTED NEW MATERIAL SUPERPOSITION OF SSSS, FFFF, FF AND
C     POSTYIELD BEHAVIOR FLAG 5 (SIMPLE SOFTENING FOR USE IN
C     MATERIAL SUPERPOSITION - DS, ISTB, 11/2018
C
C====================================================================
C
C     Variables provided by Abaqus:
C
C     STRAN(NTENS)  - strains, accounting for rigid rotation
C     DSTRAN(NTENS) - strain increments
C     TIME(1)       - step time at beginning of increment
C     TIME(2)       - total time at beginning of increment
C     DTIME         - time increment 
C     TEMP          - temperature at beginning
C     DTEMP         - increment of temperature
C     PREDEF        - array of predefined field variables
C     DPRED         - increment in predefind field variables
C     CMNAME        - user-defined material name
C     NDI           - number of direct shear components
C     NSHR          - number of engineering shear stress comps.
C     NTENS         - size of the stress|strain array (NDI+NSHR)
C     NSTATV        - number of stat variables 
C     PROPS(NPROPS) - user-specified material constants
C     NPROPS        - user-defined number of material conts.
C     COORDS(3)     - array of coordinates of point (current for NL)
C     DROT(3,3)     - rotation increment matrix (for state vars)
C     CELENT        - char. element length
C     DFGRD0(3,3)   - deformation gradient at beginning of incr.
C     DFGRD1(3,3)   - deformation gradient at end of increment
C     NOEL          - element number
C     NPT           - integration point number
C     LAYER         - layer of composite shell|solid
C     KSPT          - section point within layer
C     KSTEP         - step number
C     KINC          - increment number
C
C     PERSONALIZED LOADING
C     EPSO          - isotropic homogeneous hydrostatic compresive reference strain
C     R             - BVTV exponent for optimization
C     OF            - value of objective function in each element
C     
C
C====================================================================
C
C     Variables to be defined:
C
C     STRESS(NTENS)       - Cauchy stress at end of increment
C     DDSDDE(NTENS,NTENS) - Jacobian Stiffness matrix 
C     STATEV(NSTATV)      - Solution dependent state variables
C
C====================================================================
C
C     State variables:  
C
C     IF NUMBER IS CHANGED, CHANGE AS WELL IN LINE 1417!
C     SDV 1:     CUMULATED PLASTIC STRAIN
C     SDV 2:     SCALAR DAMAGE VARIABLE
C     SDV 3-8:   PLASTIC STRAIN VECTOR
C     SDV 9-14:  NOMINAL STRESS VECTOR
C     SDV 15:    BVTVC
C     SDV 16:    BVTVT
C     SDV 22:    OF VALUE
C
C====================================================================
C
C     User-defined properties:
C
C     PROPS(1)     - BVTVC (Greyvalue-based BVTVd of the cortical phase of the element)
C     PROPS(2)     - BVTVT (Greyvalue-based BVTVd of the trabecular phase of the element)
C     PROPS(3)     - PBVC (Volume of the cortical phase/Volume of the element)
C     PROPS(4)     - PBVT (Volume of the trabecular phase/Volume of the element)
C     PROPS(5)     - M1
C     PROPS(6)     - M2
C     PROPS(7)     - M3
C     PROPS(8-16)  - Stiffness matrix
C     PROPS(17-25) - FFFF tensor
C     PROPS(26-28) - FF tensor
C
C====================================================================
C
C     Hard wired flags:
C
C     VISCOSITY: 0 - RATE-INDEPENDENT (DEFAULT)
C                1 - LINEAR VISCOSITY
C                2 - EXPONENTIAL VISC
C                3 - LOGARITHMIC VISC
C                4 - POLYNOMIAL VISC
C                5 - POWER LAW VISC 
C     POSTYIELD: 0 - PERFECT PLASTICITY
C                1 - EXPONENTIAL HARDENING
C                2 - SIMPLE SOFTENING
C                3 - EXPONENTIAL SOFTENING
C                4 - PIECEWISE SOFTENING
C                5 - SIMPLE SOFTENING CORRECTED (DENIS)
C
C====================================================================
C
C     Input file:  
C
C     *MATERIAL, NAME=CardName
C     *DEPVAR
C     18
C     *USER MATERIAL, CONSTANTS=7, UNSYMM, TYPE=MECHANICAL
C     BVTVC, BVTVT, PBVC, PBVT, M3, M2, M1
C
C====================================================================

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,
     & DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,
     & PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,
     & COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,
     & KSPT,KSTEP,KINC)
C
      IMPLICIT NONE
C
C     INPUT/OUTPUT VARIABLES:
      INTEGER NTENS,NDI,NSHR
      INTEGER NSTATV
      INTEGER NPROPS
C
C     INPUT/OUTPUT VARIABLES THAT HAVE TO BE DEFINED BY THE UMAT
C     IN ALL SITUATIONS:
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV)
      DOUBLE PRECISION DDSDDE(NTENS,NTENS)
      DOUBLE PRECISION SSE,SPD,SCD
C
C     ONLY IN FULLY COUPLED TEMPERATURE-DISPLACEMENT ANALYSIS:
      DOUBLE PRECISION RPL,DDSDDT(NTENS),DRPLDE(NTENS),DRPLDT
C
C     OUTPUT VARIABLES THAT CAN BE UPDATED WITHIN THE UMAT (TIME INC):
      DOUBLE PRECISION PNEWDT
C
C     INPUT VARIABLES: SUPPLIED BY ABAQUS:
      DOUBLE PRECISION STRAN(NTENS),DSTRAN(NTENS)
      DOUBLE PRECISION TIME(2),DTIME
      DOUBLE PRECISION TEMP,DTEMP
      DOUBLE PRECISION PREDEF(1),DPRED(1)
      CHARACTER*8      CMNAME
      DOUBLE PRECISION PROPS(NPROPS)
      DOUBLE PRECISION COORDS(3)
      DOUBLE PRECISION DROT(3,3)
      DOUBLE PRECISION CELENT
      DOUBLE PRECISION DFGRD0(3,3),DFGRD1(3,3)
      INTEGER          NOEL,NPT,LAYER,KSPT,KSTEP,KINC
C
C     EXTERNAL FUNCTIONS
      DOUBLE PRECISION RADK,DRADK,DAM,DDAM,VECNORM 
      DOUBLE PRECISION VISCY,VDYDK,TSFU
C
C     VISCOPLASTIC DAMAGE UMAT VARIABLES
      INTEGER ITER,MAXITER,K1,K2,DENSFL,VISCFL,PYFL
      DOUBLE PRECISION MM,ETA,TOL,KONSTD,CRITD
      DOUBLE PRECISION BVTV,BVTVC,BVTVT,PBVC,PBVT
      DOUBLE PRECISION MM1,MM2,MM3,SCA1,SCA2, SCA3
      DOUBLE PRECISION E0,V0,MU0,EAA,VA0,MUA0,KS,LS
      DOUBLE PRECISION SIGD0P,SIGD0N,ZETA0,TAUD0,PP,QQ
      DOUBLE PRECISION SIGDAP,SIGDAN,ZETAA0,TAUDA0,S0,SA
      DOUBLE PRECISION KSLOPE,KMAX,KWIDTH,RDY
      DOUBLE PRECISION KMIN,GMIN,ND,EXPS
      DOUBLE PRECISION SSSS(NTENS,NTENS)
      DOUBLE PRECISION SSSST(NTENS,NTENS)
      DOUBLE PRECISION SSSSC(NTENS,NTENS)
      DOUBLE PRECISION CCCC(NTENS,NTENS)
      DOUBLE PRECISION CCCCT(NTENS,NTENS)
      DOUBLE PRECISION CCCCC(NTENS,NTENS)
      DOUBLE PRECISION FFFFT(NTENS,NTENS)
      DOUBLE PRECISION FFFFC(NTENS,NTENS)
      DOUBLE PRECISION FFFF(NTENS,NTENS)
      DOUBLE PRECISION FFFFTI(NTENS,NTENS)
      DOUBLE PRECISION FFFFCI(NTENS,NTENS)
      DOUBLE PRECISION FFFFTCI(NTENS,NTENS)
      DOUBLE PRECISION FFT(3,3)
      DOUBLE PRECISION FFTI(3,3)
      DOUBLE PRECISION FFC(3,3)
      DOUBLE PRECISION FFCI(3,3)
      DOUBLE PRECISION FFTCI(3,3)
      DOUBLE PRECISION FF(NTENS)
      DOUBLE PRECISION FFM(3, 3)
      DOUBLE PRECISION DDSY(NTENS,NTENS),DSY(NTENS)
      DOUBLE PRECISION ETOT1(NTENS),ETOT0(NTENS),SS0(NTENS)
      DOUBLE PRECISION KAPPA0,DKAPPAI,DKAPPA1,DDK1,KAPPA1
      DOUBLE PRECISION STR(NTENS),SSI(NTENS),SS1(NTENS),DSS1(NTENS)
      DOUBLE PRECISION EPLAS0(NTENS),EPLAS1(NTENS)
      DOUBLE PRECISION YSTR,YS,DYDS(NTENS),DYDK,DMG0
      DOUBLE PRECISION DMG,DDMG,RAD,DRAD,FFS(NTENS),SFFS,SRA(NTENS)
      DOUBLE PRECISION NP(NTENS),DNPDS(NTENS,NTENS),DNPDK(NTENS)
      DOUBLE PRECISION RR(NTENS),DRRDS(NTENS,NTENS),DRRDK(NTENS)
      DOUBLE PRECISION HI,DHDS(NTENS),DHDK,NORMRR,ABSY
      DOUBLE PRECISION SSSA(NTENS,NTENS),TANM(NTENS,NTENS)
C
C     PRIMAL CPP ALGORITHM AND LINE SEARCH VARIABLES
      INTEGER ITERL,FLAG
      DOUBLE PRECISION J11(6,6),J22,J12(6),J21(6)
      DOUBLE PRECISION RRR(7),RRI(7),XXI(7),XX1(7),DDD(7),JJJ(7,7)
      DOUBLE PRECISION INVJ(7,7),DDJ(7,7),AUX,DMK,UPLIM
      DOUBLE PRECISION ETAL,BETA,MKI,MK1,ALPHA1,ALPHA2,ALPHA
C     VARIABLES FOR MATERIAL SUPERPOSITION
      INTEGER NROTC
      INTEGER NROTT
      INTEGER I, J
      DOUBLE PRECISION VC(6,6)
      DOUBLE PRECISION DC(6)
      DOUBLE PRECISION VT(6,6)
      DOUBLE PRECISION DT(6)
      DOUBLE PRECISION FFFFINTC(6,6)
      DOUBLE PRECISION FFFFINTT(6,6)
      DOUBLE PRECISION FFFFSUP(6,6)
      DOUBLE PRECISION FFFFSUPMUL(6,6)
C
C     PERSONALIZED LOADINF
      DOUBLE PRECISION EPS0, R, OF
C      
C     INTERFACE FOR ARRAY VALUED EXTERNAL FUNCTION VECDYAD
      INTERFACE
      FUNCTION VECDYAD(BB,CC)
      DOUBLE PRECISION VECDYAD(6,6)
      DOUBLE PRECISION BB(6),CC(6)
      END FUNCTION VECDYAD
C
C     INTERFACE FOR VECTOR VALUED EXTERNAL FUNCTION VDYDS
      FUNCTION VDYDS(VISCFL,ETA,MM,DTIME,DKAPPA,DHDS,HI)
      DOUBLE PRECISION VDYDS(6)
      DOUBLE PRECISION ETA,MM,DTIME,DKAPPA,DHDS(6),HI
      INTEGER VISCFL
      END FUNCTION VDYDS
C
C     INTERFACE FOR VECTOR DOT PRODUCT
      FUNCTION DOTP(DD,EE)
      DOUBLE PRECISION DOTP
      DOUBLE PRECISION DD(6)
      DOUBLE PRECISION EE(6)
      END FUNCTION DOTP
      END INTERFACE
C     _______________________________________________________________
C
C     VARIABLE INITIALISATION
C     _______________________________________________________________
C     TENSOR INITIALISATION
      SSSST  = 0.0D0
      SSSSC  = 0.0D0
      SSSS   = 0.0D0
      CCCCT  = 0.0D0
      CCCCC  = 0.0D0
      CCCC   = 0.0D0
      FFFFT  = 0.0D0
      FFFFTI = 0.0D0
      FFFFC  = 0.0D0
      FFFFCI = 0.0D0
      FFFFTCI= 0.0D0
      FFFF   = 0.0D0
      DDSY   = 0.0D0
      DNPDS  = 0.0D0
      DRRDS  = 0.0D0
      SSSA   = 0.0D0
      TANM   = 0.0D0
      JJJ    = 0.0D0
      INVJ   = 0.0D0
      DDJ    = 0.0D0
      J11    = 0.0D0
      VC     = 0.0D0
      DC     = 0.0D0
      VT     = 0.0D0
      DT     = 0.0D0
      FFFFINTC= 0.0D0
      FFFFINTT= 0.0D0
      FFFFSUP= 0.0D0
C
C     VECTOR INITIALISATION
      FFT    = 0.0D0
      FFTI   = 0.0D0
      FFC    = 0.0D0
      FFCI   = 0.0D0
      FFTCI  = 0.0D0
      FF     = 0.0D0
      FFM    = 0.0D0
      ETOT1  = 0.0D0
      ETOT0  = 0.0D0
      SS0    = 0.0D0
      STR    = 0.0D0
      SSI    = 0.0D0
      SS1    = 0.0D0
      DSS1   = 0.0D0
      EPLAS0 = 0.0D0
      EPLAS1 = 0.0D0
      DYDS   = 0.0D0
      FFS    = 0.0D0
      NP     = 0.0D0
      DNPDK  = 0.0D0
      RR     = 0.0D0
      DRRDK  = 0.0D0
      DHDS   = 0.0D0
      RRR    = 0.0D0
      RRI    = 0.0D0
      XXI    = 0.0D0
      XX1    = 0.0D0
      DDD    = 0.0D0
      J12    = 0.0D0
      J21    = 0.0D0
C
C     SCALAR INITIALISATION
      ITER    = 0
      ITERL   = 0
      FLAG    = 1
      K1      = 0
      K2      = 0
      MM      = 0.0D0
      ETA     = 0.0D0
      RDY     = 0.0D0
      KAPPA0  = 0.0D0
      DKAPPAI = 0.0D0
      DKAPPA1 = 0.0D0
      DDK1    = 0.0D0
      KAPPA1  = 0.0D0
      YSTR    = 0.0D0
      YS      = 0.0D0
      DYDK    = 0.0D0
      DMG     = 0.0D0
      DMG0    = 0.0D0
      DDMG    = 0.0D0
      RAD     = 0.0D0
      DRAD    = 0.0D0
      SFFS    = 0.0D0
      HI      = 0.0D0
      DHDK    = 0.0D0
      NORMRR  = 0.0D0
      ABSY    = 1.0D0
      KSLOPE  = 1.0D0
      KMAX    = 1.0D0
      KWIDTH  = 1.0D0
      KONSTD  = 0.0D0
      CRITD   = 0.0D0
      S0      = 0.0D0
      ETAL    = 0.0D0
      BETA    = 0.0D0
      MKI     = 0.0D0
      MK1     = 0.0D0
      ALPHA1  = 0.0D0
      ALPHA2  = 0.0D0
      ALPHA   = 0.0D0
      AUX     = 0.0D0
      UPLIM   = 0.0D0
      DMK     = 0.0D0
      J22     = 0.0D0
      BVTV    = 0.0D0
      BVTVC   = 0.0D0
      BVTVT   = 0.0D0
      PBVC    = 0.0D0
      PBVT    = 0.0D0
      MM1     = 0.0D0
      MM2     = 0.0D0
      MM3     = 0.0D0
      SCA1    = 0.0D0
      SCA2    = 0.0D0
      SCA3    = 0.0D0
      I       = 0.0D0
      J       = 0.0D0
      NROTC   = 0.0D0
      NROTT   = 0.0D0
C
C     PERSONALIZED LOADING
      EPS0    = 0.001D0
      R       = 0.0D0
      OF      = 0.0D0
C     _______________________________________________________________
C
C     ANALYSE PARAMETERS
C     _______________________________________________________________
C     TOLERANCE OF NUMERICAL ERROR  
      TOL = 1.0D-12
C
C     MAXIMUM ITERATIONS ALLOWED  
      MAXITER = 100
C
C     INPUTFILE PARAMETERS (Read from .inp file)
      BVTVC = PROPS(1)
      BVTVT = PROPS(2)
      PBVC  = PROPS(3)
      PBVT  = PROPS(4)
      MM1   = PROPS(5)
      MM2   = PROPS(6)
      MM3   = PROPS(7)
C      
C     CORRECT ELEMENTS WITH PBV =! 0 BUT BVTV = 0
      IF(PBVC.GT.0.0D0.AND.BVTVC.LT.0.01D0) THEN
        BVTVC = 0.01D0
      ENDIF
      IF(PBVT.GT.0.0D0.AND.BVTVT.LT.0.01D0) THEN
        BVTVT = 0.01D0
      ENDIF
C    
C     VISCOPLASTICITY PARAMETERS
      ETA = 1.D-4
      MM  = 1.0D0
C
C     LINE SEARCH PARAMETERS
      BETA = 1.D-4
      ETAL = 0.25D0
C
C     DAMAGE PARAMETERS (WOLFRAM J BIOMECH 2011)
      KONSTD = 8.0075D0
      CRITD  = 0.85D0
C
C     YIELD/STRENGTH RATIO
C     Adapted Denis 2018: RDY = 0.7D0+0.29D0*(PHIC/(PHIC+PHIT)) BEFORE: RDY = 0.7D0
C     Only for accurate HFE! PYFL = 5!
      RDY = 0.7D0
C
C     TO BE FITTED
      SCA1 = 1.3D0
      SCA2 = 1.52D0
      KMAX = 0.01D0
      SCA3 = 1.0D0
C    
C     _______________________________________________________________
C
C     MATERIAL PARAMETERS
C     _______________________________________________________________
C     FABRIC- AND DENSITY-BASED ORTHOTROPIC BONE, MAIN DIRECTION 3
C     ___________________________________________________________________
C     FOR TRABECULAR BONE
C     FLAGS (See above YIELD/STRENGTH RATIO)
      VISCFL = 0
      PYFL   = 2
C     POSTYIELD PARAMETERS (KMAX to be scaled)
      KSLOPE = 1000.D0
      KWIDTH = 8.0D0
      KMIN   = 0.1D0
      GMIN   = -2.0D0
      EXPS   = 300.0D0
      ND     = 2.0D0
C_______________________________________________________________
C_____MIXED PHASE ELEMENTS______________________________________
C_______________________________________________________________
C
      IF (PBVT.GT.0.0D0.AND.PBVC.GT.0.0D0) THEN
C       ELASTICITY PARAMETERS (FABRGWTB PMUBC, PHZ, 2013) NEW VERSION 2017
        E0  = 9759.0D0*(12000.0D0/9759.0D0)*SCA1
        V0  = 0.2278D0
        MU0 = 3117.0D0*(12000.0D0/9759.0D0)*SCA1
        KS  = 1.91D0
        LS  = 1.1D0
C       STRENGTH PARAMETERS (FABRGWTB, PHZ, 2013)
        SIGD0P = 57.69D0*SCA2
        SIGD0N = 73.1D0*SCA2
        ZETA0  = 0.28D0
        TAUD0  = 29.61D0*SCA2
        PP     = 1.82D0
        QQ     = 0.98D0
C       STIFFNESS TENSOR SSSS
        SSSST(1,1) = E0*(1.0D0-V0)/(1.0D0+V0)/(1.0D0-2.0D0*V0)*MM1
     &             **(2.0D0*LS)*(BVTVT**KS)
        SSSST(2,2) = E0*(1.0D0-V0)/(1.0D0+V0)/(1.0D0-2.0D0*V0)*MM2
     &             **(2.0D0*LS)*(BVTVT**KS)
        SSSST(3,3) = E0*(1.0D0-V0)/(1.0D0+V0)/(1.0D0-2.0D0*V0)*MM3
     &             **(2.0D0*LS)*(BVTVT**KS)
        SSSST(4,4) = 2.0D0*MU0*(MM1**LS)*(MM2**LS)*(BVTVT**KS)
        SSSST(5,5) = 2.0D0*MU0*(MM3**LS)*(MM1**LS)*(BVTVT**KS)
        SSSST(6,6) = 2.0D0*MU0*(MM2**LS)*(MM3**LS)*(BVTVT**KS)
        SSSST(2,1) = E0*V0/(1.0D0+V0)/(1.0D0-2.0D0*V0)*(MM1**LS)
     &             *(MM2**LS)*(BVTVT**KS)
        SSSST(3,1) = E0*V0/(1.0D0+V0)/(1.0D0-2.0D0*V0)*(MM3**LS)
     &             *(MM1**LS)*(BVTVT**KS)
        SSSST(3,2) = E0*V0/(1.0D0+V0)/(1.0D0-2.0D0*V0)*(MM2**LS)
     &             *(MM3**LS)*(BVTVT**KS)
        SSSST(1,2) = SSSST(2,1)
        SSSST(1,3) = SSSST(3,1)
        SSSST(2,3) = SSSST(3,2)
C
C       COMPLIANCE TENSOR CCCC
        CCCCT(1,1) = 1.0D0/(E0*(MM1**(2.0D0*LS))*(BVTVT**KS))
        CCCCT(2,2) = 1.0D0/(E0*(MM2**(2.0D0*LS))*(BVTVT**KS))
        CCCCT(3,3) = 1.0D0/(E0*(MM3**(2.0D0*LS))*(BVTVT**KS))
        CCCCT(4,4) = 0.5D0/(MU0*(MM1**LS)*(MM2**LS)*(BVTVT**KS))
        CCCCT(5,5) = 0.5D0/(MU0*(MM3**LS)*(MM1**LS)*(BVTVT**KS))
        CCCCT(6,6) = 0.5D0/(MU0*(MM2**LS)*(MM3**LS)*(BVTVT**KS))
        CCCCT(2,1) = -V0/(E0*(MM1**LS)*(MM2**LS)*(BVTVT**KS))
        CCCCT(3,1) = -V0/(E0*(MM1**LS)*(MM3**LS)*(BVTVT**KS))
        CCCCT(3,2) = -V0/(E0*(MM2**LS)*(MM3**LS)*(BVTVT**KS))
        CCCCT(1,2) = CCCCT(2,1)
        CCCCT(1,3) = CCCCT(3,1)
        CCCCT(2,3) = CCCCT(3,2)
C
C       QUADRIC FOURTH ORDER TENSOR FFFF
        S0        = (SIGD0P+SIGD0N)/2.0D0/SIGD0P/SIGD0N
        FFFFT(1,1) = S0**2/((BVTVT**PP)*MM1**(2.0D0*QQ))**2
        FFFFT(2,2) = S0**2/((BVTVT**PP)*MM2**(2.0D0*QQ))**2
        FFFFT(3,3) = S0**2/((BVTVT**PP)*MM3**(2.0D0*QQ))**2
        FFFFT(1,2) = -(ZETA0*(MM1/MM2)**(2.0D0*QQ))
     &             *S0**2/((BVTVT**PP)*MM1**(2.0D0*QQ))**2
        FFFFT(1,3) = -(ZETA0*(MM1/MM3)**(2.0D0*QQ))
     &             *S0**2/((BVTVT**PP)*MM1**(2.0D0*QQ))**2
        FFFFT(2,3) = -(ZETA0*(MM2/MM3)**(2.0D0*QQ))
     &             *S0**2/((BVTVT**PP)*MM2**(2.0D0*QQ))**2
        FFFFT(2,1) = FFFFT(1,2)
        FFFFT(3,1) = FFFFT(1,3)
        FFFFT(3,2) = FFFFT(2,3)
        FFFFT(4,4) = 0.5D0/((TAUD0*(BVTVT**PP)
     &             *(MM1*MM2)**QQ)**2)
        FFFFT(5,5) = 0.5D0/((TAUD0*(BVTVT**PP)
     &             *(MM1*MM3)**QQ)**2)
        FFFFT(6,6) = 0.5D0/((TAUD0*(BVTVT**PP)
     &             *(MM2*MM3)**QQ)**2)
C
C     QUADRIC SECOND ORDER TENSOR FF
        FFT(1,1) = -(SIGD0P-SIGD0N)/2.0D0/SIGD0P/SIGD0N
     &         /((BVTVT**PP)*MM1**(2.0D0*QQ))
        FFT(2,2) = -(SIGD0P-SIGD0N)/2.0D0/SIGD0P/SIGD0N
     &         /((BVTVT**PP)*MM2**(2.0D0*QQ))
        FFT(3,3) = -(SIGD0P-SIGD0N)/2.0D0/SIGD0P/SIGD0N
     &         /((BVTVT**PP)*MM3**(2.0D0*QQ))
C
C       ELASTICITY PARAMETERS (TIRBCT, JJS, 2013) CORTEX
        E0   = 15079.6D0*(16000.0D0/24578.5D0)*SCA1
        EAA  = 24578.5D0*(16000.0D0/24578.5D0)*SCA1
        V0   = 0.4620D0
        VA0  = 0.354D0
        MU0 = 3117.0D0*(12000.0D0/9759.0D0)*SCA1
        MUA0 = 6578.01D0*(16000.0D0/24578.5D0)*SCA1
        KS   = 1.0D0
C
C       STRENGTH PARAMETERS (TIRBCT, JJS, 2013) CORTEX
        SIGD0P = 56.54D0*(16000.0D0/24578.5D0)*SCA2
        SIGD0N = 201.28D0*(16000.0D0/24578.5D0)*SCA2
        SIGDAP = 176.40D0*(16000.0D0/24578.5D0)*SCA2
        SIGDAN = 268.00D0*(16000.0D0/24578.5D0)*SCA2
        ZETA0  = 0.0074D0
        ZETAA0 = 1.4045D0
        TAUDA0 = 82.55D0*(16000.0D0/24578.5D0)*SCA2
        PP     = 1.0D0
C
C
C       TRANSVERSELY ISOTROPIC COMPLIANCE TENSOR CCCC
        MU0       = E0/2.0D0/(1.0D0+V0) 
        CCCCC(1,1) = 1.0D0/E0
        CCCCC(2,2) = 1.0D0/E0
        CCCCC(3,3) = 1.0D0/EAA
        CCCCC(4,4) = 0.5D0/MU0
        CCCCC(5,5) = 0.5D0/MUA0
        CCCCC(6,6) = 0.5D0/MUA0
        CCCCC(2,1) = -V0/E0
        CCCCC(3,1) = -VA0/EAA
        CCCCC(3,2) = -VA0/EAA
        CCCCC(1,2) = CCCCC(2,1)
        CCCCC(1,3) = CCCCC(3,1)
        CCCCC(2,3) = CCCCC(3,2)
        CCCCC      = CCCCC/BVTVC**KS
C
C       TRANSVERSELY ISOTROPIC STIFFNESS TENSOR SSSS
        CALL MIGS(CCCCC,6,SSSSC)
C
C       TRANSVERSELY ISOTROPIC QUADRIC FOURTH ORDER TENSOR FFFF
        S0        = (SIGD0P+SIGD0N)/2.0D0/SIGD0P/SIGD0N
        SA        = (SIGDAP+SIGDAN)/2.0D0/SIGDAP/SIGDAN
        TAUD0     = DSQRT(0.5D0/S0**2/(1.0D0+ZETA0))
        FFFFC(1,1) = S0**2
        FFFFC(2,2) = S0**2
        FFFFC(3,3) = SA**2
        FFFFC(1,2) = -ZETA0*S0**2
        FFFFC(1,3) = -ZETAA0*SA**2
        FFFFC(2,3) = -ZETAA0*SA**2
        FFFFC(2,1) = FFFFC(1,2)
        FFFFC(3,1) = FFFFC(1,3)
        FFFFC(3,2) = FFFFC(2,3)
        FFFFC(4,4) = 0.5D0/TAUD0**2
        FFFFC(5,5) = 0.5D0/TAUDA0**2
        FFFFC(6,6) = 0.5D0/TAUDA0**2
        FFFFC      = FFFFC/BVTVC**PP/BVTVC**PP
C 
C       TRANSVERSELY ISOTROPIC QUADRIC SECOND ORDER TENSOR FF
        FFC(1,1) = -(SIGD0P-SIGD0N)/2.0D0/SIGD0P/SIGD0N
        FFC(2,2) = -(SIGD0P-SIGD0N)/2.0D0/SIGD0P/SIGD0N
        FFC(3,3) = -(SIGDAP-SIGDAN)/2.0D0/SIGDAP/SIGDAN
        FFC    = FFC/BVTVC**PP
C
C     MATERIAL SUPERPOSITION
C
        SSSS = SSSSC*PBVC+SSSST*PBVT
        CALL MIGS(SSSS, 6, CCCC)
C
C	STEP 1) FFFF INT
C	----------------
	CALL JACOBI(FFFFC, 6, DC, VC, NROTC)
	DO I=1, 6
	  FFFFINTC=FFFFINTC+(1.0D0/DSQRT(DC(I))*VECDYAD(VC(1:6,I), VC(1:6,I)))/DOTP(VC(1:6,I), VC(1:6,I))
	END DO
C	
	CALL JACOBI(FFFFT, 6, DT, VT, NROTC)
	DO I=1, 6
	  FFFFINTT=FFFFINTT+((1.0D0/DSQRT(DT(I))*VECDYAD(VT(1:6,I), VT(1:6,I))))/DOTP(VT(1:6,I), VT(1:6,I))
	END DO
C
C	STEP 2) SUPERPOSITION
C	----------------
	FFFFSUP = PBVC*FFFFINTC+PBVT*FFFFINTT
	FFFFSUPMUL = MATMUL(FFFFSUP,FFFFSUP)
	CALL MIGS(FFFFSUPMUL,6,FFFF)
        CALL MIGS(FFT,3,FFTI)
        CALL MIGS(FFC,3,FFCI)
        FFTCI = FFCI*PBVC+FFTI*PBVT
        CALL MIGS(FFTCI, 3, FFM)
        FF(1) = FFM(1,1)
        FF(2) = FFM(2,2)
        FF(3) = FFM(3,3)
        FF(4) = 0.0D0
        FF(5) = 0.0D0
        FF(6) = 0.0D0
C     IF ONLY ONE PHASE IS PRESENT
C     ----------------------------------------------------------------------
C     TRABECULAR SSSS, FFFF, FF
C     POSTYIELD PARAMETERS TRAB: KS=1000/KW=8/KMIN=0.1/GMIN=-2/EXPS=300/ND=2
C     KSLOPE CHANGED Denis 22.10.18 (300)
C
C_______________________________________________________________
C_____TRABECULAR ELEMENTS_______________________________________
C_______________________________________________________________
C
      ELSE IF (PBVT.GT.0.0D0.AND.PBVC.EQ.0.0D0) THEN
C      
C       BVTV was BVTVC*PBVC before. But as we now superimpose the stiffness matrices, I think we should change this to BVTV=BVTVC!
        BVTV = BVTVT
C        
C       ELASTICITY PARAMETERS (FABRGWTB PMUBC, PHZ, 2013) NEW VERSION 2017 TRABECULAR BONE
        E0  = 9759.0D0*(12000.0D0/9759.0D0)*SCA1
        V0  = 0.2278D0
        MU0 = 3117.0D0*(12000.0D0/9759.0D0)*SCA1
        KS  = 1.91D0
        LS  = 1.1D0
C        
C       STRENGTH PARAMETERS (FABRGWTB, PHZ, 2013) TRABECULAR BONE
        SIGD0P = 57.69D0*SCA2
        SIGD0N = 73.1D0*SCA2
        ZETA0  = 0.28D0
        TAUD0  = 29.61D0*SCA2
        PP     = 1.82D0
        QQ     = 0.98D0
C        
C       STIFFNESS TENSOR SSSS
        SSSST(1,1) = E0*(1.0D0-V0)/(1.0D0+V0)/(1.0D0-2.0D0*V0)*MM1
     &             **(2.0D0*LS)*(BVTV**KS)
        SSSST(2,2) = E0*(1.0D0-V0)/(1.0D0+V0)/(1.0D0-2.0D0*V0)*MM2
     &             **(2.0D0*LS)*(BVTV**KS)
        SSSST(3,3) = E0*(1.0D0-V0)/(1.0D0+V0)/(1.0D0-2.0D0*V0)*MM3
     &             **(2.0D0*LS)*(BVTV**KS)
        SSSST(4,4) = 2.0D0*MU0*(MM1**LS)*(MM2**LS)*(BVTV**KS)
        SSSST(5,5) = 2.0D0*MU0*(MM3**LS)*(MM1**LS)*(BVTV**KS)
        SSSST(6,6) = 2.0D0*MU0*(MM2**LS)*(MM3**LS)*(BVTV**KS)
        SSSST(2,1) = E0*V0/(1.0D0+V0)/(1.0D0-2.0D0*V0)*(MM1**LS)
     &             *(MM2**LS)*(BVTV**KS)
        SSSST(3,1) = E0*V0/(1.0D0+V0)/(1.0D0-2.0D0*V0)*(MM3**LS)
     &             *(MM1**LS)*(BVTV**KS)
        SSSST(3,2) = E0*V0/(1.0D0+V0)/(1.0D0-2.0D0*V0)*(MM2**LS)
     &             *(MM3**LS)*(BVTV**KS)
        SSSST(1,2) = SSSST(2,1)
        SSSST(1,3) = SSSST(3,1)
        SSSST(2,3) = SSSST(3,2)
C
C       COMPLIANCE TENSOR CCCC
        CCCCT(1,1) = 1.0D0/(E0*(MM1**(2.0D0*LS))*(BVTV**KS))
        CCCCT(2,2) = 1.0D0/(E0*(MM2**(2.0D0*LS))*(BVTV**KS))
        CCCCT(3,3) = 1.0D0/(E0*(MM3**(2.0D0*LS))*(BVTV**KS))
        CCCCT(4,4) = 0.5D0/(MU0*(MM1**LS)*(MM2**LS)*(BVTV**KS))
        CCCCT(5,5) = 0.5D0/(MU0*(MM3**LS)*(MM1**LS)*(BVTV**KS))
        CCCCT(6,6) = 0.5D0/(MU0*(MM2**LS)*(MM3**LS)*(BVTV**KS))
        CCCCT(2,1) = -V0/(E0*(MM1**LS)*(MM2**LS)*(BVTV**KS))
        CCCCT(3,1) = -V0/(E0*(MM1**LS)*(MM3**LS)*(BVTV**KS))
        CCCCT(3,2) = -V0/(E0*(MM2**LS)*(MM3**LS)*(BVTV**KS))
        CCCCT(1,2) = CCCCT(2,1)
        CCCCT(1,3) = CCCCT(3,1)
        CCCCT(2,3) = CCCCT(3,2)
C
C       QUADRIC FOURTH ORDER TENSOR FFFF
        S0        = (SIGD0P+SIGD0N)/2.0D0/SIGD0P/SIGD0N
        FFFFT(1,1) = S0**2/((BVTV**PP)*MM1**(2.0D0*QQ))**2
        FFFFT(2,2) = S0**2/((BVTV**PP)*MM2**(2.0D0*QQ))**2
        FFFFT(3,3) = S0**2/((BVTV**PP)*MM3**(2.0D0*QQ))**2
        FFFFT(1,2) = -(ZETA0*(MM1/MM2)**(2.0D0*QQ))
     &             *S0**2/((BVTV**PP)*MM1**(2.0D0*QQ))**2
        FFFFT(1,3) = -(ZETA0*(MM1/MM3)**(2.0D0*QQ))
     &             *S0**2/((BVTV**PP)*MM1**(2.0D0*QQ))**2
        FFFFT(2,3) = -(ZETA0*(MM2/MM3)**(2.0D0*QQ))
     &             *S0**2/((BVTV**PP)*MM2**(2.0D0*QQ))**2
        FFFFT(2,1) = FFFFT(1,2)
        FFFFT(3,1) = FFFFT(1,3)
        FFFFT(3,2) = FFFFT(2,3)
        FFFFT(4,4) = 0.5D0/((TAUD0*(BVTV**PP)
     &             *(MM1*MM2)**QQ)**2)
        FFFFT(5,5) = 0.5D0/((TAUD0*(BVTV**PP)
     &             *(MM1*MM3)**QQ)**2)
        FFFFT(6,6) = 0.5D0/((TAUD0*(BVTV**PP)
     &             *(MM2*MM3)**QQ)**2)
C
C     QUADRIC SECOND ORDER TENSOR FF
        FFT(1,1) = -(SIGD0P-SIGD0N)/2.0D0/SIGD0P/SIGD0N
     &         /((BVTV**PP)*MM1**(2.0D0*QQ))
        FFT(2,2) = -(SIGD0P-SIGD0N)/2.0D0/SIGD0P/SIGD0N
     &         /((BVTV**PP)*MM2**(2.0D0*QQ))
        FFT(3,3) = -(SIGD0P-SIGD0N)/2.0D0/SIGD0P/SIGD0N
     &         /((BVTV**PP)*MM3**(2.0D0*QQ))
C
C       DUMMY SUPERPOSITION FOR ONLY TRABECULAR BONE
        SSSS = SSSST*PBVT
        CALL MIGS(SSSS, 6, CCCC)
        FFFF = FFFFT*PBVT
	    FFT  = FFT*PBVT
        FF(1) = FFT(1,1)
        FF(2) = FFT(2,2)
        FF(3) = FFT(3,3)
        FF(4) = 0.0D0
        FF(5) = 0.0D0
        FF(6) = 0.0D0
C
C_______________________________________________________________
C_____CORTICAL ELEMENTS_________________________________________
C_______________________________________________________________
C
C     CORTICAL SSSS, FFFF, FF
      ELSE IF (PBVC.GT.0.0D0.AND.PBVT.EQ.0.0D0) THEN
C     DENSITY-BASED TRANSVERSELY ISOTROPIC COMPACT BONE, MAIN DIRECTION 3
C     ___________________________________________________________________
C
C       BVTV was BVTVC*PBVC before. But as we now superimpose the stiffness matrices, I think we should change this to BVTV=BVTVC!
        BVTV = BVTVC
C
C       ELASTICITY PARAMETERS (TIRBCT, JJS, 2013) CORTEX
        E0   = 15079.6D0*(16000.0D0/24578.5D0)*SCA1
        EAA  = 24578.5D0*(16000.0D0/24578.5D0)*SCA1
        V0   = 0.4620D0
        VA0  = 0.354D0
        MU0 = 3117.0D0*(12000.0D0/9759.0D0)*SCA1
        MUA0 = 6578.01D0*(16000.0D0/24578.5D0)*SCA1
        KS   = 1.0D0
C
C       STRENGTH PARAMETERS (TIRBCT, JJS, 2013) CORTEX
        SIGD0P = 56.54D0*(16000.0D0/24578.5D0)*SCA2
        SIGD0N = 201.28D0*(16000.0D0/24578.5D0)*SCA2
        SIGDAP = 176.40D0*(16000.0D0/24578.5D0)*SCA2
        SIGDAN = 268.00D0*(16000.0D0/24578.5D0)*SCA2
        ZETA0  = 0.0074D0
        ZETAA0 = 1.4045D0
        TAUDA0 = 82.55D0*(16000.0D0/24578.5D0)*SCA2
        PP     = 1.0D0
C
C
C       TRANSVERSELY ISOTROPIC COMPLIANCE TENSOR CCCC
        MU0       = E0/2.0D0/(1.0D0+V0) 
        CCCCC(1,1) = 1.0D0/E0
        CCCCC(2,2) = 1.0D0/E0
        CCCCC(3,3) = 1.0D0/EAA
        CCCCC(4,4) = 0.5D0/MU0
        CCCCC(5,5) = 0.5D0/MUA0
        CCCCC(6,6) = 0.5D0/MUA0
        CCCCC(2,1) = -V0/E0
        CCCCC(3,1) = -VA0/EAA
        CCCCC(3,2) = -VA0/EAA
        CCCCC(1,2) = CCCCC(2,1)
        CCCCC(1,3) = CCCCC(3,1)
        CCCCC(2,3) = CCCCC(3,2)
        CCCCC      = CCCCC/BVTV**KS
C
C       TRANSVERSELY ISOTROPIC STIFFNESS TENSOR SSSS
        CALL MIGS(CCCCC,6,SSSSC)
C
C       TRANSVERSELY ISOTROPIC QUADRIC FOURTH ORDER TENSOR FFFF
        S0        = (SIGD0P+SIGD0N)/2.0D0/SIGD0P/SIGD0N
        SA        = (SIGDAP+SIGDAN)/2.0D0/SIGDAP/SIGDAN
        TAUD0     = DSQRT(0.5D0/S0**2/(1.0D0+ZETA0))
        FFFFC(1,1) = S0**2
        FFFFC(2,2) = S0**2
        FFFFC(3,3) = SA**2
        FFFFC(1,2) = -ZETA0*S0**2
        FFFFC(1,3) = -ZETAA0*SA**2
        FFFFC(2,3) = -ZETAA0*SA**2
        FFFFC(2,1) = FFFFC(1,2)
        FFFFC(3,1) = FFFFC(1,3)
        FFFFC(3,2) = FFFFC(2,3)
        FFFFC(4,4) = 0.5D0/TAUD0**2
        FFFFC(5,5) = 0.5D0/TAUDA0**2
        FFFFC(6,6) = 0.5D0/TAUDA0**2
        FFFFC      = FFFFC/BVTV**PP/BVTV**PP
C 
C       TRANSVERSELY ISOTROPIC QUADRIC SECOND ORDER TENSOR FF
        FFC(1,1) = -(SIGD0P-SIGD0N)/2.0D0/SIGD0P/SIGD0N
        FFC(2,2) = -(SIGD0P-SIGD0N)/2.0D0/SIGD0P/SIGD0N
        FFC(3,3) = -(SIGDAP-SIGDAN)/2.0D0/SIGDAP/SIGDAN
        FFC    = FFC/BVTV**PP
C        
C       DUMMY SUPERPOSITION FOR ONLY TRABECULAR BONE
C
        SSSS = SSSSC*PBVC
	    CALL MIGS(SSSS, 6, CCCC)
        FFFF = FFFFC*PBVC
	    FFC  = FFC*PBVC
        FF(1) = FFC(1,1)
        FF(2) = FFC(2,2)
        FF(3) = FFC(3,3)
        FF(4) = 0.0D0
        FF(5) = 0.0D0
        FF(6) = 0.0D0
      ENDIF
C
C     _______________________________________________________________
      ETOT1 = STRAN+DSTRAN
      ETOT0 = STRAN
C 
C     SET ALL STATE VARIABLES TO ZERO AT THE BEGINNING OF THE SIMULATION
      IF (KSTEP.EQ.(1).AND.KINC .EQ.(1)) THEN
        DO K1=1,18
          STATEV(K1) = 0.0D0
        END DO
      END IF
C     
      KAPPA0 = STATEV(1)
C
      IF (KAPPA0.LT.0.0D0) THEN
        KAPPA0 = 0.0D0
      END IF 
C
C     RECOVER AND ROTATE PLASTIC STRAIN 
      CALL ROTSIG(STATEV(3),DROT,EPLAS0,2,NDI,NSHR)
C 
C     CONVERT INITIAL STRAINS TO VOITH-MANDL NOTATION
      DO K1=4,6
        ETOT1(K1)  = ETOT1(K1)/DSQRT(2.0D0)
        ETOT0(K1)  = ETOT0(K1)/DSQRT(2.0D0)
        EPLAS0(K1) = EPLAS0(K1)/DSQRT(2.0D0)
      END DO
C
C     STRESS AT THE BEGINNING OF THE CURRENT INCREMENT
      SS0 = (1.0D0-DAM(KAPPA0,KONSTD,CRITD))
     &    *MATMUL(SSSS,ETOT0-EPLAS0)
C
C     ELASTIC TRIAL STRESS
      STR = (1.0D0-DAM(KAPPA0,KONSTD,CRITD))
     &    *MATMUL(SSSS,ETOT1-EPLAS0)
C 
C     YIELD CRITERION WITH TRIAL STRESS
      RAD  = RADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,KAPPA0,KMIN
     &     ,GMIN,ND,EXPS)
      DRAD = DRADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,KAPPA0,KMIN
     &     ,GMIN,ND,EXPS)
      FFS  = MATMUL(FFFF,STR)
      SFFS = DOT_PRODUCT(STR,FFS)
C 
      YSTR = DSQRT(SFFS)+DOT_PRODUCT(FF,STR)-RAD
C     _______________________________________________________________
C
C     ELASTIC CASE
C     _______________________________________________________________
      IF (YSTR.LE.TOL) THEN
C
C       UPDATE STATE VARIABLE
        KAPPA1 = KAPPA0
        EPLAS1 = EPLAS0 
C
C       ELASTIC STRESS UPDATE
        SS1 = STR
C
C       ELASTIC JACOBIAN
        TANM = (1.0D0-DAM(KAPPA0,KONSTD,CRITD))*SSSS
C     _______________________________________________________________
C
C     INELASTIC CASE
C     _______________________________________________________________
      ELSE IF (YSTR.GT.TOL) THEN
C    
        ITER    = 0
        EPLAS1  = EPLAS0
        DKAPPA1 = 0.0D0
        KAPPA1  = KAPPA0+DKAPPA1
        SS1     = STR
C 
        DMG  = DAM(KAPPA1,KONSTD,CRITD)
        DMG0 = DAM(KAPPA0,KONSTD,CRITD)
        DDMG = DDAM(KAPPA1,KONSTD,CRITD)
        RAD  = RADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,KAPPA1,KMIN,GMIN,ND
     &        ,EXPS)
        DRAD = DRADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,KAPPA1,KMIN,GMIN,ND
     &        ,EXPS)
C 
        FFS  = MATMUL(FFFF,SS1)
        SFFS = DOT_PRODUCT(SS1,FFS)
C
        DSY  = 1.0D0/DSQRT(SFFS)*FFS+FF
        DDSY = -1.0D0*(SFFS)**(-1.5D0)*VECDYAD(FFS,FFS)
     &       +1.0D0/DSQRT(SFFS)*FFFF
C
        HI = VECNORM(DSY)
        NP = DSY/HI
C 
        YS = DSQRT(SFFS)+DOT_PRODUCT(FF,SS1)-RAD
     &     +VISCY(VISCFL,ETA,MM,DTIME,DKAPPA1,HI)
        RR = -MATMUL(CCCC,SS1-STR)/(1.0D0-DMG)-(DMG-DMG0)/(1.0D0-DMG)
     &     *(ETOT1-EPLAS0)-DKAPPA1*NP
C   
        ABSY   = DABS(YS)
        NORMRR = VECNORM(RR)
C
C     _______________________________________________________________
C
C         NEWTON CLOSEST POINT PROJECTION W/O LINE SEARCH
C     _______________________________________________________________ 
        DO WHILE ((NORMRR.GT.TOL.OR.ABSY.GT.TOL).AND.ITER.LE.10)
C 
          ITER = ITER+1 
C 
          SSI     = SS1
          DKAPPAI = DKAPPA1
C 
          DHDS  = 1.0D0/DSQRT(DOT_PRODUCT(DSY,DSY))*MATMUL(DDSY,DSY)
          DHDK  = 0.0D0
          DYDS  = DSY+VDYDS(VISCFL,ETA,MM,DTIME,DKAPPAI,DHDS,HI)
          DYDK  = -DRAD+VDYDK(VISCFL,ETA,MM,DTIME,DKAPPAI,DHDK,HI)
          DNPDS = (DDSY*HI-VECDYAD(DSY,DHDS))/HI**2
          DNPDK = 0.0D0
C
          DRRDS = -CCCC/(1.0D0-DMG)-DKAPPAI*DNPDS
          DRRDK = -DDMG/(1.0D0-DMG)**2*(MATMUL(CCCC,SSI-STR)+(1.0D0
     &          -DMG0)*(ETOT1-EPLAS0))-NP-DKAPPAI*DNPDK
C 
          CALL MIGS(-DRRDS,6,SSSA)
C 
          DDK1 = -(YS/VECNORM(DYDS)+DOT_PRODUCT(NP,MATMUL(SSSA,RR)))
     &         /(DOT_PRODUCT(NP,MATMUL(SSSA,DRRDK))
     &         +DYDK/VECNORM(DYDS))
C 
          DSS1 = MATMUL(SSSA,(RR+DRRDK*DDK1))
C 
          SS1 = SSI+DSS1
C 
          DKAPPA1 = DKAPPAI+DDK1
C 
C         IF DKAPPA1 BECOMES SMALLER THAN 0 AND THEREFORE INADMISSIBLE,
C         RESTART THE NEWTON SCHEME WITH A DIFFERENT STARTING POINT
          IF (DKAPPA1.LT.0.0D0) THEN
            SS1     = STR
            DKAPPA1 = 0.0D0
          END IF
C 
          KAPPA1 = KAPPA0+DKAPPA1
C
C         CHECK RESIDUAL AND YIELD CRITERION AT I+1 FOR ITERATION BREAK
          DMG  = DAM(KAPPA1,KONSTD,CRITD)
          DDMG = DDAM(KAPPA1,KONSTD,CRITD)
          RAD  = RADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,KAPPA1,KMIN
     &         ,GMIN,ND,EXPS) 
          DRAD = DRADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,KAPPA1,KMIN
     &         ,GMIN,ND,EXPS) 
          FFS  = MATMUL(FFFF,SS1)
          SFFS = DOT_PRODUCT(SS1,FFS)
C
          DSY  = 1.0D0/DSQRT(SFFS)*FFS+FF
          DDSY = -1.0D0*(SFFS)**(-1.5D0)*VECDYAD(FFS,FFS)
     &         +1.0D0/DSQRT(SFFS)*FFFF
C
          HI = VECNORM(DSY)
          NP = DSY/HI
C 
          YS = DSQRT(SFFS)+DOT_PRODUCT(FF,SS1)-RAD
     &       +VISCY(VISCFL,ETA,MM,DTIME,DKAPPA1,HI)
          RR = -MATMUL(CCCC,SS1-STR)/(1.0D0-DMG)-(DMG-DMG0)
     &       /(1.0D0-DMG)*(ETOT1-EPLAS0)-DKAPPA1*NP
C 
          NORMRR = VECNORM(RR)
          ABSY   = DABS(YS)
C 
        END DO 
C     _______________________________________________________________
C
C         PRIMAL CPPA FOR GLOBAL CONVERGENCE
C     _______________________________________________________________ 
        IF (ITER.GT.10.AND.(NORMRR.GT.TOL.OR.ABSY.GT.TOL)) THEN 
C
          ITER    = 0
          EPLAS1  = EPLAS0
          DKAPPA1 = 0.0D0
          KAPPA1  = KAPPA0+DKAPPA1
          SS1     = STR
C 
          DMG  = DAM(KAPPA1,KONSTD,CRITD)
          DMG0 = DAM(KAPPA0,KONSTD,CRITD)
          DDMG = DDAM(KAPPA1,KONSTD,CRITD)
          RAD  = RADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,KAPPA1,KMIN
     &         ,GMIN,ND,EXPS)
          DRAD = DRADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,KAPPA1,KMIN
     &         ,GMIN,ND,EXPS)
C 
          FFS  = MATMUL(FFFF,SS1)
          SFFS = DOT_PRODUCT(SS1,FFS)
C
          DSY  = 1.0D0/DSQRT(SFFS)*FFS+FF
          DDSY = -1.0D0*(SFFS)**(-1.5D0)*VECDYAD(FFS,FFS)
     &         +1.0D0/DSQRT(SFFS)*FFFF
C
          HI = VECNORM(DSY)
          NP = DSY/HI
C 
          YS = DSQRT(SFFS)+DOT_PRODUCT(FF,SS1)-RAD
     &       +VISCY(VISCFL,ETA,MM,DTIME,DKAPPA1,HI)
          RR = -MATMUL(CCCC,SS1-STR)/(1.0D0-DMG)-(DMG-DMG0)
     &       /(1.0D0-DMG)*(ETOT1-EPLAS0)-DKAPPA1*NP
C   
          ABSY   = DABS(YS)
          NORMRR = VECNORM(RR)
C
          DO K1 = 1,6
            RRI(K1) = RR(K1)
          END DO
            RRI(7) = YS
C
          DO WHILE ((NORMRR.GT.TOL.OR.ABSY.GT.TOL).AND.ITER.LE.MAXITER)
C 
            ITER = ITER+1 
C 
            SSI     = SS1
            DKAPPAI = DKAPPA1
C 
            DHDS  = 1.0D0/DSQRT(DOT_PRODUCT(DSY,DSY))*MATMUL(DDSY,DSY)
            DHDK  = 0.0D0
            DYDS  = DSY+VDYDS(VISCFL,ETA,MM,DTIME,DKAPPAI,DHDS,HI)
            DYDK  = -DRAD+VDYDK(VISCFL,ETA,MM,DTIME,DKAPPAI,DHDK,HI)
            DNPDS = (DDSY*HI-VECDYAD(DSY,DHDS))/HI**2
            DNPDK = 0.0D0
C
            DRRDS = -CCCC/(1.0D0-DMG)-DKAPPAI*DNPDS
            DRRDK = -DDMG/(1.0D0-DMG)**2*(MATMUL(CCCC,SSI-STR)+
     &            (1.0D0-DMG0)*(ETOT1-EPLAS0))-NP-DKAPPAI*DNPDK
C 
C           BUILD SYSTEM RESIDUAL VECTOR
            DO K1 = 1,6
              RRR(K1) = RR(K1)
            END DO
            RRR(7) = YS
C
            J11 = DRRDS
            J22 = DYDK
            J12 = DRRDK
            J21 = DYDS
C
C           BUILD SYSTEM JACOBIAN
            DO K1 = 1,6
              DO K2 = 1,6
                JJJ(K1,K2) = J11(K1,K2)
              END DO
              JJJ(7,K1) = J21(K1)
              JJJ(K1,7) = J12(K1)
            END DO
            JJJ(7,7) = J22
C
C           INVERT SYSTEM JACOBIAN
            CALL MIGS(JJJ,7,INVJ)
C
            AUX = 0.0D0
            IF (DKAPPA1.EQ.0.0D0) THEN
              AUX = DOT_PRODUCT(NP,(SS1-STR))
            END IF
C
C           DETERMINE SYSTEM UPDATES
            IF (AUX.LE.0.0D0) THEN
C
              FLAG = 1
C
              DDD=MATMUL(-INVJ,RRR)
C
            ELSE
C
              FLAG = 0
C
              DDJ = MATMUL(INVJ,TRANSPOSE(INVJ))
              DO K1  =  1,6
                DDJ(7,K1) = 0.0D0
                DDJ(K1,7) = 0.0D0
              END DO
C
              DDD = MATMUL(-MATMUL(DDJ,TRANSPOSE(JJJ)),RRR)
C
            END IF  
C
C           RECOVER STRESS AND KAPPA UPDATE
            DO K1 = 1,6
              DSS1(K1) = DDD(K1)
            END DO
            DDK1 = DDD(7)            
C 
            SS1 = SSI+DSS1
            DKAPPA1 = DKAPPAI+DDK1
C
            DO K1 = 1,6
              XXI(K1) = SS1(K1)
            END DO
            XXI(7) = DKAPPA1
C 
C           ENFORCE CONSTRAIN DKAPPA >= 0.0D0
            IF (DKAPPA1 .LT. 0.0D0) THEN
              DKAPPA1 = 0.0D0
            END IF
C 
            KAPPA1 = KAPPA0+DKAPPA1
C
C           CHECK RESIDUAL AND YIELD CRITERION AT I+1 FOR ITERATION BREAK
            DMG  = DAM(KAPPA1,KONSTD,CRITD)
            DDMG = DDAM(KAPPA1,KONSTD,CRITD)
            RAD  = RADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,KAPPA1,KMIN
     &           ,GMIN,ND,EXPS)
            DRAD = DRADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,KAPPA1,KMIN
     &           ,GMIN,ND,EXPS)
            FFS  = MATMUL(FFFF,SS1)
            SFFS = DOT_PRODUCT(SS1,FFS)
C
            DSY  = 1.0D0/DSQRT(SFFS)*FFS+FF
            DDSY = -1.0D0*(SFFS)**(-1.5D0)*VECDYAD(FFS,FFS)
     &           +1.0D0/DSQRT(SFFS)*FFFF
C
            HI = VECNORM(DSY)
            NP = DSY/HI
C 
C           LINE SEARCH MERIT FUNCTION AND DERIVATIVE OF ITERATION I
            MKI = 0.5D0*(DOT_PRODUCT(RR,RR)+YS**2)
C
            IF (FLAG.EQ.1) THEN
              DMK = -2.0D0*MKI
            ELSE
              DMK = DOT_PRODUCT(RRR,MATMUL(JJJ,DDD))
            END IF
C           
C           RESIDUAL AND LINE SEARCH MERIT FUNCTION OF ITERATION I+1
            YS = DSQRT(SFFS)+DOT_PRODUCT(FF,SS1)-RAD
     &         +VISCY(VISCFL,ETA,MM,DTIME,DKAPPA1,HI)
            RR = -MATMUL(CCCC,SS1-STR)/(1.0D0-DMG)-(DMG-DMG0)
     &         /(1.0D0-DMG)*(ETOT1-EPLAS0)-DKAPPA1*NP
C
            DO K1 = 1,6
              RRR(K1) = RR(K1)
            END DO
            RRR(7) = YS
C
            MK1 = 0.5D0*(DOT_PRODUCT(RR,RR)+YS**2)
C
C           LINE SEARCH ALGORITHM FOR CONSTRAINED PROBLEMS (ARMERO 2002)
            IF (FLAG.EQ.1.AND.DKAPPA1.GE.0.0D0) THEN
              UPLIM = (1.0D0-2.0D0*BETA*ALPHA)*MKI
            ELSE
              UPLIM = MKI+BETA*DOT_PRODUCT(RRR,MATMUL(JJJ,XX1-XXI))
            END IF
C
            IF (MK1.GT.UPLIM) THEN
C
              ITERL = 0
              ALPHA = 1.0D0
              XX1 = XXI
C
              DO WHILE (MK1.GT.UPLIM.AND.ALPHA.GE.BETA.AND.ITERL.LT.20)
C
                ITERL = ITERL+1
C
                ALPHA1 = ETAL*ALPHA
                ALPHA2 = -ALPHA**2*DMK/2.0D0/(MK1-MKI-ALPHA*DMK)
C
                IF (ALPHA1.GE.ALPHA2) THEN
                  ALPHA = ALPHA1
                ELSE
                  ALPHA = ALPHA2
                END IF
C
                DKAPPA1 = DKAPPAI+ALPHA*DDK1
C
                IF (DKAPPA1.LT.0.0D0) THEN
                  DKAPPA1 = 0.0D0
                END IF

                KAPPA1 = KAPPA0+DKAPPA1
                SS1    = SSI+ALPHA*DSS1
C
                DO K1 = 1,6
                  XX1(K1) = SS1(K1)
                END DO
                XX1(7) = DKAPPA1
C
                DMG  = DAM(KAPPA1,KONSTD,CRITD)
                RAD  = RADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,KAPPA1,KMIN
     &               ,GMIN,ND,EXPS)
                FFS  = MATMUL(FFFF,SS1)
                SFFS = DOT_PRODUCT(SS1,FFS)
                DSY  = 1.0D0/DSQRT(SFFS)*FFS+FF
                HI   = VECNORM(DSY)
                NP   = DSY/HI
C   
                YS = DSQRT(SFFS)+DOT_PRODUCT(FF,SS1)-RAD
     &             +VISCY(VISCFL,ETA,MM,DTIME,DKAPPA1,HI)
                RR = -MATMUL(CCCC,SS1-STR)/(1.0D0-DMG)-(DMG-DMG0)
     &             /(1.0D0-DMG)*(ETOT1-EPLAS0)-DKAPPA1*NP
C
                DO K1 = 1,6
                  RRR(K1) = RR(K1)
                END DO
                RRR(7) = YS
C
                MK1 = 0.5D0*(DOT_PRODUCT(RR,RR)+YS**2)
C
                IF (FLAG.EQ.1.AND.DKAPPA1.GE.0.0D0) THEN
                  UPLIM = (1.0D0-2.0D0*BETA*ALPHA)*MKI
                ELSE
                  UPLIM = MKI+BETA
     &                  *DOT_PRODUCT(RRI,MATMUL(JJJ,XX1-XXI))
                END IF
C       
              END DO
C
              DMG  = DAM(KAPPA1,KONSTD,CRITD)
              DDMG = DDAM(KAPPA1,KONSTD,CRITD)
              RAD  = RADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,KAPPA1,KMIN
     &             ,GMIN,ND,EXPS)
              DRAD = DRADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,KAPPA1,KMIN
     &             ,GMIN,ND,EXPS)
              FFS  = MATMUL(FFFF,SS1)
              SFFS = DOT_PRODUCT(SS1,FFS)
C
              DSY  = 1.0D0/DSQRT(SFFS)*FFS+FF
              DDSY = -1.0D0*(SFFS)**(-1.5D0)*VECDYAD(FFS,FFS)
     &             +1.0D0/DSQRT(SFFS)*FFFF
C
              HI = VECNORM(DSY)
              NP = DSY/HI
C         
              YS = DSQRT(SFFS)+DOT_PRODUCT(FF,SS1)-RAD
     &           +VISCY(VISCFL,ETA,MM,DTIME,DKAPPA1,HI)
              RR = -MATMUL(CCCC,SS1-STR)/(1.0D0-DMG)-(DMG-DMG0)
     &           /(1.0D0-DMG)*(ETOT1-EPLAS0)-DKAPPA1*NP
C
              DO K1 = 1,6
                RRI(K1) = RR(K1)
              END DO
              RRI(7) = YS
C
            END IF
C 
            NORMRR = VECNORM(RR)
            ABSY = DABS(YS)
C 
          END DO
        END IF
C 
C       TANGENT STIFFNESS MATRIX
        DMG  = DAM(KAPPA1,KONSTD,CRITD)
        DDMG = DDAM(KAPPA1,KONSTD,CRITD)
        RAD  = RADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,KAPPA1,KMIN,GMIN
     &       ,ND,EXPS)
        DRAD = DRADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,KAPPA1,KMIN
     &       ,GMIN,ND,EXPS)
        FFS  = MATMUL(FFFF,SS1)
        SFFS = DOT_PRODUCT(SS1,FFS)
C
        DSY  = 1.0D0/DSQRT(SFFS)*FFS+FF
        DDSY = -1.0D0*(SFFS)**(-1.5D0)*VECDYAD(FFS,FFS)
     &       +1.0D0/DSQRT(SFFS)*FFFF
C
        HI = VECNORM(DSY)
        NP = DSY/HI
C         
        YS = DSQRT(SFFS)+DOT_PRODUCT(FF,SS1)-RAD
     &     +VISCY(VISCFL,ETA,MM,DTIME,DKAPPA1,HI)
        RR = -MATMUL(CCCC,SS1-STR)/(1.0D0-DMG)-(DMG-DMG0)
     &     /(1.0D0-DMG)*(ETOT1-EPLAS0)-DKAPPA1*NP
C
        DHDS  = 1.0D0/DSQRT(DOT_PRODUCT(DSY,DSY))*MATMUL(DDSY,DSY)
        DHDK  = 0.0D0
        DYDS  = DSY+VDYDS(VISCFL,ETA,MM,DTIME,DKAPPAI,DHDS,HI)
        DYDK  = -DRAD+VDYDK(VISCFL,ETA,MM,DTIME,DKAPPAI,DHDK,HI)
        DNPDS = (DDSY*HI-VECDYAD(DSY,DHDS))/HI**2
        DNPDK = 0.0D0
C
        DRRDS = -CCCC/(1.0D0-DMG)-DKAPPAI*DNPDS
        DRRDK = -DDMG/(1.0D0-DMG)**2*(MATMUL(CCCC,SSI-STR)+(1.0D0
     &        -DMG0)*(ETOT1-EPLAS0))-NP-DKAPPAI*DNPDK
C 
        CALL MIGS(-DRRDS,6,SSSA)
C 
        TANM = SSSA-MATMUL(MATMUL(SSSA,VECDYAD(DRRDK,NP)),SSSA)
     &       /(DOT_PRODUCT(NP,MATMUL(SSSA,DRRDK))+DYDK/VECNORM(DYDS))
C 
C       FINAL PLASTIC STRAIN
        EPLAS1 = ETOT1-MATMUL(CCCC,SS1)/(1.0D0
     &         -DAM(KAPPA1,KONSTD,CRITD))
C 
      END IF
C
C     CHECK FOR CONVERGENCE OF BACK PROJECTION ALGORITHM
      IF (YSTR.GE.0.0D0.AND.ITER.GT.MAXITER) THEN
        WRITE(*,*) 'HIGH NUMBER OF ITERATIONS NEEDED'
        WRITE(*,*) 'Step : ',KSTEP
        WRITE(*,*) 'Inc : ',KINC
        WRITE(*,*) 'Elem : ',NOEL
        WRITE(*,*) 'Point : ',NPT
        WRITE(*,*) 'THE INCREMENT SIZE WILL BE REDUCED BY 50%'
        PNEWDT = 0.5
        KAPPA1 = KAPPA0
        EPLAS1 = EPLAS0
        SS1    = STR
        TANM   = (1.0D0-DAM(KAPPA0,KONSTD,CRITD))*SSSS
      END IF
C
      DO K1 = 1,6
        STATEV(8+K1) = SS1(K1)
      END DO
C     _______________________________________________________________
C
C     RETURN VARIABLES TO ABAQUS 
C     _______________________________________________________________
C
C     DISSIPATION (PLASTICITY AND DAMAGE)
      SPD = SPD+DOT_PRODUCT((SS0+SS1)/2.0D0,EPLAS1-EPLAS0)+0.5D0*
     &    (0.5D0*(DDAM(KAPPA1,KONSTD,CRITD)
     &    +DDAM(KAPPA0,KONSTD,CRITD)))
     &    *DOT_PRODUCT((0.5D0*(ETOT1-EPLAS1+ETOT0-EPLAS0)),
     &    MATMUL(SSSS,(0.5D0*(ETOT1-EPLAS1+ETOT0-EPLAS0))))
     &    *(KAPPA1-KAPPA0)
C
C     CREEP DISSIPATION
      SCD = 0.0D0
C
C     CONVERT TANGENT STIFFNESS OPERATOR MATRIX TO ABAQUS CONVENTION 
      DO K1 = 4,6
        DO K2 = 4,6
          TANM(K1,K2)   = TANM(K1,K2)/2.0D0
          TANM(K1-3,K2) = TANM(K1-3,K2)/DSQRT(2.0D0)
          TANM(K1,K2-3) = TANM(K1,K2-3)/DSQRT(2.0D0)
        END DO
      END DO
C     _______________________________________________________________
C
C     CONVERT STRESSES AND STRAINS TO ABAQUS CONVENTION 
C     _______________________________________________________________
      DO K1=4,6
        SS1(K1)    = SS1(K1)/DSQRT(2.0D0)
        EPLAS1(K1) = EPLAS1(K1)*DSQRT(2.0D0)
        EPLAS0(K1) = EPLAS0(K1)*DSQRT(2.0D0)
      END DO 
C
C    _______________________________________________________________    
C
C     COMPUTE OF VALUE OF OPTIMIZATION FUNCTION
C     _______________________________________________________________
C
      OF = (1.0D0/EPS0**2)*DOTP(ETOT1,ETOT1)+
     &    (2.0D0/EPS0)*(ETOT1(1)+ETOT1(2)+ETOT1(3))+3.0D0
      
C     _______________________________________________________________    
C
C     UPDATE OF FIELD VARIABLES
C     _______________________________________________________________
C 
C     STRESS TENSOR
      STRESS = SS1
C    
C     TANGENT STIFFNESS OPERATOR 
      DDSDDE = TANM
C
C     STRAIN ENERGY
      SSE = SSE+DOT_PRODUCT((SS0+SS1)/2.0D0,DSTRAN)
C
C     PLASTIC STRAIN TENSOR
      DO K1 = 3,8
        STATEV(K1) = EPLAS1(K1-2)
      END DO
C
C     CUMULATED PLASTIC STRAIN      
      STATEV(1) = KAPPA1
C
C     DAMAGE
      STATEV(2) = DAM(KAPPA1,KONSTD,CRITD)
C
C     BVTV and PBV
      STATEV(15) = BVTVC
      STATEV(16) = BVTVT
      STATEV(17) = PBVC
      STATEV(18) = PBVT
      STATEV(19) = RDY
      STATEV(20) = MM3
      STATEV(21) = (BVTVC*PBVC+BVTVT*PBVT)*MM
      STATEV(22) = OF
C      STATEV(22) = SSSS(3,3)

C 
      IF (DKAPPA1.LT.0.0D0) THEN
        WRITE(*,*) 'Inadmissible DKappa occured'
        WRITE(*,*) 'Step : ',KSTEP
        WRITE(*,*) 'Inc : ',KINC
        WRITE(*,*) 'Elem : ',NOEL
        WRITE(*,*) 'Point : ',NPT
        WRITE(*,*) 'KAPPA1 : ',KAPPA1
        WRITE(*,*) 'DKAPPA1 : ',DKAPPA1
        WRITE(*,*) 'DAM(KAPPA) : ',STATEV(2)
      END IF
C 
      RETURN
      END
C     _______________________________________________________________
C
C     MATERIAL MODEL FUNCTIONS
C     _______________________________________________________________
C
C     DAMAGE FUNCTION DAM(k)
      DOUBLE PRECISION FUNCTION DAM(KAPPA,KONSTD,CRITD)
      IMPLICIT NONE
      DOUBLE PRECISION KAPPA,KONSTD,CRITD
C     MEAN KONSTK OF WOLFRAM (JBIOMECH 2011) FOR COMPRESSION AND 
C     TENSION, AXIAL AND TRANSVERSE: 8.0075, CRITD = 0.85
      DAM=CRITD*(1.0D0-DEXP(-KONSTD*KAPPA))
      RETURN
      END
C
C     DERIVATIVE OF DAMAGE FUNCTION dDAM/dk
C 
      DOUBLE PRECISION FUNCTION DDAM(KAPPA,KONSTD,CRITD)
      IMPLICIT NONE
      DOUBLE PRECISION KAPPA,KONSTD,CRITD
C     MEAN KONSTK OF WOLFRAM (JBIOMECH 2011) FOR COMPRESSION AND 
C     TENSION, AXIAL AND TRANSVERSE: 8.0075, CRITD = 0.85
      DDAM = CRITD*KONSTD*DEXP(-KONSTD*KAPPA)
      RETURN
      END
C
C     POSTYIELD FUNCTION r(k)
      DOUBLE PRECISION FUNCTION RADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,
     &                              KAPPA,KMIN,GMIN,ND,EXPS)
      IMPLICIT NONE
      DOUBLE PRECISION KSLOPE,KMAX,KWIDTH,RDY,KAPPA,OFFS
      DOUBLE PRECISION KMIN,GMIN,ND,EXPS
      INTEGER PYFL
C 
C     PERFECT PLASTICITY
      IF (PYFL.EQ.0) THEN
        RADK = 1.0D0
C     EXPONENTIAL HARDENING
      ELSEIF (PYFL.EQ.1) THEN
        RADK = RDY+(1.0D0-RDY)*(1.0D0-DEXP(-KSLOPE*KAPPA))
C     SIMPLE SOFTENING FUNCTION
      ELSEIF (PYFL.EQ.2) THEN
        OFFS = 1.0D0/KWIDTH 
        RADK = RDY+(1.0D0-RDY)*(DEXP(-((KAPPA-KMAX)**2)
     &       /(KWIDTH*KMAX**2))-DEXP(-OFFS-KSLOPE*KAPPA))
C     CONTINUOUS SOFTENING FUNCTION WITH GMIN, KSLOPE = 1000.0
      ELSEIF (PYFL.EQ.3) THEN
        OFFS   = -DLOG(DEXP(-KMAX**2/(KMAX**2+0.1D0*KMIN**2))
     &         +GMIN*(1.0D0-DEXP(-KMAX**2/(KMAX**2+0.1D0*KMIN**2))))
        KSLOPE = 1000.0D0
        RADK   = RDY+(1.0D0-RDY)*(DEXP(-((KAPPA-KMAX)**2)
     &         /(KMAX**2+0.1D0*KMIN**2))
     &         +GMIN*(1.0D0-DEXP(-((KAPPA-KMAX)**2)/(KMAX**2
     &         +0.1D0*KMIN**2)))-DEXP(-OFFS-KSLOPE*KAPPA))
C     PIECEWISE SOFTENING FUNCTION
      ELSEIF (PYFL.EQ.4) THEN
        IF (KAPPA.LT.KMAX) THEN
          RADK = 1.0D0-((KMAX-KAPPA)/KMAX)**(EXPS*KMAX)
        ELSEIF ((KAPPA.GE.KMAX).AND.(KAPPA.LT.(KMIN+KMAX)/2.0D0)) THEN
          RADK = 1.0D0-((1.0D0-GMIN)/2.0D0)*((2.0D0*(KAPPA-KMAX))
     &         /(KMIN-KMAX))**ND
        ELSEIF ((KAPPA.GE.(KMIN+KMAX)/2.0D0).AND.(KAPPA.LT.KMIN)) THEN
          RADK = GMIN+((1.0D0-GMIN)/2.0D0)*((2.0D0*(KMIN-KAPPA))
     &         /(KMIN-KMAX))**ND
        ELSE 
          RADK = GMIN
        END IF
        RADK = RDY+(1.0D0-RDY)*RADK
C     SIMPLE SOFTENING CORRECTED
      ELSEIF (PYFL.EQ.5) THEN
	OFFS = 1.0D0/KWIDTH 
        RADK = RDY+(1.0D0-RDY)*(DEXP(-((KAPPA-KMAX)**2)
     &       /(KWIDTH*KMAX**2))-DEXP(-OFFS-KSLOPE*KAPPA/(1-RDY)))
      END IF
C
      RETURN
      END
C
C     DERIVATIVE OF POSTYIELD FUNCTION dr/dk
      DOUBLE PRECISION FUNCTION DRADK(PYFL,KSLOPE,KMAX,KWIDTH,RDY,
     &                               KAPPA,KMIN,GMIN,ND,EXPS)
      IMPLICIT NONE
      DOUBLE PRECISION KSLOPE,KMAX,KWIDTH,RDY,KAPPA,OFFS
      DOUBLE PRECISION KMIN,GMIN,ND,EXPS
      INTEGER PYFL
C
C     PERFECT PLASTICITY
      IF (PYFL.EQ.0) THEN
        DRADK = 0.0D0
C     EXPONENTIAL HARDENING
      ELSEIF (PYFL.EQ.1) THEN
        DRADK = (1.0D0-RDY)*KSLOPE*DEXP(-KSLOPE*KAPPA)
C     SIMPLE SOFTENING FUNCTION WITH 
      ELSEIF (PYFL.EQ.2) THEN
        OFFS  = 1.0D0/KWIDTH 
        DRADK = (1.0D0-RDY)*((-KAPPA+KMAX)*DEXP(-(KAPPA-KMAX)**2
     &        /(KWIDTH*KMAX**2))/(0.5D0*KWIDTH*KMAX**2)
     &        +KSLOPE*DEXP(-OFFS-KSLOPE*KAPPA))
C     CONTINUOUS SOFTENING FUNCTION WITH GMIN, KSLOPE = 1000.0
      ELSEIF (PYFL.EQ.2) THEN
        OFFS   = -DLOG(EXP(-KMAX**2/(KMAX**2+0.1D0*KMIN**2))+GMIN*
     &         (1.0D0-DEXP(-KMAX**2/(KMAX**2+0.1D0*KMIN**2))))
        KSLOPE = 1000.0D0
        DRADK  = (1.0D0-RDY)*(-2.0D0*(KAPPA-KMAX)
     &         /(KMAX**2+0.1D0*KMIN**2)
     &         *DEXP(-((KAPPA-KMAX)**2)/(KMAX**2+0.1D0*KMIN**2))
     &         +2.0D0*GMIN*(KAPPA-KMAX)/(KMAX**2+0.1D0*KMIN**2)
     &         *DEXP(-((KAPPA-KMAX)**2)/(KMAX**2+0.1D0*KMIN**2))
     &         +KSLOPE*DEXP(-OFFS-KSLOPE*KAPPA))
C     PIECEWISE SOFTENING FUNCTION
      ELSEIF (PYFL.EQ.4) THEN
        IF (KAPPA.LT.KMAX) THEN
          DRADK = EXPS*((KMAX-KAPPA)/KMAX)**(EXPS*KMAX-1.0D0)
        ELSEIF ((KAPPA.GE.KMAX).AND.(KAPPA.LT.(KMIN+KMAX)/2.0D0))THEN
          DRADK = ND*(GMIN-1.0D0)/(KMIN-KMAX)*((2.0D0*(KAPPA-KMAX))
     &          /(KMIN-KMAX))**(ND-1.0D0)
        ELSEIF ((KAPPA.GE.(KMIN+KMAX)/2.0D0).AND.(KAPPA.LT.KMIN))THEN
          DRADK = ND*(1.0D0-GMIN)/(KMIN-KMAX)*((2.0D0*(KMIN-KAPPA))
     &          /(KMIN-KMAX))**(ND-1.0D0)
        ELSE 
          DRADK = 0.0D0
        END IF
        DRADK = (1.0D0-RDY)*DRADK
C     SIMPLE SOFTENING CORRECTED
      ELSEIF (PYFL.EQ.5) THEN
        OFFS  = 1.0D0/KWIDTH 
        DRADK = (1.0D0-RDY)*((-KAPPA+KMAX)*DEXP(-(KAPPA-KMAX)**2
     &        /(KWIDTH*KMAX**2))/(0.5D0*KWIDTH*KMAX**2)
     &        +KSLOPE/(1-RDY)*DEXP(-OFFS-KSLOPE*KAPPA))
      END IF
C 
      RETURN
      END
C
C     VISCOSITY FUNCTION
      DOUBLE PRECISION FUNCTION VISCY(VISCFL,ETA,MM,DTIME,DKAPPA,HI)
      IMPLICIT NONE
      DOUBLE PRECISION ETA,MM,DTIME,DKAPPA,HI
      INTEGER VISCFL
C     RATE-INDEPENDENT
      IF (VISCFL.EQ.0) THEN
        VISCY = 0.0D0
C     LINEAR FLOW RULE
      ELSEIF (VISCFL.EQ.1) THEN
        VISCY = -ETA/DTIME*DKAPPA/HI
C     EXPONENTIAL FLOW RULE
      ELSEIF (VISCFL.EQ.2) THEN
        VISCY = -LOG(1.0D0+ETA/DTIME*DKAPPA/HI)/MM
C     LOGARITHMIC FLOW RULE
      ELSEIF (VISCFL.EQ.3) THEN
        VISCY = -(EXP(ETA/DTIME*DKAPPA/HI)-1.0D0)/MM
C     POLYNOMIAL FLOW RULE
      ELSEIF (VISCFL.EQ.4) THEN
        VISCY = 0.5D0*MM-DSQRT((MM**2)/4.0D0+ETA/DTIME*DKAPPA/HI)
C     POWER LAW FLOW RULE
      ELSEIF (VISCFL.EQ.5) THEN
        VISCY = -SIGN(1.0D0,DKAPPA)*(ETA/DTIME
     &        *DABS(DKAPPA)/HI)**(1.0D0/MM)
      END IF
      RETURN
      END
C
      FUNCTION VDYDS(VISCFL,ETA,MM,DTIME,DKAPPA,DHDS,HI)
      IMPLICIT NONE
      DOUBLE PRECISION VDYDS(6)
      DOUBLE PRECISION ETA,MM,DTIME,DKAPPA,DHDS(6),HI
      INTEGER VISCFL
C     RATE-INDEPENDENT
      IF (VISCFL.EQ.0) THEN
        VDYDS = 0.0D0
C     LINEAR FLOW RULE
      ELSEIF (VISCFL.EQ.1) THEN
        VDYDS = 1.0D0/HI**2*DHDS*ETA/DTIME*DKAPPA
C     EXPONENTIAL FLOW RULE
      ELSEIF (VISCFL.EQ.2) THEN
        VDYDS = ETA/DTIME*DKAPPA/MM*DHDS/(HI**2)
     &        *1.0D0/(1.0D0+ETA/DTIME*DKAPPA/HI)
C     LOGARITHMIC FLOW RULE
      ELSEIF (VISCFL.EQ.3) THEN
        VDYDS = ETA/DTIME*DKAPPA/MM*DHDS/(HI**2)
     &        *EXP(ETA/DTIME*DKAPPA/HI)
C     POLYNOMIAL FLOW RULE
      ELSEIF (VISCFL.EQ.4) THEN
        VDYDS = 0.5D0*DHDS*ETA/DTIME*DABS(DKAPPA)/(HI**2)
     &        /DSQRT((MM**2)/4.0D0+ETA/DTIME*DKAPPA/HI)
C     POWER LAW FLOW RULE
      ELSEIF (VISCFL.EQ.5) THEN
        VDYDS = 1.0D0/MM/HI*DHDS*SIGN(1.0D0,DKAPPA)
     &        *(ETA/DTIME*DABS(DKAPPA)/HI)**(1.0D0/MM)
      END IF
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION VDYDK(VISCFL,ETA,MM,DTIME,DKAPPA,
     &                          DHDK,HI)
      IMPLICIT NONE
      DOUBLE PRECISION ETA,MM,DTIME,DKAPPA,DHDK,HI
      INTEGER VISCFL
C     RATE-INDEPENDENT
      IF (VISCFL.EQ.0) THEN
        VDYDK = 0.0D0
C     LINEAR FLOW RULE
      ELSEIF (VISCFL.EQ.1) THEN
        VDYDK = -(HI-DKAPPA*DHDK)/HI**2*ETA/DTIME
C     EXPONENTIAL FLOW RULE
      ELSEIF (VISCFL.EQ.2) THEN
        VDYDK = -ETA/DTIME/MM*(HI-DKAPPA*DHDK)/HI**2
     &        *1.0D0/(1.0D0+ETA/DTIME*DKAPPA/HI)
C     LOGARITHMIC FLOW RULE
      ELSEIF (VISCFL.EQ.3) THEN
        VDYDK = -ETA/DTIME/MM*(HI-DKAPPA*DHDK)/HI**2
     &        *EXP(ETA/DTIME*DKAPPA/HI)
C     POLYNOMIAL FLOW RULE
      ELSEIF (VISCFL.EQ.4) THEN
        VDYDK = -0.5D0*(HI-DKAPPA*DHDK)/HI**2*(ETA/DTIME)
     &        /DSQRT((MM**2)/4.0D0+ETA/DTIME*DKAPPA/HI)
C     POWER LAW FLOW RULE
      ELSEIF (VISCFL.EQ.5) THEN
        VDYDK = -1.0D0/MM*SIGN(1.0D0,DKAPPA)*(SIGN(1.0D0,DKAPPA)*HI
     &        -DABS(DKAPPA)*DHDK)/HI**2*(ETA/DTIME)**(1.0D0/MM)
     &        *(DABS(DKAPPA)/HI)**(1.0D0/MM-1.0)
      END IF
      RETURN
      END
C     _______________________________________________________________
C
C     MATHEMATICAL FUNCTIONS
C     _______________________________________________________________
C     
C     INTRINSIC FUNCTIONS:
C
C     MATMUL(AAAA,BBBB)
C     MULTIPLICATION MATRIX WITH MATRIX
C
C     MATMUL(AAAA,BB)
C     MULTIPLICATION MATRIX WITH VECTOR
C
C     DOT_PRODUCT(AA,BB)
C     SCALAR PRODUCT OF TWO VECTORS
C 
C     EXTERNAL FUNCTIONS:
C     _______________________________________________________________
C 
      FUNCTION VECDYAD(BB,CC)
      IMPLICIT NONE
C     DYADIC PRODUCT AB OF VECTORS BB AND CC WITH DIMENSION 6x1  
      DOUBLE PRECISION BB(6),CC(6)
      DOUBLE PRECISION VECDYAD(6,6)
      INTEGER I,J
C 
      VECDYAD = 0.0D0
C 
      VECDYAD = SPREAD(BB,2,SIZE(CC))*SPREAD(CC,1,SIZE(BB))
C 
      RETURN
      END
C      
C     _______________________________________________________________
C
      FUNCTION DOTP(DD,EE)
      IMPLICIT NONE
C      DOT PRODUCT DE OF VECTORS DD AND EE WITH DIMENSTION 6X1
      DOUBLE PRECISION DD(6), EE(6)
      DOUBLE PRECISION DOTP
      DOTP = 0.0D0
      DOTP = DD(1)*EE(1)+DD(2)*EE(2)+DD(3)*EE(3)+DD(4)*EE(4)+DD(5)*EE(5)+DD(6)*EE(6)
      END FUNCTION DOTP
      
C     _______________________________________________________________
C 
      DOUBLE PRECISION FUNCTION VECNORM(DD)
      IMPLICIT NONE
C     NORM OF VECTOR AA WITH DIMENSION 6x1
      DOUBLE PRECISION DD(6)
C 
      VECNORM = DSQRT(DOT_PRODUCT(DD,DD))
C 
      RETURN
      END
C     _______________________________________________________________
C 
C     SUBROUTINES FOR MATRIX INVERSION
C     _______________________________________________________________
C 
      SUBROUTINE MIGS(A,N,X)
C
C     Subroutine to invert matrix A(N,N) with the inverse stored
C     in X(N,N) in the output. A is stored in STOA and is resituted
C     as outpout
C 
      IMPLICIT NONE
C 
      DOUBLE PRECISION A(N,N),STOA(N,N),B(N,N),X(N,N)
      INTEGER INDX(N),N,I,J,K
C 
      DO 140 I = 1,N
        DO 130 J = 1,N
          STOA(I,J) = A(I,J)
          B(I,J)    = 0.0D0
  130   CONTINUE
  140 CONTINUE
      DO 150 I = 1,N
        B(I,I) = 1.0D0
  150 CONTINUE
C 
      CALL ELGS(A,N,INDX)
C 
      DO 180 I = 1,N-1
        DO 170 J = I+1,N
          DO 160 K = 1,N
            B(INDX(J),K) = B(INDX(J),K)
     &                    -A(INDX(J),I)*B(INDX(I),K)
  160     CONTINUE
  170   CONTINUE
  180 CONTINUE
C 
      DO 210 I = 1,N
        X(N,I) = B(INDX(N),I)/A(INDX(N),N)
        DO 200 J = N-1,1,-1
          X(J,I) = B(INDX(J),I)
          DO 190 K = J+1,N
            X(J,I) = X(J,I)-A(INDX(J),K)*X(K,I)
  190     CONTINUE
          X(J,I) =  X(J,I)/A(INDX(J),J)
  200   CONTINUE
  210 CONTINUE
C
C Restitution of A
C 
      DO 230 I = 1,N
        DO 220 J = 1,N
          A(I,J) = STOA(I,J)
  220   CONTINUE
  230 CONTINUE
C 
      RETURN
      END
C     _______________________________________________________________
C     ----------------------------------------------------------------
      SUBROUTINE ELGS(A,N,INDX)
C
C     Subroutine to perform the partial-pivoting Gaussian elimination
C     A(N,N) is the original matrix in the input and transformed
C     matrix plus the pivoting element ratios below the diagonal in
C     the output.  INDX(N) records the pivoting order.
C 
      IMPLICIT NONE
C 
      DOUBLE PRECISION A(N,N),C(N)
      INTEGER N,INDX(N)
      DOUBLE PRECISION C1,PI1,PI,PJ
      INTEGER K,ITMP,I,J
C
C     Initialize the index
C 
      DO 240 I = 1,N
        INDX(I) = I
  240 CONTINUE
C
C     Find the rescaling factors, one from each row
C 
        DO 260 I = 1,N
          C1 =  0.0
          DO 250 J = 1,N
            C1 = DMAX1(C1,DABS(A(I,J)))
  250     CONTINUE
          C(I) = C1
  260   CONTINUE
C
C     Search the pivoting (largest) element from each column
C 
      DO 300 J = 1,N-1
        PI1 = 0.0
        DO 270 I = J,N
          PI = DABS(A(INDX(I),J))/C(INDX(I))
          IF (PI.GT.PI1) THEN
            PI1 = PI
            K = I
          ELSE
          END IF
  270   CONTINUE
C
C     Interchange the rows via INDX(N) to record pivoting order
C 
        ITMP = INDX(J)
        INDX(J) = INDX(K)
        INDX(K) = ITMP
        DO 290 I = J+1,N
          PJ = A(INDX(I),J)/A(INDX(J),J)
C
C     Record pivoting ratios below the diagonal
C 
          A(INDX(I),J) = PJ
C
C     Modify other elements accordingly
C 
          DO 280 K = J+1,N
            A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
  280     CONTINUE
  290   CONTINUE
  300 CONTINUE
C 
      RETURN
      END
C
C     ----------------------------------------------------------------
      SUBROUTINE Jacobi(A,N,D,V,NROT)
      IMPLICIT NONE
      integer N,NROT, ip, iq, j, i, ialloc
      double precision A(1:N,1:N),D(1:N),V(1:N,1:N)
      double precision, pointer :: B(:), Z(:)
      double precision c,g,h,s,sm,t,tau,theta,tresh
      allocate(B(1:100),stat=ialloc)
      allocate(Z(1:100),stat=ialloc)

	do ip=1, N    !initialize V to identity matrix
	  do iq=1, N
	    V(ip,iq)=0.d0
	  end do
	    V(ip,ip)=1.d0
	end do
	do ip=1, N
	  B(ip)=A(ip,ip)
	  D(ip)=B(ip)
	  Z(ip)=0.d0
	end do
	NROT=0
	do i=1, 50
	  sm=0.d0
	  do ip=1, N-1     !sum off-diagonal elements
	    do iq=ip+1, N
	      sm=sm+DABS(A(ip,iq))
	    end do
	  end do
	  if(sm==0.d0) return  !normal return
	  if(i.lt.4) then
	    tresh=0.2d0*sm**2
	  else
	    tresh=0.d0
	  end if
	  do ip=1, N-1
	    do iq=ip+1, N
	      g=100.d0*ABS(A(ip,iq))
      ! after 4 sweeps, skip the rotation if the off-diagonal element is small
	      if((i.gt.4).and.(ABS(D(ip))+g.eq.ABS(D(ip))).and.
     &        (ABS(D(iq))+g.eq.ABS(D(iq)))) then
		A(ip,iq)=0.d0
	      else if(DABS(A(ip,iq)).gt.tresh) then
		h=D(iq)-D(ip)
		if(DABS(h)+g.eq.DABS(h)) then
		  t=A(ip,iq)/h
		else
		  theta=0.5d0*h/A(ip,iq)
		  t=1.d0/(DABS(theta)+DSQRT(1.d0+theta**2))
		  if(theta.lt.0.d0) t=-t
		end if
		c=1.d0/DSQRT(1.d0+t**2)
		s=t*c
		tau=s/(1.d0+c)
		h=t*A(ip,iq)
		Z(ip)=Z(ip)-h
		Z(iq)=Z(iq)+h
		D(ip)=D(ip)-h
		D(iq)=D(iq)+h
		A(ip,iq)=0.d0
		do j=1, ip-1
		  g=A(j,ip)
		  h=A(j,iq)
		  A(j,ip)=g-s*(h+g*tau)
		  A(j,iq)=h+s*(g-h*tau)
		end do
		do j=ip+1, iq-1
		  g=A(ip,j)
		  h=A(j,iq)
		  A(ip,j)=g-s*(h+g*tau)
		  A(j,iq)=h+s*(g-h*tau)
		end do
		do j=iq+1, N
		  g=A(ip,j)
		  h=A(iq,j)
		  A(ip,j)=g-s*(h+g*tau)
		  A(iq,j)=h+s*(g-h*tau)
		end do
		do j=1, N
		  g=V(j,ip)
		  h=V(j,iq)
		  V(j,ip)=g-s*(h+g*tau)
		  V(j,iq)=h+s*(g-h*tau)
		end do
		NROT=NROT+1
	      end if !if ((i.gt.4)...
	    end do !main iq loop
	  end do !main ip loop
	  do ip=1, N
	    B(ip)=B(ip)+Z(ip)
	    D(ip)=B(ip)
	    Z(ip)=0.d0
	  end do
	end do !main i loop
!	pause ' 50 iterations !'
	return
      END
C      
      subroutine printmatrix6(AAA)
          implicit none
          double precision :: AAA(6,6)
          PRINT *, ''
          PRINT *, AAA(1,1), AAA(1,2), AAA(1,3), AAA(1,4), AAA(1,5), AAA(1,6)
          PRINT *, AAA(2,1), AAA(2,2), AAA(2,3), AAA(2,4), AAA(2,5), AAA(2,6)
          PRINT *, AAA(3,1), AAA(3,2), AAA(3,3), AAA(3,4), AAA(3,5), AAA(3,6)
          PRINT *, AAA(4,1), AAA(4,2), AAA(4,3), AAA(4,4), AAA(4,5), AAA(4,6)
          PRINT *, AAA(5,1), AAA(5,2), AAA(5,3), AAA(5,4), AAA(5,5), AAA(5,6)
          PRINT *, AAA(6,1), AAA(6,2), AAA(6,3), AAA(6,4), AAA(6,5), AAA(6,6)
          PRINT *, ''
      end
