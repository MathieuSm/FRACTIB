C================================================================================================
C
C                          UMAT FOR ISOTROPIC ELASTIC CONSTITUTIVE LAW
C
C     FOUND ON https://simplifiedfem.wordpress.com/about/tutorial-write-a-simple-umat-in-abaqus/
C
C                  ADAPTED BY MATHIEU SIMON - ISTB - UNIVERSITY OF BERN - 2021
C
C================================================================================================
C
C
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
C
      INCLUDE 'ABA_PARAM.INC'
C
C
      CHARACTER*80 CMNAME
      DIMENSION DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     1 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     2 PROPS(NPROPS),COORDS(3),DROT(3,3), JSTEP(4)
C     
C     
      INTEGER NSTATV, KSTEP, KINC
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV)
      DOUBLE PRECISION DFGRD0(3,3),DFGRD1(3,3)
C
C
C     ELASTIC USER SUBROUTINE
      PARAMETER (ONE=1.0D0, TWO=2.0D0)
C
C
C    Get engineering variables
        E=PROPS(1)
        ANU=PROPS(2)
C
C
C    Compute lamé parameters
        ALAMBDA=E/(ONE+ANU)/(ONE-TWO*ANU)
        BLAMBDA=(ONE-ANU)
        CLAMBDA=(ONE-TWO*ANU)
C
C
C    Initialize material stiffness matrix
            DO I=1,NTENS
             DO J=1,NTENS
             DDSDDE(I,J)=0.0D0
             ENDDO
            ENDDO
C
C
C    Update stiffness matrx
                DDSDDE(1,1)=(ALAMBDA*BLAMBDA)
                DDSDDE(2,2)=(ALAMBDA*BLAMBDA)
                DDSDDE(3,3)=(ALAMBDA*BLAMBDA)
                DDSDDE(4,4)=(ALAMBDA*CLAMBDA)
                DDSDDE(5,5)=(ALAMBDA*CLAMBDA)
                DDSDDE(6,6)=(ALAMBDA*CLAMBDA)
                DDSDDE(1,2)=(ALAMBDA*ANU)
                DDSDDE(1,3)=(ALAMBDA*ANU)
                DDSDDE(2,3)=(ALAMBDA*ANU)
                DDSDDE(2,1)=(ALAMBDA*ANU)
                DDSDDE(3,1)=(ALAMBDA*ANU)
                DDSDDE(3,2)=(ALAMBDA*ANU)
C
C
C     Update stress tensor
         DO I=1,NTENS
            DO J=1,NTENS
            STRESS(I)=STRESS(I)+DDSDDE(I,J)*DSTRAN(J)
            ENDDO
	     ENDDO
C
C        
C     Deformation gradient
      STATEV(1) = DFGRD1(1,1)
      STATEV(2) = DFGRD1(1,2)
      STATEV(3) = DFGRD1(1,3)
      STATEV(4) = DFGRD1(2,1)
      STATEV(5) = DFGRD1(2,2)
      STATEV(6) = DFGRD1(2,3)
      STATEV(7) = DFGRD1(3,1)
      STATEV(8) = DFGRD1(3,2)
      STATEV(9) = DFGRD1(3,3)
C
C
        RETURN
        END
