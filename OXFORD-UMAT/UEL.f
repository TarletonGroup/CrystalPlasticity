C     FIRST WORKING VERSION
C     FEB. 08TH, 2024
C     ERALP DEMIR
C     UNIVERSITY OF OXFORD
C     ENGINEERING SCIENCE DEPARTMENT
C     UEL - FOR COHESIVE ELEMENTS
C
C     FRACTURE MODES: MODE-I / MIXED MODE
C     FOR ELEMENT TYPES: 
C     COH2D2 / COH2D3 
C     COH3D6 / COH3D8 / COH3D12 / COH3D16
C     COMPATIBLE WITH THE ELEMENTS
C     (INCLUDING THEIR REDUCED INTEGRATION FORMS):
C     CPS3 / CPE3 / CPS4 / CPE4 / CPS6 / CPE6 / CPS8 / CPE8
C     C3D4 / C3D6 / C3D8 / C3D10 / C3D15 / C3D20
C
C     NOTES:
C     1) CONNECTIVITY OF THE COHESIVE ELEMENTS 
C     MUST BE SUCH THAT THE BOTTOM NODES FOLLOWED
C     BY THE TOP NODES NORMAL TO THE COHESIVE PLANE OR
C     THE PLANE AT WHICH CRACK OPENING (MODE-I) OCCURS
C     2) UMAT IS USED HERE FOR ONLY VISUALIZATION PURPOSE
C     3) STATE VARIABLES OF UEL INCLUDE: 
C         - DAMAGE VARIABLE AT THE IPS.
C     
C
C     INPUTS 
C     *: REQUIRED ONLY FOR MIXED MODE; "0" OTHERWISE.
C     `: REQURED FOR 2D PROBLEMS ONLY 
C     PROPS(1) / MODE / FRACTURE MODE / (-): 
C     = "1" FOR OPENING MODE ONLY
C     = "2" FOR OPENING MODE ONLY
C     = "3" FOR MIXED MODE 
C     PROPS(2) / KI / PENALTY STIFFNESS FOR OPENING MODE / (N/MM)
C     PROPS(3) / KII / PENALTY STIFFNESS FOR SHEARING MODE / (N/MM)
C     PROPS(4) / SI / MAXIMUM NORMAL STRESS BEFORE FAILURE / (MPA)
C    *PROPS(5) / SII / MAXIMUM SHEAR STRESS BEFORE FAILURE / (MPA)
C     PROPS(6) / GCI / FRACTURE TOUGHNESS FOR OPENING MODE / (J)
C    *PROPS(7) / GCII / FRACTURE TOUGHNESS FOR SHEARING MODE / (J)
C    *PROPS(8) / ETA / BENZEGGAGH-KENANE (B-K) EXPONENT FOR MIXED MODE / (-)
C    `PROPS(9) / HEIGHT / HEIGHT OF 2D ELEMENTS / (MM) 
C      
C      
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELC,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
C      INCLUDE 'ABA_PARAM.INC'
C
      IMPLICIT NONE
C
C     USED VARIABLES
      INTEGER, INTENT(IN) :: NDOFEL,NRHS,NSVARS,NPROPS,
     1 MCRD,NNODE,JTYPE,KSTEP,KINC,MLVARX,NJPROP,
     2 JELC
      REAL(8), INTENT(IN) :: PROPS(NPROPS),
     1 COORDS(MCRD,NNODE),U(NDOFEL),
     2 PNEWDT,DTIME
      REAL(8), INTENT(INOUT) :: SVARS(NSVARS),RHS(MLVARX,NRHS)
      REAL(8), INTENT(OUT) :: AMATRX(NDOFEL,NDOFEL)
C     UNUSED VARIABLES
      INTEGER, INTENT(IN) :: LFLAGS(7),JPROPS(NJPROP),
     1 MDLOAD,NDLOAD,NPREDF
      REAL(8), INTENT(IN) :: DU(NDOFEL,NRHS),
     1 V(NDOFEL),A(NDOFEL),
     2 TIME(2),PARAMS(*),JDLTYP(MDLOAD,*),
     3 ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),PERIOD
      REAL(8), INTENT(OUT) :: ENERGY(8)
C
C
C     DEFAULT NUMERICAL PARAMETER(S)
C     --------------------------------------------
C
C     NUMERICAL INTEGRATION METHOD
C     "1": GAUSS-QUADRATURE
C     "2": NEWTON-COTES
C     "3": NEWTON-COTES CONSISTENT WITH ABAQUS (CORNER NODES ONLY)
C     DEFAULT SETTING IS "1" FOR GAUSS-QUADRATURE
      INTEGER, PARAMETER :: INTMTD = 1
C
C     NUMBER OF STATE VARIABLES PER IP
C     1. DAMAGE
      INTEGER, PARAMETER :: NSVPN = 1
C     --------------------------------------------
C
C
C
C
C     VARIABLES USED
C     NUMBER OF INTEGRATION POINTS (IP) PER COHESIVE ELEMENT
      INTEGER :: NSQPT
C     NUMBER OF NODES PER SURFACE
      INTEGER :: NSNODE
C
	INTEGER :: IP, I, J, K
C     INTEGRATION POINT COORDINATES IN ISOPARAMETRIC SPACE
      REAL(8), DIMENSION(INT(NSVARS/NSVPN),MCRD-1) :: IPCOORDS
C     INTEGRATION POINT WEIGHTS
      REAL(8), DIMENSION(INT(NSVARS/NSVPN)) :: WT
C
C     MODEL PARAMETERS
C     FAILURE MODES
C     "1": OPENING MODE (MODE I) ONLY
C     "2": MIXED MODE
C     DEFAULT SETTING IS "1" FOR MODE I ONLY
      INTEGER :: MODE
C     NORMAL STIFFNESS (MPA)
      REAL(8) :: KI
C     SHEAR STIFFNESS (MPA)
      REAL(8) :: KII
C     MAXIMUM STRESS BEFORE FAILURE (MPA)
C     NORMAL STRENGTH
      REAL(8) :: SI
C     SHEAR STRENGTH
      REAL(8) :: SII
C     FRACTURE TOUGHNESS FOR OPENING (J)
      REAL(8) :: GCI
C     FRACTURE TOUGHNESS FOR SHEAR (J)
      REAL(8) :: GCII
C     B-K EXPONENT FOR MIXED MODE (-)
C     IF ETA<0, QUADRATIC FAILURE LAW
      REAL(8) :: ETA
C     HEIGHT OF 2D ELEMENTS (MM)
      REAL(8) :: HEIGHT
C
C
C
C     DEBUGGER
      INTEGER :: DEBUG, DEBUGWAIT
C
C     RESET UNUSED OUTPUTS
      ENERGY=0.
     
C
C
C     DEBUGGER
      DEBUG=0
      DO WHILE (DEBUG.EQ.1)
          DEBUGWAIT=1
      END DO
C
C
C     READ PROPS
      CALL READPROPS(NPROPS,PROPS,
     + MODE,KI,KII,SI,SII,GCI,GCII,ETA,HEIGHT)
C
C     DETERMINE THE ELEMENT PROPERTIES
      CALL ELPROP(INTMTD,NDOFEL,MCRD,
     + NSQPT,NSNODE)
C
C
C      WRITE(*,*) 'UEL'
C      WRITE(*,*) 'JELC',JELC
C      WRITE(*,*) 'KINC', KINC
C
C
C
C     CALCULATE INTEGRATION POINT COORDINATES AND WEIGHTS
      CALL INTEGRATIONPOINTS(INTMTD,NDOFEL,NSQPT,MCRD,
     + IPCOORDS,WT)
C
C
C     AT THE VERY FIRST STEP
C     OUTPUT IMPORTANT ENTRIES
      IF ((KINC.EQ.1).AND.(KSTEP.EQ.1)) THEN
C
C
C
C         RESET SVD (STRESS, STRAIN, DAMAGE)
          K=0
          DO I =1,NSVARS
              SVARS(K)=0.
          END DO
C
C         WRITE/OUTPUT AND INITIAL CHECK  
          IF (JELC.EQ.1) THEN
C
              IF (NSVARS.LT.(NSQPT*NSVPN)) THEN
                  WRITE(*,*) 'INCREASE THE NUMBER OF SDV TO ', 
     + NSQPT*NSVPN
                  CALL XIT
              END IF
C
C
              WRITE(*,*) 'NDOFEL',NDOFEL
              WRITE(*,*) 'MLVARX',MLVARX
              WRITE(*,*) 'NNODE',NNODE
              WRITE(*,*) 'NSVARS',NSVARS
              WRITE(*,*) 'NPROPS',NPROPS
              WRITE(*,*) 'MCRD',MCRD
              WRITE(*,*) 'NRHS',NRHS
C
C
              DO I =1,NPROPS
                  WRITE(*,*) 'PROPS',I, PROPS(I)
              END DO
C
C
              DO I=1,NNODE
                  WRITE(*,*) 'COORDS',I, COORDS(1:MCRD,I)
              END DO
C
              DO I=1,NSVARS
                  WRITE(*,*) 'SVARS',I, SVARS(I)
              END DO
C
              WRITE(*,*) '*********************'
C
          END IF
C
C
      END IF     
C
C
C
C
C
CC     RESET RHS
C	DO I=1,MLVARX
C	    RHS(I,1)=0.
C      END DO
C
C
C      DO I=1,NNODE
C          WRITE(*,*) 'COORDS',I, COORDS(1:MCRD,I)
C      END DO
C
C
C      WRITE(*,*) 'U'
C      DO I=1,NDOFEL
C          WRITE(*,*) U(I)
C      END DO
C
C
C
C      WRITE(*,*) 'SVARS'
C      WRITE(*,*) SVARS
C
C
      CALL AMATRXANDRHS(MCRD,
     + NNODE,NSVARS,NSVPN,
     + NSNODE,NSQPT,NDOFEL,
     + MLVARX,NRHS,IPCOORDS,WT,COORDS,U,
     + MODE,KI,KII,SI,SII,GCI,GCII,
     + ETA,HEIGHT,
     + AMATRX,RHS,
     + SVARS)
C
C
C
C
C
C
C      WRITE(*,*) 'SVARS'
C      WRITE(*,*) SVARS
C	
C
C     
C     
      RETURN
      END SUBROUTINE UEL
C
C
C
      SUBROUTINE READPROPS(NPROPS,PROPS,
     + MODE,KI,KII,SI,SII,GCI,GCII,ETA,HEIGHT)
      IMPLICIT NONE
C     LENGTH OF PROPS VECTOR
      INTEGER, INTENT(IN) :: NPROPS
C     PROPS VECTOR
      REAL(8), DIMENSION(NPROPS), INTENT(IN) :: PROPS
C     OUTPUTS
      INTEGER, INTENT(OUT) :: MODE 
      REAL(8), INTENT(OUT) :: KI
      REAL(8), INTENT(OUT) :: KII
      REAL(8), INTENT(OUT) :: SI
      REAL(8), INTENT(OUT) :: SII
      REAL(8), INTENT(OUT) :: GCI
      REAL(8), INTENT(OUT) :: GCII
      REAL(8), INTENT(OUT) :: ETA
      REAL(8), INTENT(OUT) :: HEIGHT
C
C
C     FAILURE MODE ("1": MODE I ONLY / "2": MIXED MODE)
      MODE=INT(PROPS(1))
C     NORMAL PENALTY STIFFNESS (N/MM)
      KI = PROPS(2)
C     SHEAR PENALTY STIFFNESS (N/MM)
      KII = PROPS(3)
C     CHECK FOR KII (MAY BE SKIPPED BY THE USER FOR MIXED MODE ENTRY
C     BUT, IT IS REQUIRED FOR MODE I TOO.)
      IF (KII.EQ.0.) KII=KI
C     INTERFACE NORMAL STRENGTH (MPA)
      SI = PROPS(4)
C     INTERFACE SHEAR STRENGTH (MPA)
      SII = PROPS(5)
C     FRACTURE TOUGHNESS FOR OPENING (J)
      GCI = PROPS(6)
C     FRACTURE TOUGHNESS FOR OPENING (J)
      GCII = PROPS(7)
C     B-K EXPONENT FOR MIXED MODE (-)
      ETA = PROPS(8)
C     HEIGHT OF 2D ELEMETS (MM)
      HEIGHT = PROPS(9)
C     SET TO UNITY IF FORGOTTON
      IF (HEIGHT.EQ.0.) HEIGHT=1.
C
      RETURN
      END SUBROUTINE READPROPS
C
C
C
C
      SUBROUTINE AMATRXANDRHS(MCRD,
     + NNODE,NSVARS,NSVPN,NSNODE,
     + NSQPT,NDOFEL,MLVARX,NRHS,
     + GP,WT,COORDS,U,MODE,
     + KI,KII,SI,SII,GCI,GCII,
     + ETA,HEIGHT,
     + AMATRX,RHS,
     + SVARS)
      IMPLICIT NONE
C     INPUTS
      INTEGER, INTENT(IN) :: MCRD
      INTEGER, INTENT(IN) :: NNODE
      INTEGER, INTENT(IN) :: NSVARS
      INTEGER, INTENT(IN) :: NSVPN
      INTEGER, INTENT(IN) :: NSNODE
      INTEGER, INTENT(IN) :: NSQPT
      INTEGER, INTENT(IN) :: NDOFEL
      INTEGER, INTENT(IN) :: MLVARX
      INTEGER, INTENT(IN) :: NRHS
      REAL(8), DIMENSION(NSQPT,MCRD-1), INTENT(IN) :: GP
      REAL(8), DIMENSION(NSQPT), INTENT(IN) :: WT
      REAL(8), DIMENSION(MCRD,NNODE), INTENT(IN) :: COORDS
      REAL(8), DIMENSION(NDOFEL), INTENT(IN) :: U
C     MATERIAL PARAMETERS
      INTEGER, INTENT(IN) :: MODE
      REAL(8), INTENT(IN) :: KI
      REAL(8), INTENT(IN) :: KII
      REAL(8), INTENT(IN) :: SI
      REAL(8), INTENT(IN) :: SII
      REAL(8), INTENT(IN) :: GCI
      REAL(8), INTENT(IN) :: GCII
      REAL(8), INTENT(IN) :: ETA
      REAL(8), INTENT(IN) :: HEIGHT
C
C     OUTPUTS
C     STIFFNESS
      REAL(8), DIMENSION(NDOFEL,NDOFEL), INTENT(OUT) :: AMATRX
      REAL(8), DIMENSION(MLVARX,NRHS), INTENT(INOUT) :: RHS
      REAL(8), DIMENSION(NSVARS), INTENT(INOUT) :: SVARS
C
C
      REAL(8) :: I3(3)
      REAL(8) :: N(NSNODE), DN(MCRD-1,NSNODE)
      REAL(8) :: B(MCRD,NDOFEL), FINT(NDOFEL)
      REAL(8) :: DEFCOORDS(MCRD,NNODE)
      REAL(8) :: G, H, AREA, DAMAGE, DAMAGE0
      REAL(8) :: DNDA(NSNODE), DNDB(NSNODE)
      REAL(8) :: R(NNODE*MCRD,NNODE*MCRD), NVEC(MCRD)
      REAL(8) :: TJAC(MCRD,MCRD), RJS, LOCDEFCOORDS(MCRD,NNODE)
      REAL(8) :: DELTA(MCRD,NSNODE)
      REAL(8) :: DELTAU(MCRD)
      REAL(8) :: T(MCRD), DTDDELTA(MCRD,MCRD)
      REAL(8) :: LOCFEL(NDOFEL), FEL(NDOFEL)
      REAL(8) :: LOCKEL(NDOFEL,NDOFEL), KEL(NDOFEL,NDOFEL)
      REAL(8), DIMENSION(NSQPT) :: DAMAGEIP, DAMAGEIP0
      INTEGER :: I, J, K, J1, J2, J3, IP
C
C
C     READ VARIABLES
      DO IP = 1, NSQPT
C         1. DAMAGE VARIABLE 
          J1 = (IP-1)*NSVPN+1
          DAMAGEIP0(IP)=SVARS(J1)
      END DO
C
C
C     DOES NOT WORK WELL - USE UNDEFORMED
C     DEFORMED COORDINATES
C     BASED ON THE DIMENSION OF THE PROBLEM
      IF (MCRD.EQ.2) THEN
C         2D
          DO I=1,NNODE
              J2=MCRD*I
              J1=J2-1
              DEFCOORDS(2,I)=COORDS(2,I) + U(J2)
              DEFCOORDS(1,I)=COORDS(1,I) + U(J1)
          END DO      
      ELSEIF (MCRD.EQ.3) THEN
C         3D
          DO I=1,NNODE
              J3=MCRD*I
              J2=J3-1
              J1=J2-1
              DEFCOORDS(3,I)=COORDS(3,I) + U(J3)
              DEFCOORDS(2,I)=COORDS(2,I) + U(J2)
              DEFCOORDS(1,I)=COORDS(1,I) + U(J1)
          END DO
      END IF
C
C
C
C
C     ELEMENT AREA
      AREA=0.
C
      AMATRX = 0.
      DO I=1,MLVARX
          RHS(I,1) = 0.
      END DO
C
C
C
C     STRAIN-DISPLACEMENT
      DO IP=1,NSQPT
C
C
C
          DAMAGE0 = DAMAGEIP0(IP)
C
C          WRITE(*,*) 'IP', IP
C
C
C         BASED ON THE TYPE OF THE COHESIVE ELEMENT
C         ELEMENT TYPE SPECIFIC PROPERTIES
          SELECT CASE(NDOFEL)
C         C2D4
          CASE(8)
C             QUAD POINT COORDINATES (1D LINE FOR 2D COHESIVE ELEMENT)
              G=GP(IP,1)
C             SHAPE FUNCTIONS AND ITS DERIVATIVES OF
C             LINE TRACTION FOR 2D PROBLEMS
              CALL COH2D4_N_DN(NSNODE,MCRD-1,G,N,DN)
C         C2D6
          CASE(12)
C             QUAD POINT COORDINATES (1D LINE FOR 2D COHESIVE ELEMENT)
              G=GP(IP,1)
C             SHAPE FUNCTIONS AND ITS DERIVATIVES OF
C             LINE TRACTION FOR 2D PROBLEMS
              CALL COH2D6_N_DN(NSNODE,MCRD-1,G,N,DN)
C         C3D6
          CASE(18)
C             QUAD POINT COORDINATES (2D SURFACE FOR 3D COHESIVE ELEMENT)
              G=GP(IP,1)
              H=GP(IP,2)
C             SHAPE FUNCTIONS AND ITS DERIVATIVES OF
C             SURFACE TRACTION FOR 3D
              CALL COH3D6_N_DN(NSNODE,MCRD-1,G,H,N,DN)
C         C3D8
          CASE(24)
C             QUAD POINT COORDINATES (2D SURFACE FOR 3D COHESIVE ELEMENT)
              G=GP(IP,1)
              H=GP(IP,2)
C             SHAPE FUNCTIONS AND ITS DERIVATIVES OF
C             SURFACE TRACTION FOR 3D
              CALL COH3D8_N_DN(NSNODE,MCRD-1,G,H,N,DN)
C         C3D12
          CASE(36)
C             QUAD POINT COORDINATES (2D SURFACE FOR 3D COHESIVE ELEMENT)
              G=GP(IP,1)
              H=GP(IP,2)
C             SHAPE FUNCTIONS AND ITS DERIVATIVES OF
C             SURFACE TRACTION FOR 3D
              CALL COH3D12_N_DN(NSNODE,MCRD-1,G,H,N,DN)
C         C3D16
          CASE(48)
C             QUAD POINT COORDINATES (2D SURFACE FOR 3D COHESIVE ELEMENT)
              G=GP(IP,1)
              H=GP(IP,2)
C             SHAPE FUNCTIONS AND ITS DERIVATIVES OF
C             SURFACE TRACTION FOR 3D
              CALL COH3D16_N_DN(NSNODE,MCRD-1,G,H,N,DN)
          END SELECT
C
C
C
C         MID SURFACE DISPLACCENT - NODE DISPLACMENT COUPLING
          CALL BMATRX(NDOFEL,MCRD,NNODE,NSNODE,N,B)
C
C
C
C
C
C
C         TRANSFORMATION MATRIX
          CALL RMATRX(MCRD,NNODE,
     + NSNODE,DEFCOORDS,HEIGHT,DN,
     + R,TJAC,NVEC,RJS,LOCDEFCOORDS)
C
          IF (RJS.LT.0.) THEN
              WRITE(*,*) 'COHESIVE ELEMENT CONNECTIVITY IS WRONG!'
              WRITE(*,*) 'EXITING!'
              CALL XIT
          END IF
C
C
C          WRITE(*,*) 'NVEC: ', NVEC
C          WRITE(*,*) 'RJS', RJS
C
C          WRITE(*,*) 'TJAC: '
C          DO I=1,MCRD
C              WRITE(*,*) (TJAC(I,J),J=1,MCRD)
C          END DO
C
C
C
C         MID-SURFACE OPENING DISPLACEMENTS
          DO I=1,NSNODE
C
              DO J=1,MCRD
                  DELTA(J,I) = 
     + LOCDEFCOORDS(J,I+NSNODE)
     + - LOCDEFCOORDS(J,I)
              END DO
C
          END DO
C
C         OPENING DISPLACEMENTS AT THE INTEGRATION POINTS
          DELTAU=0.
          DO J=1,MCRD
              DO I=1,NSNODE
                  DELTAU(J) = DELTAU(J) + N(I)*DELTA(J,I)
              END DO
          END DO
C
C
C          WRITE(*,*) 'DELTAU'
C          WRITE(*,*) DELTAU
C
C         COHESIVE LAW
C
C         MODE I ONLY
          IF (MODE.EQ.1) THEN
C
              CALL COHESIVELAW_OPENINGMODE(
     + MCRD,DELTAU,DAMAGE0,KI,KII,SI,GCI,
     + DAMAGE,T,DTDDELTA)
C
C         MODE II ONLY
          ELSEIF (MODE.EQ.2) THEN
C
              CALL COHESIVELAW_SHEARMODE(
     + MCRD,DELTAU,DAMAGE0,KI,KII,SII,GCII,
     + DAMAGE,T,DTDDELTA)
C
C
C         MIXED MODE
          ELSEIF (MODE.EQ.3) THEN
C
              CALL COHESIVELAW_MIXEDMODE(
     + MCRD,DELTAU,DAMAGE0,KI,KII,SI,SII,
     + GCI,GCII,ETA,DAMAGE,T,DTDDELTA)
C
          ENDIF
C
C          WRITE(*,*) 'T'
C          WRITE(*,*) T
C
C          WRITE(*,*) 'DAMAGE-BEFORE'
C          WRITE(*,*) DAMAGE0
C
          DAMAGEIP(IP)=DAMAGE
C
C          WRITE(*,*) 'DAMAGE-AFTER'
C          WRITE(*,*) DAMAGE
C
C
C          WRITE(*,*) 'DTDDELTA'
C          DO I=1,MCRD
C              WRITE(*,*) (DTDDELTA(I,J),J=1,MCRD)
C          END DO
C
C         ELEMENT FORCE
          LOCFEL = MATMUL(TRANSPOSE(B),T)
C
C          WRITE(*,*) 'B'
C          DO I=1,NDOFEL
C              WRITE(*,*) (B(J,I),J=1,MCRD)
C          END DO
C
C
C
C         TRANSFORM THE ELEMENT FORCE TO GLOBAL
          FEL = MATMUL(TRANSPOSE(R),LOCFEL)
C
C
C          WRITE(*,*) 'FEL'
C          WRITE(*,*) FEL
C
C         CACULATE NODAL FORCES
          FINT = WT(IP)*RJS*FEL
C
C
C
C         INTERNAL FORCE
C
C         BASED ON THE DIMENSION OF THE PROBLEM
          IF (MCRD.EQ.2) THEN
C             2D
              DO J=1,NNODE
                  J2=MCRD*J
                  J1=J2-1
                  RHS(J1,1) = RHS(J1,1) - FINT(J1)
                  RHS(J2,1) = RHS(J2,1) - FINT(J2)
              END DO
          ELSEIF (MCRD.EQ.3) THEN
C             3D
              DO J=1,NNODE
                  J3=MCRD*J
                  J2=J3-1
                  J1=J2-1
                  RHS(J1,1) = RHS(J1,1) - FINT(J1)
                  RHS(J2,1) = RHS(J2,1) - FINT(J2)
                  RHS(J3,1) = RHS(J3,1) - FINT(J3)
              END DO
          END IF
C
C
C         CALCULATE LOCAL ELEMENT STIFFNESS
          LOCKEL=MATMUL(TRANSPOSE(B),MATMUL(DTDDELTA,B))
C
C
C
C         TRANSFORM THE STIFFNESS TO THE GLOBAL REFERENCE
          KEL = MATMUL(TRANSPOSE(R),MATMUL(LOCKEL,R))
C
C
C
C
C
C         STIFFNESS
          AMATRX = AMATRX + WT(IP)*RJS*KEL
C
C
C
C          WRITE(*,*) '*****************************'
C
C
C
C
C
      END DO
C
C
C     SURFACE AREA OF THE MID-SURFACE
      AREA = AREA + RJS
C
C
C
C
C
C      WRITE(*,*) 'AMATRX'
C      DO I=1,NDOFEL
C          WRITE(*,*) (AMATRX(I,J), J=1,NDOFEL)
C      END DO
C
C
C      WRITE(*,*) 'RHS'
C      DO I=1,NDOFEL
C          WRITE(*,*) RHS(I,1)
C      END DO
C
C
C
C     OVERWRITE VARIABLES
      DO IP = 1, NSQPT
C         1. DAMAGE VARIABLE 
          J1 = (IP-1)*NSVPN+1
          IF (DAMAGEIP(IP).GT.DAMAGEIP0(IP)) THEN
              SVARS(J1)=DAMAGEIP(IP)
          END IF
      END DO
C
C
      RETURN
      END SUBROUTINE AMATRXANDRHS
C
C
C
C     TRANSFORMATION MATRIX FROM THE GLOBAL ELEMENT
C     TO THE LOCAL MID-SURFACE/LINE REFERENCE
      SUBROUTINE RMATRX(MCRD,NNODE,NSNODE,DEFCOORDS,HEIGHT,DN,
     + R,T,NVEC,RJS,LOCDEFCOORDS)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MCRD
      INTEGER, INTENT(IN) :: NNODE
      INTEGER, INTENT(IN) :: NSNODE
      REAL(8), DIMENSION(MCRD,NNODE), INTENT(IN) :: DEFCOORDS
      REAL(8), INTENT(IN) :: HEIGHT
      REAL(8), DIMENSION(MCRD-1,NSNODE), INTENT(IN) :: DN
C     OVERALL TRANSFORMATION MATRIX
      REAL(8), DIMENSION(NNODE*MCRD,NNODE*MCRD), INTENT(OUT) :: R
C     JACOBIAN
      REAL(8), DIMENSION(MCRD,MCRD), INTENT(OUT) :: T
C     UNIT MID-PLANE/LINE NORMAL
      REAL(8), DIMENSION(MCRD), INTENT(OUT) :: NVEC
C     SURFACE OR LINE JACOBIAN (AREA OR LENGTH)
      REAL(8), INTENT(OUT) :: RJS
      REAL(8), DIMENSION(MCRD,NNODE), INTENT(OUT) :: LOCDEFCOORDS
C
C     VARIABLES WITHIN
      REAL(8) :: DNDA(NSNODE), DNDB(NSNODE)
      REAL(8) :: SVEC(MCRD), TVEC(MCRD)
      REAL(8) :: DET, TINV(MCRD,MCRD)
      REAL(8) :: SDEFCOORDS(MCRD,NSNODE)
      INTEGER :: I, J, K, IND
C
      SDEFCOORDS=0.
      DO I=1,MCRD
C         AVERAGE THE TOP AND BOTTOM NODE COORDINATES
          DO J=1,NSNODE
              SDEFCOORDS(I,J) = 
     + (DEFCOORDS(I,J)+DEFCOORDS(I,J+NSNODE))/2.
          END DO
      END DO
C
C      WRITE(*,*) 'SDEFCOORDS'
C      DO I=1,NSNODE
C          WRITE(*,*) (SDEFCOORDS(J,I),J=1,MCRD)
C      END DO
C
C
C     DIMENSIONS
C     3D
      IF (MCRD.EQ.3) THEN
          DO I=1,NSNODE
              DNDA(I) = DN(1,I)
              DNDB(I) = DN(2,I)
          END DO
C
          CALL SURFJAC(MCRD,NSNODE,DNDA,DNDB,SDEFCOORDS,
     + SVEC,TVEC,NVEC,RJS)
C
          DO I=1,MCRD
              T(1,I)=SVEC(I)
              T(2,I)=TVEC(I)
              T(3,I)=NVEC(I)
          END DO
C
!          CALL INV3X3(T,TINV,DET)
C
C     2D
      ELSEIF (MCRD.EQ.2) THEN
          DO I=1,NSNODE
              DNDA(I) = DN(1,I)
          END DO
C
          CALL LINEJAC(MCRD,NSNODE,
     + DNDA,SDEFCOORDS,HEIGHT,
     + SVEC,NVEC,RJS)
C
          DO I=1,MCRD
              T(1,I)=SVEC(I)
              T(2,I)=NVEC(I)
          END DO
C
!          CALL INV2X2(T,TINV,DET)
C
      END IF
C
C      WRITE(*,*) 'TJAC'
C      DO I=1,MCRD
C          WRITE(*,*) (T(I,J),J=1,MCRD)
C      END DO
C
C      WRITE(*,*) 'RJS', RJS
C
C     DEFORMED COORDINATES IN THE LOCAL REFERENCE
      LOCDEFCOORDS=0.
      DO I=1,NNODE
          DO J=1,MCRD
              DO K=1,MCRD
                  LOCDEFCOORDS(J,I) =
     + LOCDEFCOORDS(J,I) + T(J,K) * DEFCOORDS(K,I)
              END DO
          END DO
      END DO
C
C
C
C
      R = 0.
      DO I = 1, NNODE
C
          IND = MCRD * (I-1)
C
          DO J = 1, MCRD
              DO K = 1,MCRD
                  R(J+IND, K+IND) = T(J,K)
              END DO
          END DO
C
      END DO
C
C
C
C      WRITE(*,*) 'LOCDEFCOORDS'
C      DO I=1,NNODE
C          WRITE(*,*) (LOCDEFCOORDS(J,I),J=1,MCRD)
C      END DO
C
      RETURN
      END SUBROUTINE RMATRX
C
C
C
C
C     COHESIVE MODEL FOR ONLY-OPENING MODE FAILURE
      SUBROUTINE COHESIVELAW_OPENINGMODE(
     + MCRD,DELTA,DAMAGE0,KI,KII,SI,GCI,
     + DAMAGE,T,DTDDELTA)
C     COMMENT OUT THE IMPLICIT NONE
C     OTHERWISE MACAULEY FUNCTION DOES NOT WORK
C      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MCRD
      REAL(8), DIMENSION(MCRD), INTENT(IN) :: DELTA
      REAL(8), INTENT(IN) :: DAMAGE0
C     MATERIAL PARAMETERS
      REAL(8), INTENT(IN) :: KI
      REAL(8), INTENT(IN) :: KII
      REAL(8), INTENT(IN) :: SI
      REAL(8), INTENT(IN) :: GCI
      REAL(8), INTENT(OUT) :: DAMAGE
      REAL(8), DIMENSION(MCRD), INTENT(OUT) :: T
      REAL(8), DIMENSION(MCRD,MCRD), INTENT(OUT) :: DTDDELTA
C
C     VARIABLES WITHIN
      REAL(8) :: DI, DII, DELTA0, DELTAF
      REAL(8) :: DELTAMAX, LAMBDA
      INTEGER :: F
      REAL(8) :: D0(MCRD,MCRD), D(MCRD,MCRD), H
      REAL(8) :: B, A(MCRD,MCRD), C(MCRD)
      INTEGER :: I, J, K
C
C     MATERIAL CONSTANTS
C
C     MAXIMUM DISPLACEMENT BEFORE FAILURE
      DELTA0 = SI / KI
C
C      WRITE(*,*) 'DELTA0', DELTA0
C
C     DISPLACEMENT AT FAILURE
      DELTAF = 2. * GCI / SI
C     
C     MODE-I
      DI = DELTA(MCRD)
      LAMBDA = SQRT(DI**2)
C
C      WRITE(*,*) 'DELTA', DELTA
C
C     MODE-II
      DII=0.
      DO I=1,MCRD-1
          DII=DII+DELTA(I)**2
      END DO
      DII = SQRT(DII)
C
C     DAMAGE THRESHOLD
      DELTAMAX = DELTA0*DELTAF/(DELTAF-DAMAGE0*(DELTAF-DELTA0))
C      DELTAMAX = DELTA0 + DAMAGE0*(DELTAF-DELTA0)
C
C      WRITE(*,*) 'DAMAGE0', DAMAGE0
C
C      WRITE(*,*) 'DELTAMAX', DELTAMAX
C
C      WRITE(*,*) 'DELTAF', DELTAF
C
C      WRITE(*,*) 'DELTA0', DELTA0
C
C
C     IF THE DAMAGE THRESHOLD IS EXCEEDED
      IF (DI.GT.DELTAMAX) THEN
C
C         LOADING
          F = 1
C
          DELTAMAX = DI
C
      ELSE
C
C         UNLOADING
          F = 0
C
      END IF
C
C
C     UPDATE DAMAGE
      DAMAGE = DELTAF*(DELTAMAX-DELTA0)/DELTAMAX/(DELTAF-DELTA0)
C
C
C     CHECK THE DAMAGE PARAMETER
      IF (DAMAGE.GT.1.) THEN
C
          DAMAGE = 1.
          F = 0
C
      END IF
C
C
C     UNDAMAGED STIFFNESS
      D0 = 0.
      DO I=1,MCRD
          IF (I.EQ.MCRD) THEN
              D0(I,I)=KI
          ELSE
              D0(I,I)=KII
          END IF
      END DO
C
C
C     DAMAGED STIFFNESS MATRIX
      D = (1.-DAMAGE)*D0
C
C
C     CLOSURE - RETAINS BACK THE OPENING STIFFNESS
      IF (DI.LT.0.) THEN
          D(MCRD,MCRD) = KI
      END IF
C
C
C
      DTDDELTA=D
C
C
C     ADD THE REMAINING TERMS
      IF (F.EQ.1) THEN
C
C
          B = MACAULEY(-DI)/DI
C
          A = 0.
          DO I=1,MCRD
              A(I,I)=1.    
          END DO
          A(MCRD,MCRD) = A(MCRD,MCRD) + B
C
C
          C=0.
          DO I=1,MCRD
              DO J=1,MCRD
                  C(I)=C(I)+A(I,J)*DELTA(J)
              END DO
          END DO
C
C
C
C         TANGENT STIFFNESS CONSTANT
C          IF (LAMBDA.GT.0.) THEN
              H = DELTAF*DELTA0/(DELTAF-DELTA0)/LAMBDA**3
C          ELSE
C              H = 0.
C          END IF
C          
C
C
          DO I=1,MCRD
              DO J=1,MCRD
C                  DO K=1,MCRD
                      DTDDELTA(I,J)=
     + DTDDELTA(I,J) - H*KI*C(I)*C(J)
C                  END DO
              END DO
          END DO
C
C         PENETRATION
          IF (DELTA(MCRD).LT.0.) THEN
              DO I=1,MCRD
                  DTDDELTA(MCRD,I)=D(MCRD,I)
                  DTDDELTA(I,MCRD)=D(I,MCRD)
              END DO
          END IF
C
C
      END IF
C
C
C     CHECK THE TANGENT STIFFNESS MATRIX
      DO I=1,MCRD
          DO J=1,MCRD
              IF (ABS(DTDDELTA(I,J)).GT.D(I,J))
     + DTDDELTA(I,J) = D(I,J)
          END DO
      END DO
C
C
C     CALCULATE TRACTION
      T=0.
      DO I=1,MCRD
C          DO J=1,MCRD
C              T(I) = T(I) + D(I,J)*DELTA(J)
              T(I) = D(I,I)*DELTA(I)
C          END DO
      END DO
C
C
C      WRITE(*,*) 'T'
C      WRITE(*,*) T
C
C
C
      RETURN
      END SUBROUTINE COHESIVELAW_OPENINGMODE
C      
C
C
C
C
C     COHESIVE MODEL FOR ONLY-SHEAR MODE FAILURE
      SUBROUTINE COHESIVELAW_SHEARMODE(
     + MCRD,DELTA,DAMAGE0,KI,KII,SII,GCII,
     + DAMAGE,T,DTDDELTA)
C     COMMENT OUT THE IMPLICIT NONE
C     OTHERWISE MACAULEY FUNCTION DOES NOT WORK
C      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MCRD
      REAL(8), DIMENSION(MCRD), INTENT(IN) :: DELTA
      REAL(8), INTENT(IN) :: DAMAGE0
C     MATERIAL PARAMETERS
      REAL(8), INTENT(IN) :: KI
      REAL(8), INTENT(IN) :: KII
      REAL(8), INTENT(IN) :: SII
      REAL(8), INTENT(IN) :: GCII
      REAL(8), INTENT(OUT) :: DAMAGE
      REAL(8), DIMENSION(MCRD), INTENT(OUT) :: T
      REAL(8), DIMENSION(MCRD,MCRD), INTENT(OUT) :: DTDDELTA
C
C     VARIABLES WITHIN
      REAL(8) :: DI, DII, DELTA0, DELTAF
      REAL(8) :: DELTAMAX, LAMBDA
      INTEGER :: F
      REAL(8) :: D0(MCRD,MCRD), D(MCRD,MCRD), H
      REAL(8) :: B, A(MCRD,MCRD), C(MCRD)
      INTEGER :: I, J, K
C
C     MATERIAL CONSTANTS
C
C     MAXIMUM DISPLACEMENT BEFORE FAILURE
      DELTA0 = SII / KII
C
C     DISPLACEMENT AT FAILURE
      DELTAF = 2. * GCII / SII
C     
C     MODE-I
      DI = DELTA(MCRD)
C
C      WRITE(*,*) 'DELTA', DELTA
C
C     MODE-II
      DII=0.
      DO I=1,MCRD-1
          DII=DII+DELTA(I)**2  
      END DO
      DII=SQRT(DII)
      LAMBDA = DII
C
C
C     DAMAGE THRESHOLD
      DELTAMAX = DELTA0*DELTAF/(DELTAF-DAMAGE0*(DELTAF-DELTA0))
C      DELTAMAX = DELTA0 + DAMAGE0*(DELTAF-DELTA0)
C
C      WRITE(*,*) 'DAMAGE0', DAMAGE0
C
C      WRITE(*,*) 'DELTAMAX', DELTAMAX
C
C      WRITE(*,*) 'DELTAF', DELTAF
C
C      WRITE(*,*) 'DELTA0', DELTA0
C
C
C     IF THE DAMAGE THRESHOLD IS EXCEEDED
      IF (DII.GT.DELTAMAX) THEN
C
C         LOADING
          F = 1
C
          DELTAMAX = LAMBDA
C
      ELSE
C
C         UNLOADING
          F = 0
C
      END IF
C
C
C     UPDATE DAMAGE
      DAMAGE = DELTAF*(DELTAMAX-DELTA0)/DELTAMAX/(DELTAF-DELTA0)
C
C
C     CHECK THE DAMAGE PARAMETER
      IF (DAMAGE.GT.1.) THEN
C
          DAMAGE = 1.
          F = 0
C
      END IF
C
C
C     UNDAMAGED STIFFNESS
      D0 = 0.
      DO I=1,MCRD
          IF (I.EQ.MCRD) THEN
              D0(I,I)=KI
          ELSE
              D0(I,I)=KII
          END IF
      END DO
C
C
C*******************************************
!C     THIS AVOIDS THE CONVERGENCE ISSUES
      D = (1.-DAMAGE)*D0
      D(MCRD,MCRD)=(1.0-DAMAGE0)*D0(MCRD,MCRD)
C*******************************************
C
C
C     CLOSURE - RETAINS BACK THE OPENING STIFFNESS
      IF (DI.LT.0.) THEN
          D(MCRD,MCRD) = KI
      END IF
C
C
C
      DTDDELTA=D
C
C
C     ADD THE REMAINING TERMS
      IF (F.EQ.1) THEN
C
C
          B = MACAULEY(-DI)/DI
C
          A = 0.
          DO I=1,MCRD
              A(I,I)=1.    
          END DO
          A(MCRD,MCRD) = A(MCRD,MCRD) + B
C
C
          C=0.
          DO I=1,MCRD
              DO J=1,MCRD
                  C(I)=C(I)+A(I,J)*DELTA(J)
              END DO
          END DO
C
C
C
C         TANGENT STIFFNESS CONSTANT
C          IF (LAMBDA.GT.0.) THEN
              H = DELTAF*DELTA0/(DELTAF-DELTA0)/LAMBDA**3
C          ELSE
C              H = 0.
C          END IF
          
C
C
          DO I=1,MCRD
              DO J=1,MCRD
C                  DO K=1,MCRD
                      DTDDELTA(I,J)=
     + DTDDELTA(I,J) - H*KII*C(I)*C(J)
C                  END DO
              END DO
          END DO
C
C         PENETRATION
          IF (DELTA(MCRD).LT.0.) THEN
              DO I=1,MCRD
                  DTDDELTA(MCRD,I)=D(MCRD,I)
                  DTDDELTA(I,MCRD)=D(I,MCRD)
              END DO
          END IF
C
C
      END IF
C
C
C     CHECK THE TANGENT STIFFNESS MATRIX
      DO I=1,MCRD
          DO J=1,MCRD
              IF (ABS(DTDDELTA(I,J)).GT.D(I,J))
     + DTDDELTA(I,J) = D(I,J)
          END DO
      END DO
C
C
C     CALCULATE TRACTION
      T=0.
      DO I=1,MCRD
C          DO J=1,MCRD
C              T(I) = T(I) + D(I,J)*DELTA(J)
              T(I) = D(I,I)*DELTA(I)
C          END DO
      END DO
C
C
C      WRITE(*,*) 'T'
C      WRITE(*,*) T
C
C
      RETURN
      END SUBROUTINE COHESIVELAW_SHEARMODE
C      
C
C
C     COHESIVE MODEL FOR MIXED MODE
      SUBROUTINE COHESIVELAW_MIXEDMODE(
     + MCRD,DELTA,DAMAGE0,
     + KI,KII,SI,SII,
     + GCI,GCII,ETA,
     + DAMAGE,T,DTDDELTA)
C     COMMENT OUT THE IMPLICIT NONE
C     OTHERWISE MACAULEY FUNCTION DOES NOT WORK
C      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MCRD
      REAL(8), DIMENSION(MCRD), INTENT(IN) :: DELTA
      REAL(8), INTENT(IN) :: DAMAGE0
C     MATERIAL PARAMETERS
      REAL(8), INTENT(IN) :: KI
      REAL(8), INTENT(IN) :: KII
      REAL(8), INTENT(IN) :: SI
      REAL(8), INTENT(IN) :: SII
      REAL(8), INTENT(IN) :: GCI
      REAL(8), INTENT(IN) :: GCII
      REAL(8), INTENT(IN) :: ETA
      REAL(8), INTENT(OUT) :: DAMAGE
      REAL(8), DIMENSION(MCRD), INTENT(OUT) :: T
      REAL(8), DIMENSION(MCRD,MCRD), INTENT(OUT) :: DTDDELTA
C
C     VARIABLES WITHIN
      REAL(8) :: DI, DII, ODI, ODII 
      REAL(8) :: DELTA0, DELTAF
      REAL(8) :: DELTAMAX
      INTEGER :: F
      REAL(8) :: D0(MCRD,MCRD), D(MCRD,MCRD), H
      REAL(8) :: B, A(MCRD,MCRD), C(MCRD)
      REAL(8) :: ALPHA, BETA, LAMBDA
      REAL(8) :: CST, PEN
      INTEGER :: I, J, K
C
C     MATERIAL CONSTANTS
C
C     OPENING DISPLACEMENTS
      ODI = SI / KI
      ODII = SII / KII
C
C     
C     MODE-I
      DI = DELTA(MCRD)
C
C      WRITE(*,*) 'DELTA', DELTA
C
C     MODE-II
      DII=0.
      DO I=1,MCRD-1
          DII=DII+DELTA(I)**2  
      END DO
      DII = SQRT(DII)
C
C
C     DETERMINE MODE MIXITY
      IF (DI.LT.1.0D-19) THEN !MODE II
          BETA = 1.
          LAMBDA = DII
      ELSE
          BETA = DII/(DI+DII)
          LAMBDA = SQRT(DI*DI+DII*DII)
      ENDIF
C
C      WRITE(*,*) 'BETA'
C      WRITE(*,*) BETA
C
C     PENALTY STIFFNESS FOR MIXED MODE
      PEN = (1.-BETA)*KI+BETA*KII
C
C     MAXIMUM DISPLACEMENT BEFORE FAILURE
C     B-K CRITERION
      CST = (BETA**2/(1.+2.*BETA**2-2.*BETA))**ETA
      DELTA0 = SQRT(ODI**2+(ODII**2-ODI**2)*CST)
      DELTAF = 2.*(GCI+(GCII-GCI)*CST)/DELTA0/PEN      

C
C
C     DAMAGE THRESHOLD
      DELTAMAX = DELTA0*DELTAF/(DELTAF-DAMAGE0*(DELTAF-DELTA0))
C
C      WRITE(*,*) 'DAMAGE0', DAMAGE0
C
C      WRITE(*,*) 'DELTAMAX', DELTAMAX
C
C      WRITE(*,*) 'DELTAF', DELTAF
C
C      WRITE(*,*) 'DELTA0', DELTA0
C
C
C     IF THE DAMAGE THRESHOLD IS EXCEEDED
      IF (LAMBDA.GT.DELTAMAX) THEN
C
C         LOADING
          F = 1
C
          DELTAMAX = LAMBDA
C
      ELSE
C
C         UNLOADING
          F = 0
C
      END IF
C
C
C     UPDATE DAMAGE
      DAMAGE = DELTAF*(DELTAMAX-DELTA0)/DELTAMAX/(DELTAF-DELTA0)
C
C
C
C     CHECK THE DAMAGE PARAMETER
      IF (DAMAGE.GT.1.) THEN
C
          DAMAGE = 1.
          F = 0
C
      END IF
C
C      WRITE(*,*) 'DAMAGE', DAMAGE
C
C     UNDAMAGED STIFFNESS
      D0 = 0.
      DO I=1,MCRD
          IF (I.EQ.MCRD) THEN
              D0(I,I)=KI
          ELSE
              D0(I,I)=KII
          END IF
      END DO
C
C
C     STIFFNESS MATRIX
C     THIS DID NOT HELP
      D = (1.-DAMAGE)*D0
C
C
C     CLOSURE - RETAINS BACK THE OPENING STIFFNESS
      IF (DI.LT.0.) THEN
          D(MCRD,MCRD) = KI
      END IF
C
C
C
      DTDDELTA=D
C
C
C     ADD THE REMAINING TERMS
      IF (F.EQ.1) THEN
C
C
          B = MACAULEY(-DI)/DI
C
          A = 0.
          DO I=1,MCRD
              A(I,I)=1.    
          END DO
          A(MCRD,MCRD) = A(MCRD,MCRD) + B
C
C
          C=0.
          DO I=1,MCRD
              DO J=1,MCRD
                  C(I)=C(I)+A(I,J)*DELTA(J)
              END DO
          END DO
C
C
C
C         TANGENT STIFFNESS CONSTANT
          H = DELTAF*DELTA0/(DELTAF-DELTA0)/LAMBDA**3
C
C
          DO I=1,MCRD
              DO J=1,MCRD
C                  DO K=1,MCRD
                      DTDDELTA(I,J)=
     + D(I,J) - H*PEN*C(I)*C(J)
C                  END DO
              END DO
          END DO
C
C         PENETRATION
          IF (DELTA(MCRD).LT.0.) THEN
              DO I=1,MCRD
                  DTDDELTA(MCRD,I)=D(MCRD,I)
                  DTDDELTA(I,MCRD)=D(I,MCRD)
              END DO
          END IF
C
C
      END IF
C
C
C     CHECK THE TANGENT STIFFNESS MATRIX
      DO I=1,MCRD
          DO J=1,MCRD
              IF (ABS(DTDDELTA(I,J)).GT.D(I,J))
     + DTDDELTA(I,J) = D(I,J)
          END DO
      END DO
C
C     CALCULATE TRACTION
      T=0.
      DO I=1,MCRD
C          DO J=1,MCRD
C              T(I) = T(I) + D(I,J)*DELTA(J)
          T(I) = D(I,I)*DELTA(I)
C          END DO
      END DO
C
C      WRITE(*,*) 'T'
C      WRITE(*,*) T
C
C
C
      RETURN
      END SUBROUTINE COHESIVELAW_MIXEDMODE
C      
C
C
C
C     SURFACE JACOBIAN CALCULATION FOR SURFACE/LINES
      SUBROUTINE SURFJAC(MCRD,NSNODE,DNDA,DNDB,CRDS,
     + SVEC,TVEC,NVEC,RJS)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MCRD
      INTEGER, INTENT(IN) :: NSNODE
      REAL(8), DIMENSION(NSNODE), INTENT(IN) :: DNDA
      REAL(8), DIMENSION(NSNODE), INTENT(IN) :: DNDB
      REAL(8), DIMENSION(MCRD,NSNODE), INTENT(IN) :: CRDS
C     TANGENT DIRECTION
      REAL(8), DIMENSION(MCRD), INTENT(OUT) :: SVEC
C     TRANSVERSE DIRECTION
      REAL(8), DIMENSION(MCRD), INTENT(OUT) :: TVEC
C     NORMAL DIRECTION
      REAL(8), DIMENSION(MCRD), INTENT(OUT) :: NVEC
C     SURFACE OR LINE JACOBIAN (AREA OR LENGTH)
      REAL(8), INTENT(OUT) :: RJS
C
C
C     VARIABLES USED WITHIN THIS SUBROUTINE
      REAL(8) :: TJAC(MCRD,MCRD), MAG
!      REAL(8) :: TEMP(MCRD), DET
      INTEGER :: I, J
      
C     3D JACOBIAN
      TJAC=0.
      DO I=1, NSNODE
C
          DO J=1,MCRD
              TJAC(1,J) = TJAC(1,J) + DNDA(I)*CRDS(J,I)
              TJAC(2,J) = TJAC(2,J) + DNDB(I)*CRDS(J,I)
          END DO
C
      END DO
C
C
C
      NVEC(1) = TJAC(1,2)*TJAC(2,3) - TJAC(2,2)*TJAC(1,3)
      NVEC(2) = TJAC(2,1)*TJAC(1,3) - TJAC(1,1)*TJAC(2,3)
      NVEC(3) = TJAC(1,1)*TJAC(2,2) - TJAC(2,1)*TJAC(1,2)
C
      RJS = 0.
      DO I=1,MCRD
C         NORM OF NVEC
          RJS = RJS + NVEC(I)**2
      END DO
      RJS=SQRT(RJS)
C
      NVEC = NVEC / RJS
C
C     FORM TJAC
      MAG = 0.
      DO I=1,MCRD
          TJAC(3,I)=NVEC(I)
C         NORM OF S
          MAG = MAG + TJAC(1,I)**2
      END DO
      MAG = SQRT(MAG)
C
C     ASSIGN S
      DO I=1,MCRD
          SVEC(I) = TJAC(1,I)/MAG
      END DO
C
C     FIND T BY CROSS PRODUCT (T=NXS)
      TVEC(1) = NVEC(2)*SVEC(3) - SVEC(2)*NVEC(3)
      TVEC(2) = NVEC(3)*SVEC(1) - SVEC(3)*NVEC(1)
      TVEC(3) = NVEC(1)*SVEC(2) - SVEC(1)*NVEC(2)
C
C      
C
      RETURN
      END SUBROUTINE SURFJAC
C
C
C
C
C
C
C     SURFACE JACOBIAN CALCULATION FOR LINES
      SUBROUTINE LINEJAC(MCRD,
     + NSNODE,DNDA,CRDS,HEIGHT,
     + SVEC,NVEC,RJS)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MCRD
      INTEGER, INTENT(IN) :: NSNODE
      REAL(8), DIMENSION(NSNODE), INTENT(IN) :: DNDA
      REAL(8), DIMENSION(MCRD,NSNODE), INTENT(IN) :: CRDS
      REAL(8), INTENT(IN) :: HEIGHT
C     TANGENT DIRECTION
      REAL(8), DIMENSION(MCRD), INTENT(OUT) :: SVEC
C     NORMAL DIRECTION
      REAL(8), DIMENSION(MCRD), INTENT(OUT) :: NVEC
C     SURFACE OR LINE JACOBIAN (AREA OR LENGTH)
      REAL(8), INTENT(OUT) :: RJS
C
C
C     VARIABLES USED WITHIN THIS SUBROUTINE
      REAL(8) :: TJAC(MCRD)!,TJAC2D(MCRD,MCRD), DET
      INTEGER :: I, J
      
C     1D JACOBIAN
      TJAC=0.
      DO I=1, NSNODE
C
C
          DO J=1,MCRD
              TJAC(J) = TJAC(J) + DNDA(I)*CRDS(J,I)
          END DO
C
C
      END DO
C
C      
      RJS = 0.
      DO I=1,MCRD
C         NORM OF N
          RJS = RJS + TJAC(I)**2
      END DO
      RJS=SQRT(RJS)
C
C      WRITE(*,*) 'RJS', RJS
C
C     ASSIGN OUTPUTS
C     DIMENSIONS COULD BE DIFFERENT!
      DO I=1,MCRD
          SVEC(I)=TJAC(I)/RJS
      END DO
      NVEC(1) = -SVEC(2)
      NVEC(2) = SVEC(1)
C
C
C     CORRECT FOR THE HEIGHT
C     BECAUSE FAILURE IS BASED ON MPA UNITS INSTEAD OF N/MM
      RJS=RJS*HEIGHT
C
      RETURN
      END SUBROUTINE LINEJAC
C
C
C
C
C     MAPPING FOR:
C     NODE DISPLACEMENTS TO MIDSRUFACE DISPLACEMENTS   
      SUBROUTINE BMATRX(NDOFEL,MCRD,NNODE,NSNODE,N,B)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDOFEL
      INTEGER, INTENT(IN) :: MCRD
      INTEGER, INTENT(IN) :: NNODE
      INTEGER, INTENT(IN) :: NSNODE
      REAL(8), DIMENSION(NSNODE), INTENT(IN) :: N
      REAL(8), DIMENSION(MCRD,MCRD*NNODE), INTENT(OUT) :: B
C
      INTEGER :: I, J
C
      B=0.

C     ELEMENT TYPE SPECIFIC PROPERTIES
      SELECT CASE(NDOFEL)
C         C2D4
          CASE(8)
              B(1,1) =-N(1)
              B(1,3) =-N(2)
              B(1,5) = N(1)
              B(1,7) = N(2)
              B(2,2) =-N(1)
              B(2,4) =-N(2)
              B(2,6) = N(1)
              B(2,8) = N(2)
C         C2D6
          CASE(12)
              B(1,1) =-N(1)
              B(1,3) =-N(2)
              B(1,5) =-N(3)
              B(1,7) = N(1)
              B(1,9) = N(2)
              B(1,11) = N(3)
              B(2,2) =-N(1)
              B(2,4) =-N(2)
              B(2,6) =-N(3)
              B(2,8) =  N(1)
              B(2,10) = N(2)
              B(2,12) = N(3)
C         C3D6
          CASE(18)
              B(1,1) =-N(1)
              B(1,4) =-N(2)
              B(1,7) =-N(3)
              B(1,10)= N(1)
              B(1,13)= N(2)
              B(1,16)= N(3)
              B(2,2) =-N(1)
              B(2,5) =-N(2)
              B(2,8) =-N(3)
              B(2,11)= N(1)
              B(2,14)= N(2)
              B(2,17)= N(3)
              B(3,3) =-N(1)
              B(3,6) =-N(2)
              B(3,9) =-N(3)
              B(3,12)= N(1)
              B(3,15)= N(2)
              B(3,18)= N(3)
C         C3D8
          CASE(24)
              B(1,1) =-N(1)
              B(1,4) =-N(2)
              B(1,7) =-N(3)
              B(1,10)=-N(4)
              B(1,13)= N(1)
              B(1,16)= N(2)
              B(1,19)= N(3)
              B(1,22)= N(4)
              B(2,2) =-N(1)
              B(2,5) =-N(2)
              B(2,8) =-N(3)
              B(2,11)=-N(4)
              B(2,14)= N(1)
              B(2,17)= N(2)
              B(2,20)= N(3)
              B(2,23)= N(4)
              B(3,3) =-N(1)
              B(3,6) =-N(2)
              B(3,9) =-N(3)
              B(3,12)=-N(4)
              B(3,15)= N(1)
              B(3,18)= N(2)
              B(3,21)= N(3)
              B(3,24)= N(4)
C         C3D12
          CASE(36)
C
              B(1,1) = -N(1)
              B(1,4) = -N(2)
              B(1,7) = -N(3)
              B(1,10) = -N(4)
              B(1,13) = -N(5)
              B(1,16) = -N(6)
              B(1,19) = N(1)
              B(1,22) = N(2)
              B(1,25) = N(3)
              B(1,28) = N(4)
              B(1,31) = N(5)
              B(1,34) = N(6)
C
              B(2,2) = -N(1)
              B(2,5) = -N(2)
              B(2,8) = -N(3)
              B(2,11) = -N(4)
              B(2,14) = -N(5)
              B(2,17) = -N(6)
              B(2,20) = N(1)
              B(2,23) = N(2)
              B(2,26) = N(3)
              B(2,29) = N(4)
              B(2,32) = N(5)
              B(2,35) = N(6)
C
              B(3,3) = -N(1)
              B(3,6) = -N(2)
              B(3,9) = -N(3)
              B(3,12) = -N(4)
              B(3,15) = -N(5)
              B(3,18) = -N(6)
              B(3,21) = N(1)
              B(3,24) = N(2)
              B(3,27) = N(3)
              B(3,30) = N(4)
              B(3,33) = N(5)
              B(3,36) = N(6)
C
C         C3D16
          CASE(48)
C
              B(1,1) = -N(1)
              B(1,4) = -N(2)
              B(1,7) = -N(3)
              B(1,10)= -N(4)
C
              B(1,13)= -N(5)
              B(1,16)= -N(6)
              B(1,19)= -N(7)
              B(1,22)= -N(8)
C
              B(1,25)= N(1)
              B(1,28)= N(2)
              B(1,31)= N(3)
              B(1,34)= N(4)
C
              B(1,37)= N(5)
              B(1,40)= N(6)
              B(1,43)= N(7)
              B(1,46)= N(8)
C
              B(2,2) = -N(1)
              B(2,5) = -N(2)
              B(2,8) = -N(3)
              B(2,11)= -N(4)
C
              B(2,14)= -N(5)
              B(2,17)= -N(6)
              B(2,20)= -N(7)
              B(2,23)= -N(8)
C
              B(2,26)= N(1)
              B(2,29)= N(2)
              B(2,32)= N(3)
              B(2,35)= N(4)
C
              B(2,38)= N(5)
              B(2,41)= N(6)
              B(2,44)= N(7)
              B(2,47)= N(8)
C
              B(3,3) = -N(1)
              B(3,6) = -N(2)
              B(3,9) = -N(3)
              B(3,12)= -N(4)
C
              B(3,15)= -N(5)
              B(3,18)= -N(6)
              B(3,21)= -N(7)
              B(3,24)= -N(8)
C
              B(3,27)= N(1)
              B(3,30)= N(2)
              B(3,33)= N(3)
              B(3,36)= N(4)
C
              B(3,39)= N(5)
              B(3,42)= N(6)
              B(3,45)= N(7)
              B(3,48)= N(8)
C
      END SELECT
C
C
C
      RETURN
      END SUBROUTINE BMATRX      
C
C
C
C
C     ASSIGN BASIC ELEMENT PROPERTIES
C     DEPENDS ON THE ELEMENT TYPE
      SUBROUTINE ELPROP(INTMTD,NDOFEL,MCRD,
     + NSQPT,NSNODE)
      IMPLICIT NONE
C     INPUTS
C     INTEGRATION METHOD
      INTEGER, INTENT(IN) :: INTMTD
C     TOTAL NUMBER OF DEGREES OF FREEDOM OF A SINGLE ELEMENT
C     NDOFEL = NODES X DOFS
      INTEGER, INTENT(IN) :: NDOFEL
C     NUMBER OF DIMENSIONS
      INTEGER, INTENT(IN) :: MCRD
C     OUTPUTS
C     NUMBER OF INTEGRATION POINTS (IP) PER COHESIVE ELEMENT
      INTEGER, INTENT(OUT) :: NSQPT
C     NUMBER OF NODES PER SURFACE
      INTEGER, INTENT(OUT) :: NSNODE
C     
C     VARIABLES USED IN THIS SUBROUTINE
C
C
C     BASED ON THE NUMERICAL INTEGRATION METHOD
C     GAUSS-QUADRATURE      
      IF (INTMTD.EQ.1) THEN
C         ELEMENT TYPE SPECIFIC PROPERTIES
          SELECT CASE(NDOFEL)
C         C2D4
          CASE(8)
              NSQPT = 2
              NSNODE = 2
C         C2D6
          CASE(12)
              NSQPT = 3
              NSNODE = 3
C         C3D6
          CASE(18)
              NSQPT = 1
              NSNODE = 3
C         C3D8
          CASE(24)
              NSQPT = 4
              NSNODE = 4
C         C3D12
          CASE(36)
              NSQPT = 3
              NSNODE = 6
C         C3D16
          CASE(48)
              NSQPT = 9
              NSNODE = 8
          END SELECT
C     NEWTON-COTES 
      ELSEIF (INTMTD.EQ.2) THEN
C         ELEMENT TYPE SPECIFIC PROPERTIES
          SELECT CASE(NDOFEL)
C         C2D4
          CASE(8)
              NSQPT = 2
              NSNODE = 2
C         C2D6
          CASE(12)
              NSQPT = 3
              NSNODE = 3
C         C3D6
          CASE(18)
              NSQPT = 3
              NSNODE = 3
C         C3D8
          CASE(24)
              NSQPT = 4
              NSNODE = 4
C         C3D12
          CASE(36)
              NSQPT = 6
              NSNODE = 6
C         C3D16
          CASE(48)
              NSQPT = 8
              NSNODE = 8
          END SELECT
C     NEWTON-COTES (CONSISTENT WITH ABAQUS)
C     QUAD POINTS AT THE CORNER NODES)
      ELSEIF (INTMTD.EQ.3) THEN
C         ELEMENT TYPE SPECIFIC PROPERTIES
          SELECT CASE(NDOFEL)
C         C2D4
          CASE(8)
              NSQPT = 2
              NSNODE = 2
C         C2D6
          CASE(12)
              NSQPT = 2
              NSNODE = 3
C         C3D6
          CASE(18)
              NSQPT = 3
              NSNODE = 3
C         C3D8
          CASE(24)
              NSQPT = 4
              NSNODE = 4
C         C3D12
          CASE(36)
              NSQPT = 3
              NSNODE = 6
C         C3D16
          CASE(48)
              NSQPT = 4
              NSNODE = 8
          END SELECT
          END IF
C
C
C
C
C    
C
      RETURN
      END SUBROUTINE ELPROP
C
C
C
C     INTEGRATION POINT COORDINATES AND WEIGHTS
      SUBROUTINE INTEGRATIONPOINTS(INTMTD,NDOFEL,NSQPT,MCRD,
     + IPCOORDS, WT)
      IMPLICIT NONE
C     INPUTS
C     INTEGRATION METHOD
      INTEGER, INTENT(IN) :: INTMTD
C     TOTAL NUMBER OF DEGREES OF FREEDOM OF A SINGLE ELEMENT
C     NDOFEL = NODES X DOFS
      INTEGER, INTENT(IN) :: NDOFEL
C     NUMBER OF INTEGRATION POINTS (IP) PER COHESIVE ELEMENT
      INTEGER, INTENT(IN) :: NSQPT
C     NUMBER OF DIMENSIONS
      INTEGER, INTENT(IN) :: MCRD
C     OUTPUTS
C     INTEGRATION POINT COORDINATES
      REAL(8), INTENT(OUT) :: IPCOORDS(NSQPT,MCRD-1)
C     INTEGRATION POINT WEIGHTS
      REAL(8), INTENT(OUT) :: WT(NSQPT)
C
C
C
C     INTEGRATION POINT COORDINATES AND WEIGHTS
C     BASED ON THE NUMERICAL INTEGRATION METHOD
C     GAUSS-QUADRATURE      
      IF (INTMTD.EQ.1) THEN
C         ELEMENT TYPE SPECIFIC PROPERTIES
          SELECT CASE(NDOFEL)
C         C2D4
          CASE(8)
              IPCOORDS(1,:) = (/ -1. /)
              IPCOORDS(2,:) = (/  1. /)
              IPCOORDS = IPCOORDS/SQRT(3.)
              WT=1.
C         C2D6
          CASE(12)
              IPCOORDS(1,:) = (/ -1. /)
              IPCOORDS(2,:) = (/  1. /)
              IPCOORDS(3,:) = (/  0. /)
              IPCOORDS = IPCOORDS*SQRT(0.6)
              WT(1)=0.55555555555
              WT(2)=0.55555555555
              WT(3)=0.88888888888
C         C3D6
          CASE(18)
              IPCOORDS(1,:) = (/ 1., 1. /)
              IPCOORDS = IPCOORDS / 3.
              WT=1./2.
C         C3D8
          CASE(24)
              IPCOORDS(1,:) = (/ -1., -1. /)
              IPCOORDS(2,:) = (/  1., -1. /)
              IPCOORDS(3,:) = (/ -1.,  1. /)
              IPCOORDS(4,:) = (/  1.,  1. /)
              IPCOORDS = IPCOORDS/SQRT(3.)
              WT=1.
C         C3D12
          CASE(36)
              IPCOORDS(1,:) = (/ 1./6., 1./6. /)
              IPCOORDS(2,:) = (/ 2./3., 1./6. /)
              IPCOORDS(3,:) = (/ 1./6., 2./3. /)
              WT=1./6.
C         C3D16
          CASE(48)
              IPCOORDS(1,:) = (/ -1., -1. /)
              IPCOORDS(2,:) = (/  0., -1. /)
              IPCOORDS(3,:) = (/  1., -1. /)
              IPCOORDS(4,:) = (/ -1.,  0. /)
              IPCOORDS(5,:) = (/  0.,  0. /)
              IPCOORDS(6,:) = (/  1.,  0. /)
              IPCOORDS(7,:) = (/ -1.,  1. /)
              IPCOORDS(8,:) = (/  0.,  1. /)
              IPCOORDS(9,:) = (/  1.,  1. /)
              IPCOORDS =IPCOORDS*SQRT(0.6)
              WT(1)=25./81.;WT(3)=25./81.
              WT(7)=25./81.;WT(9)=25./81.
              WT(2)=40./81.;WT(4)=40./81.
              WT(6)=40./81.;WT(8)=40./81.
              WT(5)=64./81.
          END SELECT
C     NEWTON-COTES 
      ELSEIF (INTMTD.EQ.2) THEN
C         ELEMENT TYPE SPECIFIC PROPERTIES
          SELECT CASE(NDOFEL)
C         C2D4
          CASE(8)
              IPCOORDS(1,:) = (/ -1. /)
              IPCOORDS(2,:) = (/  1. /)
              WT=1.
C         C2D6
          CASE(12)
              IPCOORDS(1,:) = (/ -1. /)
              IPCOORDS(2,:) = (/  1. /)
              IPCOORDS(3,:) = (/  0. /)
              WT(1)=1./3.
              WT(2)=1./3.
              WT(3)=4./3.
C         C3D6
          CASE(18)
              IPCOORDS(1,:) = (/ 0., 0. /)
              IPCOORDS(2,:) = (/ 1., 0. /)
              IPCOORDS(3,:) = (/ 0., 1. /)
              WT=1./2./3.
C         C3D8
          CASE(24)
              IPCOORDS(1,:) = (/ -1., -1. /)
              IPCOORDS(2,:) = (/  1., -1. /)
              IPCOORDS(3,:) = (/ -1.,  1. /)
              IPCOORDS(4,:) = (/  1.,  1. /)
              WT=1.
C         C3D12
          CASE(36)
              IPCOORDS(1,:) = (/ 0., 0. /)
              IPCOORDS(2,:) = (/ 1., 0. /)
              IPCOORDS(3,:) = (/ 0., 1. /)
              IPCOORDS(4,:) = (/ 0.5, 0. /)
              IPCOORDS(5,:) = (/ 0., 0.5 /)
              IPCOORDS(6,:) = (/ 0.5, 0.5 /)
              WT(1:3)=1./9.
              WT(4:6)=16./9.
C         C3D16
          CASE(48)
              IPCOORDS(1,:) = (/ -1., -1. /)
              IPCOORDS(2,:) = (/  1., -1. /)
              IPCOORDS(3,:) = (/  1.,  1. /)
              IPCOORDS(4,:) = (/ -1.,  1. /)
              IPCOORDS(5,:) = (/  0., -1. /)
              IPCOORDS(6,:) = (/  1.,  0. /)
              IPCOORDS(7,:) = (/  0.,  1. /)
              IPCOORDS(8,:) = (/ -1.,  0. /)
              WT(1:4)=1./9.
              WT(5:8)=16./9.
          END SELECT
C     NEWTON-COTES (CONSISTENT WITH ABAQUS)
C     QUAD POINTS AT THE CORNER NODES)
      ELSEIF (INTMTD.EQ.3) THEN
C         ELEMENT TYPE SPECIFIC PROPERTIES
          SELECT CASE(NDOFEL)
C         C2D4
          CASE(8)
              IPCOORDS(1,:) = (/ -1. /)
              IPCOORDS(2,:) = (/  1. /)
              WT=1.
C         C2D6
          CASE(12)
              IPCOORDS(1,:) = (/ -1. /)
              IPCOORDS(2,:) = (/  1. /)
              WT=1.
C         C3D6
          CASE(18)
              IPCOORDS(1,:) = (/ 0., 0. /)
              IPCOORDS(2,:) = (/ 1., 0. /)
              IPCOORDS(3,:) = (/ 0., 1. /)
              WT=1./2./3.
C         C3D8
          CASE(24)
              IPCOORDS(1,:) = (/ -1., -1. /)
              IPCOORDS(2,:) = (/  1., -1. /)
              IPCOORDS(3,:) = (/ -1.,  1. /)
              IPCOORDS(4,:) = (/  1.,  1. /)
              WT=1.
C         C3D12
          CASE(36)
              IPCOORDS(1,:) = (/ 0., 0. /)
              IPCOORDS(2,:) = (/ 1., 0. /)
              IPCOORDS(3,:) = (/ 0., 1. /)
              WT=1./2./3.
C         C3D16
          CASE(48)
              IPCOORDS(1,:) = (/ -1., -1. /)
              IPCOORDS(2,:) = (/  1., -1. /)
              IPCOORDS(3,:) = (/ -1.,  1. /)
              IPCOORDS(4,:) = (/  1.,  1. /)
              WT=1.
          END SELECT          
          END IF
C
C
      RETURN
      END SUBROUTINE INTEGRATIONPOINTS
C
C
C
C
C     FINITE ELEMENT SHAPE FUNCTIONS AND DERIVATIVES
C     FOR ELEMENT TYPE: COH2D4
      SUBROUTINE COH2D4_N_DN(NSNODE,NUMDIM,G,N,DN)
      IMPLICIT NONE
C     INPUTS
C     NUMBER OF NODES PER SURFACE
      INTEGER, INTENT(IN) :: NSNODE
C     DIMENSIONS OF THE ANALYSIS
      INTEGER, INTENT(IN) :: NUMDIM
C     PARAMETRIC COORDINATES
      REAL(8), INTENT(IN) :: G
C
C     OUTPUTS
C     SHAPE FUNCTIONS
      REAL(8), DIMENSION(NSNODE), INTENT(OUT) :: N
C     SHAPE FUNCTION DERIVATIVES
      REAL(8), DIMENSION(NUMDIM,NSNODE), INTENT(OUT) :: DN
C
C
C     SHAPE FUNCTIONS FOR CP4 ELEMENT
      N(1) = (1. - G) / 2.
C
      N(2) = (1. + G) / 2.
C
C
C
C
C     SHAPE FUNCTION DERIVATIVES WRT NATURAL COORDS FOR CP4 ELEMENT
C
C     DN_DG
      DN(1,1) = -1. / 2.
C
      DN(1,2) =  1. / 2.
C
C
C
C
      RETURN
      END SUBROUTINE COH2D4_N_DN
C
C
C
C
C
C     FINITE ELEMENT SHAPE FUNCTIONS AND DERIVATIVES
C     FOR ELEMENT TYPE: COH2D6
      SUBROUTINE COH2D6_N_DN(NSNODE,NUMDIM,G,N,DN)
      IMPLICIT NONE
C     INPUTS
C     NUMBER OF NODES PER SURFACE
      INTEGER, INTENT(IN) :: NSNODE
C     DIMENSIONS OF THE ANALYSIS
      INTEGER, INTENT(IN) :: NUMDIM
C     PARAMETRIC COORDINATES
      REAL(8), INTENT(IN) :: G
C
C     OUTPUTS
C     SHAPE FUNCTIONS
      REAL(8), DIMENSION(NSNODE), INTENT(OUT) :: N
C     SHAPE FUNCTION DERIVATIVES
      REAL(8), DIMENSION(NUMDIM,NSNODE), INTENT(OUT) :: DN
C
C
C     SHAPE FUNCTIONS FOR CP8 ELEMENT
      N(1) = G * (G - 1.) / 2.
C
      N(2) = G * (G + 1.) / 2.
C
      N(3) = (1 + G) * (1 - G)
C
C
C
C
C     SHAPE FUNCTION DERIVATIVES WRT NATURAL COORDS FOR CP4 ELEMENT
C
C     DN_DG
      DN(1,1) = G - 1. / 2.
C
      DN(1,2) = G + 1. / 2.
C
      DN(1,3) = -2. * G
C
C
C
C
      RETURN
      END SUBROUTINE COH2D6_N_DN
C
C
C
C
C
C
C     FINITE ELEMENT SHAPE FUNCTIONS AND DERIVATIVES
C     FOR ELEMENT TYPE: COH3D6
      SUBROUTINE COH3D6_N_DN(NNPEL,NUMDIM,G,H,N,DN)
      IMPLICIT NONE
C     INPUTS
C     NUMBER OF NODES PER ELEMENT
      INTEGER, INTENT(IN) :: NNPEL
C     DIMENSIONS OF THE ANALYSIS
      INTEGER, INTENT(IN) :: NUMDIM
C     PARAMETRIC COORDINATES
      REAL(8), INTENT(IN) :: G
      REAL(8), INTENT(IN) :: H
C
C     OUTPUTS
C     SHAPE FUNCTIONS
      REAL(8), DIMENSION(NNPEL), INTENT(OUT) :: N
C     SHAPE FUNCTION DERIVATIVES
      REAL(8), DIMENSION(NUMDIM,NNPEL), INTENT(OUT) :: DN
C
C
C     SHAPE FUNCTIONS FOR CP3 ELEMENT
      N(1) = 1. - G - H
C
      N(2) = G
C
      N(3) = H
C
C
C
C
C
C     SHAPE FUNCTION DERIVATIVES WRT NATURAL COORDS FOR CP3 ELEMENT
C     DN_DG
      DN(1,1) = -1.
C
      DN(1,2) = 1.
C
      DN(1,3) = 0.
C
C     DN_DH      
      DN(2,1) = -1.
C
      DN(2,2) = 0.
C
      DN(2,3) = 1.
C
C
C
C
      RETURN
      END SUBROUTINE COH3D6_N_DN
C
C
C
C
C
C
C
C
C
C
C     FINITE ELEMENT SHAPE FUNCTIONS AND DERIVATIVES
C     FOR ELEMENT TYPE: COH3D8
      SUBROUTINE COH3D8_N_DN(NSNODE,NUMDIM,G,H,N,DN)
      IMPLICIT NONE
C     INPUTS
C     NUMBER OF NODES PER SURFACE
      INTEGER, INTENT(IN) :: NSNODE
C     DIMENSIONS OF THE ANALYSIS
      INTEGER, INTENT(IN) :: NUMDIM
C     PARAMETRIC COORDINATES
      REAL(8), INTENT(IN) :: G
      REAL(8), INTENT(IN) :: H
C
C     OUTPUTS
C     SHAPE FUNCTIONS
      REAL(8), DIMENSION(NSNODE), INTENT(OUT) :: N
C     SHAPE FUNCTION DERIVATIVES
      REAL(8), DIMENSION(NUMDIM,NSNODE), INTENT(OUT) :: DN
C
C
C     SHAPE FUNCTIONS FOR CP4 ELEMENT
      N(1) = (1. - G) * (1. - H) / 4.
C
      N(2) = (1. + G) * (1. - H) / 4.
C
      N(3) = (1. + G) * (1. + H) / 4.
C
      N(4) = (1. - G) * (1. + H) / 4.
C
C
C
C
C     SHAPE FUNCTION DERIVATIVES WRT NATURAL COORDS FOR CP4 ELEMENT
C
C     DN_DG
      DN(1,1) = -1. * (1. - H) / 4.
C
      DN(1,2) = 1. * (1. - H) / 4.
C
      DN(1,3) = 1. * (1. + H) / 4.
C
      DN(1,4) = -1. * (1. + H) / 4.
C
C     DN_DH
      DN(2,1) = (1. - G) * -1. / 4.
C
      DN(2,2) = (1. + G) * -1. / 4.
C
      DN(2,3) = (1. + G) * 1. / 4.
 
      DN(2,4) = (1. - G) * 1. / 4.
C
C
C
C
C
      RETURN
      END SUBROUTINE COH3D8_N_DN
C
C
C
C
C
C
C
C
C
C     FINITE ELEMENT SHAPE FUNCTIONS AND DERIVATIVES
C     FOR ELEMENT TYPE: COH3D12
      SUBROUTINE COH3D12_N_DN(NNPEL,NUMDIM,G,H,N,DN)
      IMPLICIT NONE
C     INPUTS
C     NUMBER OF NODES PER ELEMENT
      INTEGER, INTENT(IN) :: NNPEL
C     DIMENSIONS OF THE ANALYSIS
      INTEGER, INTENT(IN) :: NUMDIM
C     PARAMETRIC COORDINATES
      REAL(8), INTENT(IN) :: G
      REAL(8), INTENT(IN) :: H
C
C     OUTPUTS
C     SHAPE FUNCTIONS
      REAL(8), DIMENSION(NNPEL), INTENT(OUT) :: N
C     SHAPE FUNCTION DERIVATIVES
      REAL(8), DIMENSION(NUMDIM,NNPEL), INTENT(OUT) :: DN
C
C
C     SHAPE FUNCTIONS FOR CP6 ELEMENT
      N(1) = 2. * (0.5 - G - H) * (1. - G - H)
C
      N(2) = 2. * G * (G - 0.5)
C
      N(3) = 2. * H * (H - 0.5)
C
      N(4) = 4. * G * (1. - G - H)
C
      N(5) = 4. * G * H
C
      N(6) = 4. * H * (1. - G - H)
C
C
C
C
C     SHAPE FUNCTION DERIVATIVES WRT NATURAL COORDS FOR C3D4 ELEMENT
C     DN_DG
      DN(1,1) = 2. * (-1.) * (1. - G - H) + 2. * (0.5 - G - H) * (-1.)
C
      DN(1,2) = 2. * (1.) * (G - 0.5) + 2. * G * (1.)
C
      DN(1,3) = 0.
C
      DN(1,4) = 4. * (1.) * (1. - G - H) + 4. * G * (-1.)
C
      DN(1,5) = 4. * (1.) * H
C
      DN(1,6) = 4. * H * (-1.)
C
C     DN_DH       
      DN(2,1) = 2. * (-1.) * (1. - G - H) + 2. * (0.5 - G - H) * (-1.)
C
      DN(2,2) = 0.
C
      DN(2,3) = 2. * (1.) * (H - 0.5)  + 2. * H * (1.)
C
      DN(2,4) = 4. * G * (-1.)
C
      DN(2,5) = 4. * G * (1.)
C
      DN(2,6) = 4. * (1.) * (1. - G - H) + 4. * H * (-1.)
C
C
C
C
      RETURN
      END SUBROUTINE COH3D12_N_DN
C
C
C
C
C
C     FINITE ELEMENT SHAPE FUNCTIONS AND DERIVATIVES
C     FOR ELEMENT TYPE: CPE8 OR CPS8
      SUBROUTINE COH3D16_N_DN(NSNODE,NUMDIM,G,H,N,DN)
      IMPLICIT NONE
C     INPUTS
C     NUMBER OF NODES PER SURFACE
      INTEGER, INTENT(IN) :: NSNODE
C     DIMENSIONS OF THE ANALYSIS
      INTEGER, INTENT(IN) :: NUMDIM
C     PARAMETRIC COORDINATES
      REAL(8), INTENT(IN) :: G
      REAL(8), INTENT(IN) :: H
C
C     OUTPUTS
C     SHAPE FUNCTIONS
      REAL(8), DIMENSION(NSNODE), INTENT(OUT) :: N
C     SHAPE FUNCTION DERIVATIVES
      REAL(8), DIMENSION(NUMDIM,NSNODE), INTENT(OUT) :: DN
C
C
C     SHAPE FUNCTIONS FOR CP8 ELEMENT
      N(1) = -1./4. * (1. - G) * (1. - H) * (1. + G + H)
C
      N(2) = -1./4. * (1. + G) * (1. - H) * (1. - G + H)
C
      N(3) = -1./4. * (1. + G) * (1. + H) * (1. - G - H)
C
      N(4) = -1./4. * (1. - G) * (1. + H) * (1. + G - H)
C
      N(5) = 1./2. * (1. - G) * (1. + G) * (1. - H)
C
      N(6) = 1./2. * (1. - H) * (1. + H) * (1. + G)
C
      N(7) = 1./2. * (1. - G) * (1. + G) * (1. + H)
C
      N(8) = 1./2. * (1. - H) * (1. + H) * (1. - G)
C
C
C
C
C
C     SHAPE FUNCTION DERIVATIVES WRT NATURAL COORDS FOR CP8 ELEMENT
C     DN_DG
      DN(1,1) = (-1./4.) * (-1.) * (1. - H) * (1. + G + H) +
     + (-1./4.) * (1. - G) * (1. - H) * (1.)
C
      DN(1,2) = (-1./4.) * (1.) *  (1. - H) * (1. - G + H) +
     + (-1./4.) * (1. + G) * (1. - H) * (-1.)
C
      DN(1,3) = (-1./4.) * (1.) *  (1. + H) * (1. - G - H) +
     + (-1./4.) * (1. + G) * (1. + H) * (-1.)
C
      DN(1,4) = (-1./4.) * (-1.) * (1. + H) * (1. + G - H) +
     + (-1./4.) * (1. - G) * (1. + H) * (1.)
C
      DN(1,5) = 1./2. * (-1.) * (1. + G) * (1. - H) +
     + 1./2. * (1. - G) * (1.) * (1. - H)
C
      DN(1,6) = 1./2. * (1. - H) * (1. + H) * (1.)
C
      DN(1,7) = 1./2. * (-1.) * (1. + G) * (1. + H) +
     + 1./2. * (1. - G) * (1.) * (1. + H)
C
      DN(1,8) = 1./2. * (1. - H) * (1. + H) * (-1.)
C
C
C
C     DN_DH
      DN(2,1) = (-1./4.) * (1. - G) * (-1.) * (1. + G + H)  +
     + (-1./4.) * (1. - G) * (1. - H) * (1.)
C
      DN(2,2) = (-1./4.) * (1. + G) * (-1.) * (1. - G + H) +
     + (-1./4.) * (1. + G) * (1. - H) * (1.)
C
      DN(2,3) = (-1./4.) * (1. + G) *  (1.) * (1. - G - H) + 
     + (-1./4.) * (1. + G) * (1. + H) * (-1.)
 
      DN(2,4) = (-1./4.) * (1. - G) *  (1.) * (1. + G - H) +
     + (-1./4.) * (1. - G) * (1. + H) * (-1.)
C
      DN(2,5) = 1./2. * (1. - G) * (1. + G) * (-1.)
C
      DN(2,6) = 1./2. * (-1.) * (1. + H) * (1. + G) +
     + 1./2. * (1. - H) * (1.) * (1. + G)
C
      DN(2,7) = 1./2. * (1. - G) * (1. + G) * (1.)
C
      DN(2,8) = 1./2. * (-1.) * (1. + H) * (1. - G) +
     + 1./2. * (1. - H) * (1.) * (1. - G)
C
C          
C     
      RETURN
      END SUBROUTINE COH3D16_N_DN
C
C
C
C
C	THIS SUBROUTINE INVERTS A 3X3 MATRIX
C	INPUT:	MATRIX								---	A(3,3)
C	OUTPUT:	INVERETED MATRIX, DETERMINANT		---	INVA(3,3),DET
	SUBROUTINE INV3X3(A,INVA,DET)
	IMPLICIT NONE
      REAL(8), INTENT(IN)  :: A(3,3)
      REAL(8), INTENT(OUT) :: INVA(3,3), DET
	INTEGER :: I,J
      REAL(8) :: SMALLNUM
C
      SMALLNUM = 1.D-20
C
C
C	FIRST CALCULATE THE DETERMINANT
	CALL DETER3X3(A,DET)
C	IF THE DETERMINANT IS GREATER THAN CERTAIN VALUE
	IF (ABS(DET) < SMALLNUM) THEN
		INVA=0.0D+0
	ELSE
		INVA(1,1)=((A(2,2)*A(3,3))-(A(2,3)*A(3,2)))/DET
		INVA(2,1)=-((A(2,1)*A(3,3))-(A(2,3)*A(3,1)))/DET
		INVA(3,1)=((A(2,1)*A(3,2))-(A(2,2)*A(3,1)))/DET
		INVA(1,2)=-((A(1,2)*A(3,3))-(A(1,3)*A(3,2)))/DET
		INVA(2,2)=((A(1,1)*A(3,3))-(A(1,3)*A(3,1)))/DET
		INVA(3,2)=-((A(1,1)*A(3,2))-(A(1,2)*A(3,1)))/DET
		INVA(1,3)=((A(1,2)*A(2,3))-(A(1,3)*A(2,2)))/DET
		INVA(2,3)=-((A(1,1)*A(2,3))-(A(2,1)*A(1,3)))/DET
		INVA(3,3)=((A(1,1)*A(2,2))-(A(1,2)*A(2,1)))/DET
	ENDIF
	RETURN
      END SUBROUTINE INV3X3
C
C
C
C	THIS SUBROUTINE INVERTS A 2X2 MATRIX
C	INPUT:	MATRIX								---	A(2,2)
C	OUTPUT:	INVERETED MATRIX, DETERMINANT		---	INVA(2,2),DET
	SUBROUTINE INV2X2(A,INVA,DET)
	IMPLICIT NONE
      REAL(8), INTENT(IN)  :: A(2,2)
      REAL(8), INTENT(OUT) :: INVA(2,2), DET
	INTEGER :: I, J
      REAL(8) :: SMALLNUM
C
      SMALLNUM = 1.D-20
C
C
C	FIRST CALCULATE THE DETERMINANT
	CALL DETER2X2(A,DET)
C	IF THE DETERMINANT IS GREATER THAN CERTAIN VALUE
	IF (ABS(DET).LT.SMALLNUM) THEN
		INVA=0.0D+0
	ELSE
		INVA(1,1) = A(2,2)/DET
		INVA(1,2) =-A(1,2)/DET
		INVA(2,1) =-A(2,1)/DET
          INVA(2,2) = A(1,1)/DET
	ENDIF
	RETURN
	END SUBROUTINE INV2X2      
C
C
C
C
C *************************************************
C *      THE DETERMINANT OF A 2X2 MATRIX          *
C *************************************************
      SUBROUTINE DETER2X2(DMIN,D)
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: DMIN(2,2)
      REAL(8), INTENT(OUT) :: D
C
      D=0.
      D=DMIN(1,1)*DMIN(2,2)-DMIN(1,2)*DMIN(2,1)
C
      RETURN
      END SUBROUTINE DETER2X2
C
C
C
C
C *************************************************
C *      THE DETERMINANT OF A 3X3 MATRIX          *
C *************************************************
      SUBROUTINE DETER3X3(DMIN,D)
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: DMIN(3,3)
      REAL(8), INTENT(OUT) :: D
C
      D = 0.
      D = DMIN(1,1)*DMIN(2,2)*DMIN(3,3) + 
     + DMIN(1,2)*DMIN(2,3)*DMIN(3,1) + 
     + DMIN(2,1)*DMIN(3,2)*DMIN(1,3) -
     + DMIN(1,3)*DMIN(2,2)*DMIN(3,1) -
     + DMIN(1,1)*DMIN(2,3)*DMIN(3,2) -
     + DMIN(1,2)*DMIN(2,1)*DMIN(3,3)
C
      RETURN
      END SUBROUTINE DETER3X3
C      
C
C     
C
C
C
C
      REAL(8) FUNCTION MACAULEY(X)
      IMPLICIT NONE
      REAL(8) X
      MACAULEY = (ABS(X) + X)/2.
      RETURN
      END FUNCTION MACAULEY
C
C
C