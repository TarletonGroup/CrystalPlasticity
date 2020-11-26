********************************************
** KMAT calculates the material behaviour **
**  Global stiffness matrix and stress    **
********************************************

      SUBROUTINE kmat(dtime,nsvars,usvars,xI,jelem,kint,time,F,
     + L,iphase,irradiate,C,stressvec,dstressinc,totstran,dtotstran,
     + TEMP,DTEMP,vms,pdot,pnewdt,gndon,nSys,nTwin,ns,coords,
     + TwinIntegral,nTwinStart,nTwinEnd,twinon)

      INCLUDE 'ABA_PARAM.INC'

      ! number of Abaqus state variables
      INTEGER,intent(in) :: nsvars

      ! element number
      INTEGER,intent(in) :: jelem

      ! integration point number
      INTEGER,intent(in) :: kint

      ! crystal type
      INTEGER,intent(in) :: iphase

      ! activate irradiation effect
      INTEGER,intent(in) :: irradiate

      ! GND activation flag
      INTEGER,intent(in) :: gndon

      ! number of slip systems
      INTEGER,intent(in) :: nSys

      ! total number of twin systems
      INTEGER,intent(in) :: nTwin

      ! number of screw dislocation systems
      INTEGER,intent(in) :: ns

      ! the active twins are the ones in the
      ! interval [nTwinStart,nTwinEnd] in the
      ! twin system file
      INTEGER,intent(in) :: nTwinStart,nTwinEnd

      ! twin systems activation flag
      INTEGER,intent(in) :: twinon

      ! time increment
      REAL*8,intent(in) :: dtime

      ! coordinates
      REAL*8,intent(in) :: coords(3)

      ! average of the twin volume fraction
      ! over the neighbourhood
      ! two twin systems
      REAL*8,intent(in) :: TwinIntegral(nTwin)

      ! Abaqus state variables
      REAL*8,intent(inout) :: usvars(nsvars)

      ! time: time(1) = step time; time(2) = total time
      REAL*8,intent(in) :: time(2)

      ! deformation gradient
      REAL*8,intent(in) :: F(3,3)

      ! velocity gradient
      REAL*8,intent(in) :: L(3,3)

      ! identity matrix
      REAL*8,intent(in) :: xI(3,3)

      ! stress vector (Cauchy, Abaqus notation)
      ! stressvec = (sigma11,sigma22,sigma33,tau12 tau13 tau23)
      REAL*8,intent(out) :: stressvec(6)

      ! Stiffness matrix (Jacobian, Abaqus notation)
      REAL*8,intent(out) :: C(6,6)

      ! Von Mises invariant plastic strain rate
      REAL*8,intent(out) :: pdot

      ! Von Mises stress
      REAL*8,intent(out) :: vms

******************************************
** The following parameters must be set **

      ! activate debug mode with Visual Studio
      ! 0 = off ; 1 = on
      integer, parameter :: debug = 0 

      ! select slip law
      ! 0 = Original slip rule with no GND coupling i.e. using alpha and beta
      ! 5 = Original slip rule with GND coupling 
      ! 6 = Slip rule with constants alpha and beta: alpha.sinh[ beta(tau-tauc)sgn(tau) ]
      ! 7 = Powerlaw plasticity
      integer, parameter :: kslip = 7

      ! initial temperature and temperature rate
      real*8, parameter :: Temperature = 293.0
      real*8, parameter :: ytemprate = 0.0

      ! homogenize twin model
      ! 0 = use discrete twin model
      ! 1 = use homogenised twin model
      integer, parameter :: homogtwin = 1   

      ! accuracy of the crystal plasticity Newton-Raphson loop
      ! WARNING: change only if you know what you are doing
      real*8, parameter :: xacc = 1.e-8  

      ! max number of iterations of the Newton-Raphson loop
      ! WARNING: change only if you know what you are doing
      integer, parameter :: maxNRiter = 1000  

**       End of parameters to set       **
******************************************

      ! dimension of the space
      INTEGER, parameter :: M=3,N=3

      ! dimension of Voigt vectors
      INTEGER, parameter :: KM=6,KN=6

      ! number of active slip systems considered
      integer, parameter :: L0=12 ! HCP
      integer, parameter :: L1=12 ! BCC
      integer, parameter :: L2=12 ! FCC
      integer, parameter :: L4=7  ! Olivine
      integer, parameter :: LalphaUranium=8 ! alpha-uranium

      ! stress matrix in the Newton-Raphson loop
      REAL*8 :: stressM(3,3)

      ! plastic strain increment in the Newton-Raphson loop
      REAL*8 :: plasStrainInc(3,3),plasStrainInc2(6)

      ! temporary variables
      REAL*8 :: prod(M),prod6(6,6)

      ! temporary normals and directions of slip/twin systems
      REAL*8 :: tempNorm(M),tempDir(M)

      ! plastic strain rate
      REAL*8 :: plasStrainRate(3,3)

      ! plastic velocity gradient
      REAL*8 :: Lp(3,3)

      ! cumulative plastic strain (for output)
      ! vector and scalar
      REAL*8 :: totplasstran(6), p

      ! identity matrix
      REAL*8 :: xIden6(KM,KN)

      ! elastic velocity gradient
      REAL*8 :: Le(3,3)

      ! derivatives of the plastic velocity gradient
      ! with respect to the stress components
      REAL*8 :: tmat(KM,KN)

      ! trial stress, starting point of the
      ! Newton-Raphson loop, vector and matrix
      REAL*8 :: trialstress(6),trialstressM(3,3)

      ! Jacobian of the Newton-Raphson loop
      ! and its inverse
      REAL*8 :: xjfai(KM,KN),xjfaiinv(KM,KN)

      ! residual of the Newton-Raphson loop
      ! vector and scalar
      REAL*8 :: faivalue,fai(6)

      ! stress increment of the Newton-Raphson loop
      REAL*8 :: dstressinc(6)

      ! stress at the current increment
      ! given back to Abaqus at the end of iteration
      ! matrix and vector
      REAL*8 :: xstressmdef(M,N),xstressdef(6)

      ! elastic stiffness matrix in the crystal reference frame
      REAL*8 :: xStiff(6,6)

      ! elastic stiffness matrix in the sample reference frame
      REAL*8 :: xStiffdef(6,6)

      ! elastic spin (W_e)
      REAL*8 :: elasspin(3,3)

      ! temporary array for elastic stiffness calculation
      REAL*8 :: tSig(6,6),tSigTranspose(6,6)

      ! rotation matrix from input file
      ! current and previous increment
      REAL*8 :: gmatinv(3,3),gmatinvnew(3,3),gmatinvold(3,3)

      ! deviatoric stress (for output)
      REAL*8 :: devstress(3,3)

      ! elastic compliance matrix in the crystal reference frame 
      REAL*8 :: compliance(6,6)

      ! temporary variables for plastic deformation update
      REAL*8 :: print2(3,3),print3(3,3)

      ! total strain increment
      ! vector and matrix
      REAL*8 :: dtotstran(6),tempstrain(3,3)

      ! cumulative total strain (for output)
      REAL*8 :: totstran(6)

      ! curl of the plastic deformation gradient
      REAL*8 :: curlfp(3,3)

      ! elastic Green-Lagrange strain (for output)
      REAL*8 :: EECrys(3,3)

      ! spin W (from total velocity gradient)
      REAL*8 :: spin(3,3)

      ! plastic deformation gradient
      ! current and previous increment
      ! and inverse
      REAL*8 :: Fp(3,3),fpold(3,3),Fpinv(3,3)

      ! elastic deformation gradient
      ! and inverse
      REAL*8 :: Fe(3,3),Feinv(3,3)

      ! thermal eigenstrain to model thermal expansion
      ! Voigt vector, matrix in the crystal reference frame
      ! matrix in the sample reference frame
      REAL*8 :: dstranth(6),thermat(3,3),expanse33(3,3) 

      ! rotated slip normals and directions (gmatinv)
      REAL*8,dimension(nSys,M) :: xNorm,xDir

      ! resolved shear stress on slip system
      ! and its sign
      REAL*8,dimension(nSys) :: tau
      REAL*8,dimension(nSys) :: signtau

      ! plastic shear rate on slip systems
      REAL*8,dimension(nSys) :: gammadot

      ! GND density (immobile, mobile)
      REAL*8,dimension(nSys) :: gndcut
      REAL*8,dimension(nSys) :: gndmob

      ! Burgers vectors
      REAL*8,dimension(nSys) :: burgerv

      ! critical resolved shear stress of slip systems
      REAL*8,dimension(nSys) :: tauc

      ! resolved shear stress, twin systems
      REAL*8,dimension(nTwin) :: tautwin  

      ! sign of the resolved shear stress of twin systems
      REAL*8,dimension(nTwin) :: signtautwin

      ! critical resolved shear stress of twin systems
      REAL*8,dimension(nTwin) :: tauctwin

      ! rotated twin normal and direction (gmatinv)
      REAL*8,dimension(nTwin,M) :: xTwinNorm,xTwinDir

      ! plastic shear rate due to twinning
      REAL*8 :: gammatwindot(nTwin) 

      ! GND density edge and screw systems
      REAL*8,dimension(nSys+ns) :: gndall

      ! scalar total GND density
      REAL*8 :: gndtot

      ! prefactor for SSD evolution equation
      REAL*8 :: gammast 

      REAL*8 :: LpFeinv(3,3), matrix(3,3), update(3,3)

      ! number of screw planes 
      integer :: screwplanes

      ! current temperature
      ! modified by temperature rate
      REAL*8 :: CurrentTemperature

      ! substructure dislocation density
      REAL*8 :: rhosub 

      ! forest dislocation density
      REAL*8 :: rhofor(nSys)

      ! total dislocation density 
      REAL*8 :: rhototal

      ! total SSD density
      REAL*8 :: rhossd

      ! ratio: resolved shear stress for slip / CRSS for slip
      REAL*8 :: xtau 

      ! twin volume fraction
      REAL*8 :: twinvolfrac(nTwin)

      ! total twin volume fraction 
      REAL*8 :: twinvolfractotal

      ! ratio: resolved sheat stress for twin / CRSS for twin
      REAL*8 :: xtautwin 
      
      ! rotation matrix due to twinning in the lattice system
      ! and temporary matrix
      REAL*8 :: TwinRot(nTwin,3,3) 
      REAL*8 :: TwinRotTemp(3,3)

      ! stiffness tensor after twinning
      ! and temporary matrix
      REAL*8 :: xStifftwin(6,6)
      REAL*8 :: xStifftwinTemp(6,6)

      ! c/a ratio for hcp crystals
      REAL*8 :: caratio

      ! shear modulus for Taylor dislocation law
      REAL*8 :: G12

      ! crystallographic slip
      ! needed for model with irradiation
      REAL*8 :: slip

      ! increase in tauc due to solute force
      REAL*8 :: tauSolute

      ! temporary variable for determinant
      REAL*8 :: deter

      ! flag to decide if crystal plasticity
      ! Newton-Raphson loop starts
      integer :: EnterNRLoop

      ! iteration number of the Newton-Raphson loop
      integer :: iter
      
      ! debug temporary variable
      integer :: debugWait

      integer :: i, j, k
      
C     *** INITIALIZE ZERO ARRAYS ***
      prod=0.
      tSig=0.
      xjfai=0.
      xjfaiinv=0.
      totstran = 0.0 
      totplasstran = 0.0
      trialstress=0.
      devstress=0.
      spin=0.
      tempstrain=0.
      Fe=0.
      Fp=0.
      gmatinvnew=0.
      gmatinv=0.
      dstranth=0.
      Lp = 0.
      plasStrainRate = 0.
      plasStrainInc2=0.
      xStiff=0.0      
      xStiffdef=0.0
      C=0.0
      trialstressM=0.0
      tmat=0.0;
      plasStrainInc=0.
      stressvec=0.
      fai=0.
      dstressinc=0.
      xIden6=0.;
      xRot=0.;
      compliance=0.      
      Le = 0.
      thermat =0.
      expanse33=0.
      curlfp=0.
      gammadot=0.
      xNorm=0.
      xDir=0.
      tau=0.
      gndcut=0.
      gndall=0.
      tautwin = 0.0
      gndmob=0.
      gammatwindot = 0.0
      gndtot=0.
      xTwinNorm = 0.0
      xTwinDir = 0.0
      twinvolfrac(1:nTwin) = 0.0
      twinvolfractotal = 0.0
      caratio = 0.0
      G12 = 0.0
      gammast = 0.0
      burgerv = 0.0
      tauc = 0.0
      screwplanes = 0
      tauctwin(1:nTwin) = 0.0
      rhossd = 0.0
      pdot = 0.0

      ! define identity matrix
      DO I=1,KM; xIden6(I,I)=1.; END DO      

      ! initialize sign of the resolved shear stress
      signtau=1.
      signtautwin=1.0

      ! calculate temperature
      CurrentTemperature = Temperature + ytemprate*time(2)

      ! set materials constants
      call kMaterialParam(iphase,caratio,compliance,G12,thermat,
     + gammast,burgerv,nSys,tauc,screwplanes,CurrentTemperature,
     + tauctwin,nTwin,twinon,nTwinStart,nTwinEnd,TwinIntegral)

      ! define rotation matrices due to twinning (in the lattice system)
      TwinRot = 0.0
      if (twinon == 1) then
        CALL ktwinrot(nTwin,TwinRot)
      end if
      TwinRotTemp = 0.0

      ! find stiffness matrix
      CALL lapinverse(compliance,6,info,xStiff)

C     *** INITIALIZE USER ARRAYS FROM STATE VARIABLES ***

      ! get rotation matrix
      DO i=1,3
        DO j=1,3
          gmatinv(i,j) = usvars(j+(i-1)*3)
        END DO
      END DO
      gmatinvold = gmatinv

      ! get cumulative plastic strain rate scalar
      p = usvars(10)

      ! get cumulative plastic strain rate vector
      DO i=1,6
        totplasstran(i) = usvars(10+i)
      END DO

      ! get cumulative total strain
      DO i=1,6
        totstran(i) = usvars(16+i)
      END DO

      ! initialize stress as the values at previous increment
      DO i=1,6
        xstressdef(i) = usvars(47+i)
      END DO

      ! get scalar total GND density
      gndtot = usvars(26)  
    
      ! get SSD dislocation density
      rhossd = usvars(54)
                
      ! get immobile GND density
      DO i=1,nSys
        gndcut(i) = usvars(56+i)
      END DO

      ! model with forest and substructure dislocation densities
      ! R.J. McCabe, L. Capolungo, P.E. Marshall, C.M. Cady, C.N. Tomé
      ! Deformation of wrought uranium: Experiments and modeling
      ! Acta Materialia 58 (2010) 5447–5459
      if (iphase == 5) then
        rhototal = 0.0
        DO i=1,nSys
          rhofor(i) = usvars(56+i)
          rhototal = rhototal + rhofor(i)
        END DO
        rhosub = usvars(65)
        rhototal = rhototal + rhosub
        if (twinon == 1) then ! twin active
          ! initialize twin volume fraction
          DO i=1,nTwin
            twinvolfrac(i) = usvars(106+i)
          END DO
          DO i=1,nTwin ! calculate total twin volume fraction
            twinvolfractotal = twinvolfractotal + twinvolfrac(i)
          END DO
          twinvolfractotal = min(twinvolfractotal,1.0)
        end if ! twin active
      end if

      ! get plastic deformation gradient
      DO i=1,3
        DO j=1,3
          Fp(i,j) = usvars(80+j+((i-1)*3))
        END DO
      END DO

      ! get curl of the plastic deformation gradient
      DO i=1,3
        DO j=1,3 
         curlfp(i,j) = usvars(37+j+(i-1)*3)
        END DO
      END DO

      ! variables needed for irradiation model (softening)
      ! crystallographic slip
      slip = usvars(35)    
      ! increase in tauc due to solute force
      tauSolute = usvars(36)

      ! calculates CRSS of slip and twin systems
      call kCRSS(iphase,tauc,nSys,G12,burgerv,gndtot,irradiate,
     + tauSolute,gndcut,rhofor,rhosub,CurrentTemperature,homogtwin,
     + nTwinStart,nTwinEnd,twinvolfrac,tauctwin,nTwin,TwinIntegral,
     + twinvolfractotal,twinon)

      ! Reorient stiffness tensor if twins are present
      ! weighting using twin volume fraction
      ! calculation in the lattice reference frame
      ! initialization of the temporary variable
      xStifftwin = 0.0
      xStifftwinTemp = 0.0
      if (twinon == 1) then ! twin active
        DO k=nTwinStart,nTwinEnd
	  DO i=1,M ! get rotation matrix for this twin
            DO j=1,N
	      TwinRotTemp(i,j) = TwinRot(k,i,j)
	    END DO
          END DO
          CALL rotord4sig(TwinRotTemp,tSig) ! rotate stiffness tensor
          prod6 = matmul(tSig,xStiff)
          tSigTranspose = transpose(tSig)
          xStifftwinTemp = twinvolfrac(k)*matmul(prod6,tSigTranspose)
          xStifftwin = xStifftwin + xStifftwinTemp
        END DO ! end of twin system
      end if ! twin active
      xStiff = (1.0 - twinvolfractotal) * xStiff + xStifftwin


C     *** DIRECTIONS FROM LATTICE TO DEFORMED AND TWINNED SYSTEM ***

      CALL kdirns(gmatinv,TwinRot,iphase,nSys,nTwin,xDir,xNorm,
     + xTwinDir,xTwinNorm,caratio)


C     *** STIFFNESS FROM LATTICE TO DEFORMED SYSTEM ***

      CALL rotord4sig(gmatinv,tSig)

      prod6 = matmul(tSig,xStiff)
      tSigTranspose = transpose(tSig)

      xStiffdef = matmul(prod6,tSigtranspose)

      ! rotate the thermal eigenstrain in the sample reference system
      expanse33 = matmul(matmul(gmatinv,thermat),transpose(gmatinv))
      expanse33 = expanse33*ytemprate*dtime !dstrain = alpha*dT

      CALL kmatvec6(expanse33,dstranth)
      dstranth(4:6) = 2.0*dstranth(4:6)


C     *** DETERMINE INCREMENT IN TOTAL STRAIN (6X1) ***

      tempstrain=(L+transpose(L))*0.5*dtime
      spin=(L-transpose(L))*0.5

      CALL kmatvec6(tempstrain,dtotstran)
      dtotstran(4:6) = 2.0*dtotstran(4:6)


C     *** COMPUTE TRIAL STRESS ***

      DO i=1,6
          stressvec(i) = xstressdef(i) ! old stress
      END DO

      trialstress = stressvec+matmul(xStiffdef,dtotstran)-matmul(xStiffdef,dstranth)            
      CALL kvecmat6(trialstress,trialstressM)

      CALL kvecmat6(stressvec,stressM)
      trialstressM = trialstressM + (matmul(spin,stressM) - 
     + matmul(stressM,spin))*dtime 


C     *** CALCULATE RESOLVED SHEAR STRESS ON SLIP AND TWIN SYSTEMS  ***
C     *** NOW CALCULATED USING THE STRESS OF THE PREVIOUS INCREMENT ***
C     *** AS TRIAL STRESS				            ***

      DO I=1,nSys
        tempNorm = xNorm(I,:); tempDir = xDir(I,:)
        prod = matmul(stressM,tempNorm)
        tau(I)= dot_product(prod,tempDir)
        signtau(I) = 1.d0      
        IF(tau(I) .LT. 0.0) THEN
          tau(I) = -1.E0*tau(I) ! always positive RSS
          signtau(I) = -1.d0 ! carry info about the sign
        END IF
      END DO

      if (twinon == 1) then ! twin active
        DO I=nTwinStart,nTwinEnd ! here only homogenized stress is considered
          tempNorm = xTwinNorm(I,:); tempDir = xTwinDir(I,:)
          prod = matmul(stressM,tempNorm)
          tautwin(I) = dot_product(prod,tempDir)
          signtautwin(I) = 1.d0
          IF(tautwin(I) .LT. 0.0) THEN
            tautwin(I) = 0.0 ! twinning cannot be induced with negative RSS (no detwinning)
          END IF
        END DO
      end if ! twin active

      xtau = maxval(tau/tauc)
      xtautwin = maxval(tautwin/tauctwin)

C     *** PLASTIC DEFORMATION ***

      ! decide if Newton Raphson loop starts
      EnterNRLoop = 0
      if (xtau >= 0.5 .or. xtautwin >= 0.5) then ! stress condition
        EnterNRLoop = 1
      else
        ! twinvolfrac > 0.5 is needed for twin completion
        ! in the discrete twin model
        if (twinon == 1) then
          if (twinvolfrac(nTwinStart) > 0.5 .or. twinvolfrac(nTwinEnd) > 0.5) then
            EnterNRLoop = 1
          end if
        end if
      end if ! stress condition


      IF (EnterNRLoop == 1) THEN

      do while (debug == 1 .and. kint == 1)           
        debugwait = 0
      end do


      faivalue=1.
      iter=0
      fpold = Fp

C     *** USE NEWTON METHOD TO DETERMINE STRESS INCREMENT ***

      DO WHILE (faivalue .gt. xacc)
      
      iter=iter+1

      !============================================================================   
      !  Slip rule:
      !  Returns Lp and tmat required to define the material jacobian.
      !============================================================================

      IF (kslip == 0) THEN ! Original slip rule with no GND coupling i.e. using alpha and beta

      CALL kslip0(xNorm,xDir,tau,tauc,caratio,dtime,nSys,0.0,iphase,Lp,tmat) 
     
      ELSE IF (kslip == 5) THEN ! Original slip rule with GND coupling 
                 
      CALL kslip5ET(xNorm,xDir,tau,signtau,tauc,burgerv,rhossd,gndtot,
     +                gndall,gndcut,gndmob,dtime,nSys,ns,iphase,Lp,tmat)

      ELSE IF (kslip == 6) THEN ! updated slip rule with GND coupling
      
      CALL  kslip6ET(xNorm,xDir,tau,signtau,tauc,burgerv,dtime,nSys,
     +        iphase,irradiate,gndcut,gndtot,rhossd,Lp,tmat,gammaDot)

      ELSE IF (kslip == 7) THEN ! Powerlaw plasticity, orthorombic alpha-Uranium
      
      CALL  kslipPowerLaw(xNorm,xDir,xTwinNorm,xTwinDir,
     +        tau,tautwin,signtau,signtautwin,
     +        tauc,tauctwin,burgerv,dtime,nSys,
     +        nTwin,iphase,irradiate,gndcut,gndtot,
     +        rhossd,twinvolfrac,twinvolfractotal,
     +        Lp,tmat,gammaDot,gammatwindot,twinon,
     +        nTwinStart,nTwinEnd)

      END IF


      if(any(tmat /= tmat)) then

          call Mutexlock( 10 )   ! lock Mutex #5
         
          pnewdt = 0.5 ! if sinh( ) has probably blown up then try again with smaller dt
          write(*,*) "*** WARNING tmat  = NaN: jelem, kint, time: ", jelem, kint, time
          
          do while (debug == 1)
            debugWait = 1! wait here to attach debugger
          end do

          call MutexUnlock( 10 )   ! unlock Mutex #5
               
          return
      end if   

      !============================================================================

C     *** DETERMINE PLASTIC STRAIN INCREMENTS

      plasStrainInc = (Lp+transpose(Lp))*0.5*dtime
      CALL kmatvec6(plasStrainInc,plasStrainInc2)
      plasStrainInc2(4:6) = 2.0*plasStrainInc2(4:6)            


C     *** CALCULATE THE STRESS INCREMENT ***

      xjfai = xIden6 + matmul(xStiffdef,tmat) ! Jacobian of the Newton loop (see Dunne, Rugg, Walker, 2007)
      
      CALL lapinverse(xjfai,6,info3,xjfaiinv) ! invert Jacobian

      IF(info3 /= 0) write(6,*) "inverse failure: xjfai in kmat"
      
      ! trial stress here is in the assumption of full elastic increment
      ! C (epsilon_total - epsilon_plastic) = C ( epsilon_elastic ) = stress
      ! therefore the algorithm is finding the zero of "fai"
      ! and the variable that is updated at each iteration is "stressvec"

      fai = trialstress - stressvec - matmul(xStiffdef,plasStrainInc2)
      dstressinc = matmul(xjfaiinv,fai)

      stressvec = stressvec + dstressinc
      CALL kvecmat6(stressvec,stressM)
      faivalue = sqrt(sum(fai*fai))


C     *** UPDATE RESOLVED SHEAR STRESS ACCORDING TO NEW STRESS ***

      DO I=1,nSys
      
          tempNorm = xNorm(I,:); tempDir = xDir(I,:)    
          prod = matmul(stressM,tempNorm)
          tau(I)= dot_product(prod,tempDir)
          signtau(I) = 1.d0   
          IF(tau(I) < 0.0) THEN
            tau(I) = -1.E0*tau(I)         
            signtau(I) = -1.d0
          END IF

      END DO

      if (twinon == 1) then ! twin active
        DO I=nTwinStart,nTwinEnd ! here only homogenized stress is considered

          tempNorm = xTwinNorm(I,:); tempDir = xTwinDir(I,:)
          prod = matmul(stressM,tempNorm)
          tautwin(I) = dot_product(prod,tempDir)
          signtautwin(I) = 1.d0
          IF(tautwin(I) .LT. 0.0) THEN
            tautwin(I) = 0.0 ! twinning cannot be induced with negative RSS (no detwinning)
          END IF

        END DO
      end if ! twin active

      xtau = maxval(tau/tauc) 
      xtautwin = maxval(tautwin/tauctwin)

      IF (iter .gt. maxNRiter) THEN

          call Mutexlock( 11 )   ! unlock Mutex
          
          pnewdt = 0.5
          WRITE(*,*) "WARNING NEWTON LOOP NOT CONVERGED: jelem, kint, time:", jelem, kint, time(1)
          WRITE(*,*) "fai", fai
          WRITE(*,*) "stressM", stressM

          call MutexUnlock( 11 )   ! unlock Mutex

          return

      END IF

      !!*** THE END OF NEWTON ITERATION ***
      END DO
	
          
C     *** NOW CALCULATE THE JACOBIAN***

      C = matmul(xjfaiinv,xStiffdef)

      ! assign the updated stress
      xstressmdef = stressM


C     *** UPDATE OUTPUT VARIABLES ***

      plasStrainrate=(Lp+transpose(Lp))*0.5 
      
      
C     *** UPDATE PLASTIC DEFORMATION GRADIENT
    
      print2 = 0.; print3 = 0.
      print2 = xI - Lp*dtime      
      CALL kdeter(print2,deter)      
      IF (deter /= 0.0) THEN
         CALL lapinverse(print2,3,info4,print3)
         Fp = matmul(print3,fpold)
      ELSE
         Fp = Fp
      END IF
      
      
    !=========================================================================
C     *** DETERMINE DENSITY OF GNDs
    ! Definitely needs to be inside plasticity loop
    ! And should use the directions that have been modified according to tau!
    !=========================================================================
    ! Edge-screw separated
      
      
      IF (gndon == 0) THEN !Switching GND evolution on and off!
       gndall = 0.
       gndcut = 0.
            
      ELSE
      
          CALL kgndl2ET(curlfp,xNorm,xDir,tau,burgerv,iphase,nSys,ns,
     +       screwplanes,jelem,kint,time,gndall,gndcut,gndmob)
          
          IF (xtau < 1.0) THEN  
              write(*,*) "xtau<1, jelem,kint,time",jelem,kint,time
          END IF
      END IF

      ! sum all GNDs
      gndtot = sum(gndall)

      ! update state variables that determine hardening
      call kHardening(pdot,p,plasStrainrate,dtime,slip,gammaDot,
     + nSys,irradiate,tauSolute,gammast,rhossd,iphase,rhofor,
     + CurrentTemperature,rhosub,twinon,twinvolfrac,nTwin,
     + nTwinStart,nTwinEnd,gammaTwinDot)

    !=========================================================================

C     *** ELASTIC DEFORMATION ***     
      ELSE
          xstressmdef = trialstressM
          C = xStiffdef      
      END IF ! end of PLASTIC DEFORMATION
      
      CALL kmatvec6(xstressmdef,xstressdef) !output stress
      devstress = xstressmdef - 1./3.*trace(xstressmdef)*xI
      vms = sqrt(3./2.*(sum(devstress*devstress))) !von mises stress 


C     *** ORIENTATION UPDATE ***
      ! Assuming that all rigid body rotation is lumped into Fe and that the elastic strains are small 
      ! then the elastic spin is the antisymmetric part of We = d(Fe)/dt inv(Fe)
      ! L = We + Fe Lp inv(Fe) therefore 
      ! We = L - Fe Lp inv(Fe)
      ! G(t+dt) = G(t) + We G(t)dt dt or an implicit update is G(t+dt) = G(t)exp[We(t+dt)dt]  ~ inv[I - We(t+dt) dt] G(t) 
      
      ! We need Fe and inv(Fe) using F = Fe Fp gives Fe = F.inv(Fp)
      CALL kdeter(Fp,deter)      

      IF (deter /= 0.) THEN
         Fpinv = 0.
         CALL lapinverse(Fp,3,info5,Fpinv)
!        IF(info5 /= 0) write(6,*) "inverse failure: print3 in kmat"
         Fe = matmul(F,Fpinv)
      ELSE
         write(*,*) "Error in orientation update: finding inv(Fp)",jelem,kint, kinc
         write(*,*) "Fp, det(Fp)", Fp, deter
         write(*,*) "fpold", fpold
         write(*,*) "Lp", Lp
         Fp = fpold
         Fpinv = 0.
         CALL lapinverse(Fp,3,info5,Fpinv)
         Fe = matmul(F,Fpinv)
         
         do while (debug == 1)
           debugwait = 1
         end do

      END IF
      
      
      CALL kdeter(Fe,deter)      
      
      IF (deter /= 0.) THEN
         Feinv = 0.
         CALL lapinverse(Fe,3,info5,Feinv)
!         IF(info5 /= 0) write(6,*) "inverse failure: print3 in kmat"         
      ELSE
          write(*,*) "Error in orientation update: finding inv(Fe)",jelem,kint, kinc
		  write(*,*) "Fe", Fe
         call XIT 
      
      END IF
            
      LpFeinv = 0.; 
      LpFeinv = matmul(Lp, Feinv)
      Le = L - matmul(Fe,LpFeinv)        
      elasspin=(Le-transpose(Le))*0.5
      matrix = xI - elasspin*dtime      
      CALL kdeter(matrix,deter)      
      

C     *** if plastic deformation took place
C     *** use the elastic spin to rotate the corotational
C     *** stress tensor
      if (iter > 0) then
        print2 = (matmul(elasspin,stressM)-matmul(stressM,elasspin))
        xstressmdef = xstressmdef + print2*dtime 
        CALL kmatvec6(xstressmdef,xstressdef) !output stress
      end if           
      
      
      IF (deter /= 0.) THEN
         update = 0.
         CALL lapinverse(matrix,3,info5,update)
         IF(info5 /= 0) write(*,*) "inverse failure: print3 in kmat"
         gmatinvnew = matmul(update,gmatinvold)                  
       
      ELSE         
         gmatinvnew = gmatinvold
         write(*,*) "WARNING gmatinv not updated at jelem,kint, kinc:", jelem,kint, kinc
      END IF      

      gmatinv = gmatinvnew            
      
      if (maxval(gmatinv) > 1) then
        !write(*,*) "something very wrong with gmatinv"
        !write(*,*) "jelem,kint,kinc",jelem,kint,kinc
        !write(*,*) "maxval(gmatinv)",maxval(gmatinv)
        !call XIT
      end if

C     *** CALCULATE GREEN_LAGRANGE STRAIN FOR OUTPUT ***
C     ***         ROTATE: LATTICE COORDINATES        ***

      EECrys = matmul(transpose(Fe),Fe)
      EECrys = EECrys - xI
      EECrys = 0.5*EECrys	
      EECrys = matmul(matmul(transpose(gmatinv),EECrys),gmatinv)
             

C     *** UPDATE STATE VARIABLES *** !

      DO i=1,3
        DO j=1,3
          usvars(j+(i-1)*3) = gmatinv(i,j)
        END DO
      END DO

      usvars(10) = p
 
      DO i=1,6
        usvars(10+i) = totplasstran(i) + plasStrainInc2(i)
      END DO

      DO i=1,6
        usvars(16+i) = totstran(i) + dtotstran(i)
      END DO
 
      usvars(26) = gndtot   
      
      DO i=1,6
       usvars(47+i) = xstressdef(i)
      END DO

      usvars(32) = maxval(plasStrainrate)
      usvars(33) = pdot
      usvars(34) = xtau
      usvars(35) = slip
      usvars(36) = tauSolute
      usvars(54) = rhossd
      usvars(55) = vms
      usvars(56) = maxval(tauc)
     
      !GNDs on indiviual systems
      !max(nSys) is currently limited to 24. IF all 48 of bcc is needed, storage should be raised to match that!
      DO i=1,nSys
       usvars(56+i) = gndcut(i)
      END DO
      
      usvars(71) = Temperature

      DO i=1,3
       DO j=1,3 
        usvars(71+j+(i-1)*3) = EECrys(i,j)
        usvars(80+j+(i-1)*3) = Fp(i,j)
       END DO
      END DO   

      if (iphase == 5) then ! orthorombic alpha-Uranium model

      ! rewrite gndcut state variables
      ! for the alpha-Uranium model
        DO i=1,nSys
          usvars(56+i) = rhofor(i)
        END DO
        usvars(65) = rhosub
        DO i=1,nSys
	  usvars(98+i) = tauc(i)
        END DO
        ! output some of the plastic strain rates on the slip systems
        usvars(90) = gammaDot(1)
        usvars(91) = gammaDot(2)
        usvars(92) = gammaDot(3)
        usvars(93) = gammaDot(4)
        usvars(94) = gammatwindot(nTwinStart)
        usvars(95) = gammatwindot(nTwinEnd)
        usvars(96) = TwinIntegral(nTwinStart)
        usvars(97) = TwinIntegral(nTwinEnd)
        DO i=1,nTwin ! twin volume fractions
          usvars(106+i) = twinvolfrac(i)
        END DO
        DO i=1,nTwin ! twin CRSS
          usvars(112+i) = tauctwin(i)
        END DO

      end if

      RETURN

      END
