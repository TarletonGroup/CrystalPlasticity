C   Christos Skamniotis
C   University of Oxford
C   December 2021 
C
C   Simplified Double exponent slip rule and empirical creep law for tertiary creep:

      subroutine NickelSuperalloy(xNorm,xDir,tau,signtau,tauc,
     + burgerv,dtime,nSys,iphase,CurrentTemperature,Lp,
     + tmat,gammaDot, cubicslip,creep, usvars, nsvars)
      
      implicit none
	  
	  ! number of slip system
      integer, intent(in):: nSys
        ! activation flag for cubic slip (additional 6 systems activated when loading is along 111)
      INTEGER,intent(in) :: cubicslip
        ! activation flag for tertiary creep
      INTEGER,intent(in) :: creep

	  ! phase
      integer, intent(in):: iphase
      
       ! number of Abaqus state variables
      INTEGER,intent(in) :: nsvars

      ! slip directions and normals	  
      real*8, intent(in) :: xNorm(nSys,3),xDir(nSys,3)
	  
	  ! resolved shear stress and critical resolved shear stress
	  ! and sign of the resolved shear stress
	  ! tauc is positive by definition
      real*8, intent(in) :: tau(nSys), tauc(nSys), signtau(nSys)
	  
	  ! Burgers vectors
      real*8, intent(in) :: burgerv(nSys)
	  
	  ! time step
      real*8, intent(in) :: dtime	  
	  
	  ! Temperature in Kelvin
	  real*8, intent(in) :: CurrentTemperature
        
      ! Abaqus state variables
        REAL*8,intent(in) :: usvars(nsvars)
	  
	  ! plastic velocity gradient
      real*8, intent(out) :: Lp(3,3)
	  
	  ! and its derivative with respect to the stress
	  real*8, intent(out) :: tmat(6,6)
	  
	  ! plastic strain rate on each slip system
	  real*8, intent(out) :: gammaDot(nSys)
        
        ! Boltzmann constant (J/K)
	  real*8, parameter :: kB = 1.38e-23
        
        ! Gas constant (J*mol/K)
	  real*8, parameter :: R = 8.314462
      
        
******************************************
** The following parameters must be set **
c
*** RATE DEPENDENT PLASTICITY (thermally activated glide)

      ! Activation energy for octahedral slip  (J)
        real*8, parameter :: Foctahedral = 9.39e-19
      ! Activation energy for cubic slip  (J)
        real*8, parameter :: Fcubic = 1.17e-18
      ! reference strain rate (1/s)
	  real*8, parameter :: gammadot0 = 1.0e7    ! real 1.0e-7
      ! rate sensitivity exponents 
	  real*8, parameter :: p = 0.78      ! real 0.78
        real*8, parameter :: q = 1.15      ! real 1.15
c
*** TERTIARY CREEP (dislocation climb & damage)
   !!!!  Initial creep rate constants        
      ! Activation energy for creep (J/mol)
      real*8, parameter :: Qo = 460000.0
      ! reference rate (1/s)
      real*8, parameter :: ao = 4.0e8 
      ! stress multiplier (1/MPa)
      real*8, parameter :: bo = 3.2e-2
   !!!!  Climb/damage constants  
      ! Activation energy for damage (J/mol)
      real*8, parameter :: QD = 340000.0
      ! reference rate (1/s)
      real*8, parameter :: aD = 6000000.0
      ! stress multiplier (1/MPa)
      real*8, parameter :: bD = -5.0e-08
   !!!!  Rafting   
C      real*8, parameter :: SS = 100
C      real*8, parameter :: TT = 1000 
C      real*8, parameter :: QQ = 20000
C      real*8, parameter :: m = -3
	  
**       End of parameters to set       **
******************************************


        ! slip system index
        integer :: i
	  
	  ! Schmid tensor and its transpose
	  real*8 :: SNij(3,3), NSij(3,3)
	  
	  ! Schmid tensor and its transpose in Voigt notation
	  real*8 :: sni(6), nsi(6)
	  
	  ! higher order Schmid tensor in Voigt notation
	  real*8 :: SNNS(6,6)
	  
	  ! temporary slip normal and slip direction
	  real*8 :: tempNorm(3), tempDir(3)
	  
	  ! temporary variable to calculate the Jacobian
	  real*8 :: result1
	 
	  ! Jacobian
        real*8 :: result4(6,6)
        
        ! activation energy
	  real*8 :: dF
        
        ! RSS/CRSS ratio
	  real*8 :: tau_ratio
      
C
C
C  *** CALCULATE LP AND THE DERIVATIVE OF PLASTIC STRAIN INCREMENT WITH 
C   RESPECT TO THE STRESS DEFINED AS tmat***
C
C
        tmat = 0.0
	  Lp = 0.0
	  result4 = 0.0
	  
	  ! contribution to Lp of all slip systems
      do i=1,nSys
	  
        tau_ratio=tau(i)/tauc(i)

       if  (tau_ratio >= 0.0) then
           
        if (tau_ratio >= 1.0) then ! avoid negative values before elevating to power q
            
            gammaDot(i) = signtau(i)*gammadot0  !*exp(-dF/(kB*CurrentTemperature))
            
        else ! standard case
            
            dF=Foctahedral
            
            if (i .gt. 12) then
               dF=Fcubic
            end if   
            
      ! strain rate due to thermally activated glide (rate dependent plasticity)
      gammaDot(i) = signtau(i)*gammadot0*
     + exp(-(dF/(kB*CurrentTemperature))*(1- tau_ratio**p)**q)
            
        end if
		 
        !  add tertiary creep strain rate
        if (creep == 1) then  
           
      gammaDot(i) = gammaDot(i) +
     + signtau(i)*ao*exp(bo*tau(i) -
     + Qo/(R*CurrentTemperature)) +
     + signtau(i)*abs(usvars(89+i))*aD*exp(bD*tau(i) -
     + QD/(R*CurrentTemperature))
                            
        end if
        
        
          tempNorm = xNorm(i,:)
          tempDir = xDir(i,:)
          SNij = spread(tempDir,2,3)*spread(tempNorm,1,3)
          NSij = spread(tempNorm,2,3)*spread(tempDir,1,3)
          call KGMATVEC6(SNij,sni)         
          call KGMATVEC6(NSij,nsi) 
          SNNS = spread(sni,2,6)*spread(nsi,1,6)
		  
          ! calculate derivative d ( gammaDot(i) ) / d ( tau(i) )
          result1 = abs(gammaDot(i))
          result1 = result1 * dF/(kB*CurrentTemperature)
          result1 = result1 * q
          result1 = result1 * (1- tau_ratio**p)**(q-1.0)
          result1 = result1 * p
          result1 = result1 / tauc(i)
          result1 = result1 * tau_ratio**(p-1.0)
          
          if (creep == 1) then
          
      result1 = result1 + ao*bo*
     + exp(bo*tau(i)-Qo/(R*CurrentTemperature)) +
     + abs(usvars(89+i))*aD*bD*
     + exp(bD*tau(i)-QD/(R*CurrentTemperature))
              
          end if   
          
		  ! contribution to Jacobian
          result4 = result4 + dtime*result1*SNNS     
		  
		  ! plastic velocity gradient contribution
          Lp = Lp + gammaDot(i)*SNij 
		
        else   
          gammaDot(i) = 0.0
		  
        end if

      end do
      
      tmat = 0.5*(result4+transpose(result4))
             
      return
      end 

