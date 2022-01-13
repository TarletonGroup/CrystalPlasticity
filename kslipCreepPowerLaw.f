C   Christos Skamniotis
C   University of Oxford
C   December 2021 
C
C   Simplified power law slip rule and empirical creep law for tertiary creep:

      subroutine kslipCreepPowerLaw(xNorm,xDir,tau,signtau,tauc,
     + dtime,nSys,iphase,CurrentTemperature,Lp,
     + tmat,gammaDot,cubicslip,creep,slipsysplasstran)
      
      implicit none
	  
	  ! number of slip system
      integer, intent(in):: nSys
	  
      ! activation flag for cubic slip (additional 6 systems activated when loading is not along 001)
      INTEGER,intent(in) :: cubicslip
	  
      ! activation flag for tertiary creep
      INTEGER,intent(in) :: creep
      
	  ! accumulated plastic strain on each slip system, signed
      real*8, intent(in) :: slipsysplasstran(nSys)

	  ! phase
      integer, intent(in):: iphase

      ! slip directions and normals	  
      real*8, intent(in) :: xNorm(nSys,3), xDir(nSys,3)
	  
	  ! resolved shear stress and critical resolved shear stress
	  ! and sign of the resolved shear stress
	  ! tauc is positive by definition
      real*8, intent(in) :: tau(nSys), tauc(nSys), signtau(nSys)
	  
	  ! time step
      real*8, intent(in) :: dtime	  
	  
	  ! Temperature in Kelvin 
	  real*8, intent(in) :: CurrentTemperature
	  
	  ! plastic velocity gradient
      real*8, intent(out) :: Lp(3,3)
	  
	  ! and its derivative with respect to the stress
	  real*8, intent(out) :: tmat(6,6)
	  
	  ! plastic strain rate on each slip system
	  real*8, intent(out) :: gammaDot(nSys)
        
      ! Gas constant (J*mol/K)
	  real*8, parameter :: R = 8.314462

******************************************
** The following parameters must be set **

*** RATE DEPENDENT PLASTICITY (thermally activated glide)

      ! reference strain rate (1/s)
	  real*8, parameter :: ref_gammaDot = 7.071136E-05
	  
      ! rate sensitivity multiplier (1/Kelvin)
	  real*8, parameter :: slopeM = -0.036273
      
	  ! rate sensitivity constant (-/-)
      real*8, parameter :: constantM = 65.33827694

*** TERTIARY CREEP (dislocation climb & damage)

      ! Initial creep rate constants        
      ! Activation energy for creep (J/mol)
      real*8, parameter :: Qo = 440.0e3
	  
      ! reference rate (1/s)
      real*8, parameter :: ao = 6.0e7
	  
      ! stress multiplier (1/MPa)
      real*8, parameter :: bo = 5.0e-2

      ! Climb/damage constants  
      ! Activation energy for damage (J/mol)
      real*8, parameter :: QD = 170.0e3

      ! reference rate (1/s)
      real*8, parameter :: aD = 6.0e-1

      ! stress multiplier (1/MPa)
      real*8, parameter :: bD = 3.5e-2
	  
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
        
      ! RSS/CRSS ratio
      real*8 :: tau_ratio

	  ! rate sensitivity
	  real*8 :: mpower

C
C  *** CALCULATE LP AND THE DERIVATIVE OF PLASTIC STRAIN INCREMENT WITH 
C   RESPECT TO THE STRESS DEFINED AS tmat***
C
C
      tmat = 0.0
      Lp = 0.0
      result4 = 0.0
      mpower = slopeM * CurrentTemperature + constantM
        
	  ! contribution to Lp of all slip systems
      do i=1,nSys
	  
        tau_ratio=tau(i)/tauc(i)

        if (tau_ratio > 0.0) then
            
      ! strain rate due to thermally activated glide (rate dependent plasticity)
      gammaDot(i) = signtau(i)*ref_gammaDot*tau_ratio**mpower
            
      ! calculate derivative d ( gammaDot(i) ) / d ( tau(i) )
      result1 = ref_gammaDot*mpower*(1/tauc(i))*(tau_ratio**(mpower-1))
          
        !  add tertiary creep rate
        if (creep == 1 .and. tau_ratio > 0.05) then  
           
      gammaDot(i) = gammaDot(i) +
     + signtau(i)*ao*exp(bo*tau(i)-Qo/(R*CurrentTemperature)) +
     + signtau(i)*abs(slipsysplasstran(i))*aD*
     + exp(bD*tau(i)-QD/(R*CurrentTemperature))
          
      result1 = result1 + ao*bo*
     + exp(bo*tau(i)-Qo/(R*CurrentTemperature)) +
     + abs(slipsysplasstran(i))*aD*bD*
     + exp(bD*tau(i)-QD/(R*CurrentTemperature))
          
        end if
        
        ! calculate SNNS
          tempNorm = xNorm(i,:)
          tempDir = xDir(i,:)
          SNij = spread(tempDir,2,3)*spread(tempNorm,1,3)
          NSij = spread(tempNorm,2,3)*spread(tempDir,1,3)
          call KGMATVEC6(SNij,sni)         
          call KGMATVEC6(NSij,nsi) 
          SNNS = spread(sni,2,6)*spread(nsi,1,6)
  
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
