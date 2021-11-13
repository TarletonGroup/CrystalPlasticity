C   Nicolo Grilli
C   University of Bristol
C   Christos Skamniotis
C   University of Oxford
C   11 Novembre 2021 
C
C   Double exponent slip rule in:
C   Zhengxuan Fan & Serge Kruch (2020) A comparison of different crystal
C   plasticity finite-element models on the simulation of nickel alloys, 
C   Materials at High Temperatures, 37:5, 328-339, DOI: 10.1080/09603409.2020.1801951

      subroutine kslipDoubleExponent(xNorm,xDir,tau,signtau,tauc,burgerv,dtime,
     + nSys,iphase,CurrentTemperature,Backstress,Lp,tmat,gammaDot)
         
      implicit none
	  
	  ! number of slip system
      integer, intent(in):: nSys
	  
	  ! phase
	  integer, intent(in):: iphase

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
	  
	  ! Temperature
	  real*8, intent(in) :: CurrentTemperature

	  ! Backstress state variable
	  real*8, intent(in) :: Backstress(nSys)
	  
	  ! plastic velocity gradient
      real*8, intent(out) :: Lp(3,3)
	  
	  ! and its derivative with respect to the stress
	  real*8, intent(out) :: tmat(6,6)
	  
	  ! plastic strain rate on each slip system
	  real*8, intent(out) :: gammaDot(nSys)
	 

******************************************
** The following parameters must be set **

      ! Free energy variation (J/mol)
      real*8, parameter :: dF = 286000.0
	  
	  ! Boltzmann constant (J/K)
	  real*8, parameter :: kB = 1.38e-23
	  
	  ! reference strain rate (1/s)
	  real*8, parameter :: gammadot0 = 1.0e7
	  
	  ! ratio between shear modulus at CurrentTemperature
	  ! and at 0K
	  real*8, parameter :: mu_over_mu0 = 103.5 / 134.0
	  
	  ! strain rate sensitivity exponents 
	  real*8, parameter :: p = 0.59
      real*8, parameter :: q = 1.8	  
	  
	  ! Peierls stress at 0K (MPa)
	  real*8, parameter :: tau0 = 342.0
	  
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
	  
	  ! ratio between free energy jump and KB * T
	  real*8 :: dF_over_kBT = dF / (kB * CurrentTemperature)
	  
	  ! argument of the exponential to calculate gammaDot
	  real*8 :: gammaDot_exp_arg
	  
	  ! effective stress for slip
	  real*8 :: tau_eff
	  
	  
	  
	  
	  ! product between tauc and mu_over_mu0
	  real*8 :: tauc_mu_over_mu0
	  
	  
      	  
      real*8 :: alpha,beta,result1, rhom,,f,T,k,gamma0,b,psi,V,
     + result4(6,6),
     + 
	 
      	 
      
C
C  *** CALCULATE LP AND THE DERIVATIVE OF PLASTIC STRAIN INCREMENT WITH 
C   RESPECT TO THE STRESS DEFINED AS tmat***
C    

      tmat = 0.0
	  Lp = 0.0
	  result4 = 0.0
	  
	  ! contribution to Lp of all slip systems
      do i=1,nSys
	  
	    tau_eff = abs(tau(i) - Backstress(i)) - mu_over_mu0 * tauc(i)

        if (tau_eff >= 0.0) then

          gammaDot_exp_arg = tau_eff / (mu_over_mu0 * tau0)
		  gammaDot_exp_arg = 1.0 - gammaDot_exp_arg**p
		  gammaDot_exp_arg = gammaDot_exp_arg**q
                                                 
          gammaDot(i) = gammadot0*exp(-dF_over_kBT*gammaDot_exp_arg)*signtau(i)
         
          tempNorm = xNorm(i,:)
          tempDir = xDir(i,:)
          SNij = spread(tempDir,2,3)*spread(tempNorm,1,3)
          NSij = spread(tempNorm,2,3)*spread(tempDir,1,3)
          call KGMATVEC6(SNij,sni)         
          call KGMATVEC6(NSij,nsi) 
          SNNS = spread(sni,2,6)*spread(nsi,1,6)
          result1 = cosh(beta*signtau(i)*(tau(i) - tauc(i)))
          
          result4 = result4 + alpha*beta*dtime*result1*SNNS          
          Lp = Lp + gammaDot(i)*SNij
		  
        else
		
          gammaDot(i)=0.0
		  
        end if

      end do
      
      tmat = 0.5*(result4+transpose(result4))
             
      return
      end 