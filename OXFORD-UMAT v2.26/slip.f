!     Oct. 6th, 2022
!     Eralp Demir
!     Slip laws
!     1. sinh law
!     2. double exponent law
!     3. power law
      module slip
      implicit none
      contains
!
!
      subroutine sinhslip(Schmid_0,
     + Schmid,SchmidxSchmid,
     + tau,X,tauc,rhofor,
     + burgerv,dt,nslip,iphase,T,slipparam,
     + irradiationmodel,irradiationparam,
     + cubicslip,caratio,
     + Lp,Dp,Pmat,gammadot,
     + dgammadot_dtau,
     + dgammadot_dtauc)
!
      use userinputs, only : useaveragestatevars, maxnparam
      use globalvariables, only : KB
      use utilities,  only : gmatvec6
      use errors, only: error
      implicit none
!     Schmid tensor at the intermediate config. (or undeformed)
      real(8), intent(in) :: Schmid_0(nslip,3,3)
!     Schmid tensor
      real(8), intent(in) :: Schmid(nslip,3,3)
!     Schmid dyadic
      real(8), intent(in) :: SchmidxSchmid(nslip,6,6)
!     RSS
      real(8), intent(in) :: tau(nslip)
!     Backstress
      real(8), intent(in) :: X(nslip)
!     CRSS
      real(8), intent(in) :: tauc(nslip)
!     Forest dislocation density
      real(8), intent(in) :: rhofor(nslip)
!     Burger's vector
      real(8), intent(in) :: burgerv(nslip)
!     time increment
      real(8), intent(in) :: dt
!     number of slip systems
      integer, intent(in) :: nslip
!     phase id
      integer, intent(in) :: iphase
!     temperature
      real(8), intent(in) :: T
!     slip parameters
      real(8), intent(in) :: slipparam(maxnparam)
!     irradiation model id
      integer, intent(in) :: irradiationmodel
!     irradiation model parameters
      real(8), intent(in) :: irradiationparam(maxnparam)
!     cubic slip for fcc
      integer, intent(in) :: cubicslip
!     c/a ratio for hcp materials
      real(8), intent(in) :: caratio
!     plastic part of velocity gradient
      real(8), intent(out) :: Lp(3,3)
!     plastic stretch rate at the deformed configuration
      real(8), intent(out) :: Dp(3,3)
!     tangent matrix required for N-R iteration (at the inner loop)
      real(8), intent(out) :: Pmat(6,6)
!     slip rates
      real(8), intent(out) :: gammadot(nslip)
!     Derivative of slip rates wrto rss
      real(8), intent(out) :: dgammadot_dtau(nslip)
!     Derivative of slip rates wrto crss
      real(8), intent(out) :: dgammadot_dtauc(nslip)
!     variables used within this subroutine       
      integer :: is
      real(8)  :: alpha0, alpha, beta0, beta, rhom, rhom0,
     + DeltaF, nu0, gamma0, AV0, psi, lambda, AV, rhoav
      real(8) :: abstau, signtau
!
!
!
      gammadot = 0.
      dgammadot_dtau = 0.
      dgammadot_dtauc = 0.
!
!
!
!
!     Obtain slip parameters
!     constant alpha
      alpha0 = slipparam(1)
!     constant beta
      beta0 = slipparam(2)
!     rhom0 - mobile dislocation density
      rhom0 = slipparam(4)
!     DeltaF - activation energy to overcome Pierls barrier
      DeltaF = slipparam(5)
!     nu0 - attempt frequency
      nu0 = slipparam(6)
!     gamma0 - multiplier for activation volume
!     Unit conversion factor is for Boltmann constant (J/K) and stress (MPa)
!     The multiplier 1d-12 stands for unit conversion for beta
      gamma0 = slipparam(7)*1.d-12
!     AV0 - activation volume (if defined not=0)
!     Unit conversion factor is for Boltmann constant (J/K) and stress (MPa)
!     The multiplier 1d-12 stands for unit conversion for beta
      AV0 = slipparam(8)*1.d-12
!
!
!
!     if alpha is defined parametrically
      if (alpha0 == 0.) then
!
!         psi - fraction of mobile dislocations
!         If there is no irradiation
          if (irradiationmodel == 0) then
!
              psi = slipparam(3)
!
          elseif (irradiationmodel == 1) then
!
              psi = irradiationparam(3)
!
          elseif (irradiationmodel == 2) then
!
              psi = slipparam(3)
!
          endif
!
!
!         Scale mobile dislocation density
          rhom = psi * rhom0
!
!
      endif
!
!
!
!
      Pmat=0.; Lp = 0.; Dp = 0.
!     Loop through slip systems
      do is = 1,nslip
!
!         Absolute value of RSS
          abstau = abs(tau(is)-X(is))
!
!         Sign of RSS
          signtau = sign(1.0,tau(is)-X(is))
!
!
!         Outer multipler "alpha" calculation:
!
!         if alpha is defined parametrically
          if (alpha0 == 0.) then
!
!
              alpha = rhom*burgerv(is)**2.*nu0*exp(-DeltaF/KB/T)
!
!         alpha is defined as a constant
          else
!
              alpha = alpha0
!
          endif
!
!
!
!         Activation volume calculation:
!
!         if AV is defined parametrically
          if (AV0 == 0.) then
!
!
              if (useaveragestatevars == 0) then
!
!                 average spacing based on forest densities
                  lambda = 1./sqrt(rhofor(is))
!
!                 activation volume
                  AV = gamma0*lambda*burgerv(is)**2.
!
!
!
              elseif (useaveragestatevars == 1) then
!
!                 average forest density
                  rhoav = sum(rhofor)/nslip
!
!                 average spacing
                  lambda = 1./sqrt(rhoav)
!
!                 activation volume
                  AV = gamma0*lambda*burgerv(is)**2.
!
!
!
              endif
!
!         AV is defined as a constant
          else
!
!
              AV = AV0 * burgerv(is)**3.
!
!
!
          end if
!
!
!         THESE CHECKS SHALL BE PLACED TO INITIALIZATION!
!!         Check for zero activation volume
!          if (AV == 0.) then
!              call error(8)
!          end if 
!
!
!
!
!         Inner multiplier "beta" calculation:
!
!         if beta is defined parametrically
          if (beta0 == 0.) then
!
!
!
              beta = AV/KB/T
!
!
!
!         beta is defined as a constant
          else
!
              beta = beta0
!
!
          endif
!
!
!
!         c/a ratio correction for HCP
          if(iphase == 3) then
!             Correction by C. Hardie - 26.05.2022
!             This correction needs to be done if activation volume
!             IS NOT DEFINED PARAMETRICALLY, because it is already taken into account then!
              if (AV0 /= 0.0) then
!                 Corrected by A. Pechero - 23.02.2023
!                 same burgers magnitude for 1st-12th slip systems
!                 changes only if slip system is greater than 12th
                  if (is > 12) then
                      alpha = caratio*caratio*alpha
                      beta = caratio*caratio*beta
                  endif
              end if
          end if
!
!         This is critical threshold
          if (abstau >= tauc(is)) then
!
!             slip rate
              gammadot(is) = alpha*
     + sinh(beta*(abstau-tauc(is)))*signtau
!
!
              dgammadot_dtau(is) = alpha*beta*
     + cosh(beta*(abstau-tauc(is)))
!
!
              dgammadot_dtauc(is) = -alpha*beta*
     + cosh(beta*(abstau-tauc(is)))*signtau
!
!
              Pmat = Pmat + dt*dgammadot_dtau(is)*
     +        SchmidxSchmid(is,:,:)
!
!
              Lp = Lp + gammadot(is)*Schmid_0(is,:,:)
!
!             Update plastic velocity gradient
              Dp = Dp + gammadot(is)*Schmid(is,:,:)
!
!
!
          end if
!
          end do
!
!     Find plastic flow direction
      Dp = (Dp + transpose(Dp))/2.          
!
!     Symmetrize Pmat
      Pmat = (Pmat + transpose(Pmat))/2.
!
!
      end subroutine sinhslip
!
!
!
!
!
!
!
!
!
!
!     Zhengxuan Fan & Serge Kruch (2020) A comparison of different crystal
!     plasticity finite-element models on the simulation of nickel alloys, 
!     Materials at High Temperatures, 37:5, 328-339, DOI: 10.1080/09603409.2020.1801951      
!
      subroutine doubleexpslip(Schmid_0,
     + Schmid,SchmidxSchmid,
     + tau,X,tauc,
     + burgerv,dt,nslip,iphase,T,slipparam,
     + irradiationmodel,irradiationparam,
     + cubicslip,caratio,
     + Lp,Dp,Pmat,gammadot,dgammadot_dtau,
     + dgammadot_dtauc)
!
      use userinputs, only : useaveragestatevars, maxnparam
      use globalvariables, only : KB, smallnum
      use utilities,  only : gmatvec6
      implicit none
!     Schmid tensor at the intermediate config. (or undeformed)
      real(8), intent(in) :: Schmid_0(nslip,3,3)
!     Schmid tensor
      real(8), intent(in) :: Schmid(nslip,3,3)
!     Schmid dyadic
      real(8), intent(in) :: SchmidxSchmid(nslip,6,6)
!     RSS
      real(8), intent(in) :: tau(nslip)
!     Backstress
      real(8), intent(in) :: X(nslip)
!     CRSS
      real(8), intent(in) :: tauc(nslip)
!     Burger's vector
      real(8), intent(in) :: burgerv(nslip)
!     time increment
      real(8), intent(in) :: dt
!     number of slip systems
      integer, intent(in) :: nslip
!     phase id
      integer, intent(in) :: iphase
!     temperature
      real(8), intent(in) :: T
!     slip parameters
      real(8), intent(in) :: slipparam(maxnparam)
!     irradiation model id
      integer, intent(in) :: irradiationmodel
!     irradiation model parameters
      real(8), intent(in) :: irradiationparam(maxnparam)
!     cubic slip for fcc
      integer, intent(in) :: cubicslip
!     c/a ratio for hcp materials
      real(8), intent(in) :: caratio
!     plastic part of velocity gradient
      real(8), intent(out) :: Lp(3,3)
!     plastic stretch rate at the deformed configuration
      real(8), intent(out) :: Dp(3,3)
!     tangent matrix required for N-R iteration (at the inner loop)
      real(8), intent(out) :: Pmat(6,6)
!     slip rates
      real(8), intent(out) :: gammadot(nslip)
!     Derivative of slip rates wrto rss
      real(8), intent(out) :: dgammadot_dtau(nslip)
!     Derivative of slip rates wrto rss
      real(8), intent(out) :: dgammadot_dtauc(nslip)
!
!     variables used within this subroutine       
      integer :: is
      real(8) :: gammadot0, p, q, Foct, Fcub, DeltaF, ratio
      real(8) :: abstau, signtau
!
!
!     Obtain slip parameters
!     Reference strain rate
      gammadot0 = slipparam(1)
!     Inner exponent (p)
      p = slipparam(2)
!     Outer exponent (q)
      q = slipparam(3)
!     Activation energy for octahedral slip (J)
      Foct = slipparam(4)
!     Activation energy for cubic slip (J)
      Fcub = slipparam(5)
!
!
!
      Pmat = 0.; Lp = 0.; Dp = 0.
      gammadot = 0.
      dgammadot_dtau = 0.
      dgammadot_dtauc = 0.
!
!     Contribution to Lp of all slip systems
      do is=1,nslip
!
!
!         Absolute value of RSS
          abstau = abs(tau(is)-X(is))
!
!         Sign of RSS
          signtau = sign(1.0,tau(is)-X(is))
!
!         RSS/CRSS ratio
          ratio=abstau/tauc(is)
!
!         Avoid negative values before elevating to power q
          if (ratio >= 1.0) then
!
              gammadot(is) = signtau*gammadot0 !*exp(-dF/(kB*CurrentTemperature))
!
!         Standard case
          elseif (ratio /= 0.) then
!
!             Activation energy
              DeltaF=Foct
!
!             Cubic slip case with a different activation energy
!             If cubic slip systems are defined
              if (cubicslip == 1) then
!                 For cubic slip systems: 13-..-18
                  if (is > 12) then
                      DeltaF=Fcub
                  end if
              endif
!
!
!             Strain rate due to thermally activated glide (rate dependent plasticity)
              gammadot(is) = signtau*gammadot0*
     + exp(-(DeltaF/KB/T)*((1.-(ratio**p))**q))
!
              if (abs(gammadot(is)) > smallnum) then
!
!                 Calculate derivative d ( gammadot(i) ) / d ( tau(i) )
                  dgammadot_dtau(is) = abs(gammadot(is))*DeltaF/KB/T*q
     + *((1.- (ratio**p))**(q-1.))*p/tauc(is)*(ratio**(p-1.))
!
!
!
!                 Calculate derivative d ( gammadot(i) ) / d ( tauc(i) )
                  dgammadot_dtauc(is) = -abs(gammadot(is))*DeltaF/KB/T*q
     + *((1.- (ratio**p))**(q-1.))*p/tauc(is)*(ratio**p)*signtau
!
              else
                  gammadot(is)=0.
              end if

!
!
          endif
!
!

!
!
!
!         Contribution to Jacobian
          Pmat = Pmat +
     + dt*dgammadot_dtau(is) * SchmidxSchmid(is,:,:)
!
!         Plastic velocity gradient contribution
          Lp = Lp + gammadot(is)*Schmid_0(is,:,:)
!
!         Update plastic velocity gradient
          Dp = Dp + gammadot(is)*Schmid(is,:,:)
!
      end do
!
!     Find plastic flow direction
      Dp = (Dp + transpose(Dp))/2.          
!
!     Symmetrize Pmat
      Pmat = (Pmat + transpose(Pmat))/2.
!
!
      end subroutine doubleexpslip
!
!
!
!
!
!
!
!
!
!
!
! 
      subroutine powerslip(Schmid_0,
     + Schmid,SchmidxSchmid,
     + tau,X,tauc,
     + burgerv,dt,nslip,iphase,T,slipparam,
     + irradiationmodel,irradiationparam,
     + cubicslip,caratio,
     + Lp,Dp,Pmat,gammadot,dgammadot_dtau,
     + dgammadot_dtauc)
!
      use userinputs, only : useaveragestatevars, maxnparam
      use utilities,  only : gmatvec6
      implicit none
!     Schmid tensor at the intermediate config. (or undeformed)
      real(8), intent(in) :: Schmid_0(nslip,3,3)
!     Schmid tensor
      real(8), intent(in) :: Schmid(nslip,3,3)
!     Schmid dyadic
      real(8), intent(in) :: SchmidxSchmid(nslip,6,6)
!     RSS
      real(8), intent(in) :: tau(nslip)
!     Backstress
      real(8), intent(in) :: X(nslip)
!     CRSS
      real(8), intent(in) :: tauc(nslip)
!     Burger's vector
      real(8), intent(in) :: burgerv(nslip)
!     time increment
      real(8), intent(in) :: dt
!     number of slip systems
      integer, intent(in) :: nslip
!     phase id
      integer, intent(in) :: iphase
!     temperature
      real(8), intent(in) :: T
!     slip parameters
      real(8), intent(in) :: slipparam(maxnparam)
!     irradiation model id
      integer, intent(in) :: irradiationmodel
!     irradiation model parameters
      real(8), intent(in) :: irradiationparam(maxnparam)
!     cubic slip for fcc
      integer, intent(in) :: cubicslip
!     c/a ratio for hcp materials
      real(8), intent(in) :: caratio
!     plastic part of velocity gradient
      real(8), intent(out) :: Lp(3,3)
!     plastic stretch rate at the deformed configuration
      real(8), intent(out) :: Dp(3,3)
!     tangent matrix required for N-R iteration (at the inner loop)
      real(8), intent(out) :: Pmat(6,6)
!     slip rates
      real(8), intent(out) :: gammadot(nslip)
!     Derivative of slip rates wrto rss
      real(8), intent(out) :: dgammadot_dtau(nslip)
!     Derivative of slip rates wrto crss
      real(8), intent(out) :: dgammadot_dtauc(nslip)
!
!     Variables used within this subroutine
      integer :: is, i, j
      real(8) :: gammadot0, power_n, constant_n, slope_n, ratio
      real(8) :: abstau, signtau
!
!
!     Obtain slip parameters
!     Reference strain rate
      gammadot0 = slipparam(1)
!     rate sensitivity exponent- constant (n)
      constant_n = slipparam(2)
!     temperature dependence of exponent (dn/dT)
      slope_n = slipparam(3)
!
!
!
!
!     Temperature dependence of rate sensitivity exponent
      power_n = slope_n * T + constant_n
!
!
      Pmat = 0.; Lp = 0.; Dp = 0.
      gammadot = 0.
      dgammadot_dtau = 0.
      dgammadot_dtauc = 0.
!
!
      do is=1,nslip
!
!          
!         Absolute value of RSS
          abstau = abs(tau(is)-X(is))
!
!         Sign of RSS
          signtau = sign(1.0,tau(is)-X(is))
!      
!
          ratio = abstau / tauc(is)
!

!
          gammadot(is) = gammadot0*ratio**power_n*signtau
!
!
!
          dgammadot_dtau(is) = gammadot0*power_n
     + /tauc(is)*ratio**(power_n-1.0)
!
!
          dgammadot_dtauc(is) = -gammadot0*power_n
     + /tauc(is)*ratio**(power_n)*signtau
!
!
          Pmat = Pmat +
     + dt*dgammadot_dtau(is)*SchmidxSchmid(is,:,:)
!
!         Update plastic velocity gradient
          Lp = Lp + gammadot(is)*Schmid_0(is,:,:)
! 
!         Update plastic velocity gradient
          Dp = Dp + gammadot(is)*Schmid(is,:,:)
!
!
      end do
!
!     Find plastic flow direction
      Dp = (Dp + transpose(Dp))/2.          
!
!     Symmetrize Pmat
      Pmat = (Pmat + transpose(Pmat))/2.
!
!
!
      end subroutine powerslip
!
!
!
!
!
!
!
!
      end module slip