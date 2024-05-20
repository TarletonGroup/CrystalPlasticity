      module innerloop
      implicit none
      contains
!
!
      subroutine Dunne_innerloop(
     + Schmid_0,Schmid,
     + SchmidxSchmid,Schmidvec,
     + phaid,nslip,mattemp,Cs,
     + burgerv,cubicslip,caratio,
     + slipparam,slipmodel,
     + creepmodel,creepparam,
     + irradiationmodel,
     + irradiationparam,
     + dt,sigmatr,tauceff,
     + rhofor,X,gammasum,
     + sigma,tau,cpconv,
     + gammadot,Lp,
     + Dp,dstranp33,
     + invdpsi_dsigma,iter)
!
      use globalvariables, only : I6
!
      use userinputs, only : maxniter, tolerance, 
     + maxnparam, maxnloop, SVDinversion, maxnslip
!
      use utilities, only : matvec6,
     + nolapinverse, SVDinverse 
!
      use slip, only : sinhslip, doubleexpslip,
     + powerslip
!
      use creep, only : expcreep
!
      use errors, only : error
!
      implicit none
!
!     INPUTS    
!
!     Initial Schmid tensor
      real(8), intent(in) :: Schmid_0(nslip,3,3) 
!     Schmid tensor
      real(8), intent(in) :: Schmid(nslip,3,3)
!     Schmid dyadic
      real(8), intent(in) :: SchmidxSchmid(nslip,6,6)
!     Vectorized Schmid tensor
      real(8), intent(in) :: Schmidvec(nslip,6)
!     phase-id
      integer, intent(in) :: phaid
!     number of slip sytems
      integer, intent(in) :: nslip
!     temperature
      real(8), intent(in) :: mattemp
!     elastic compliance
      real(8), intent(in) :: Cs(6,6)
!     Burgers vectors
      real(8), intent(in) :: burgerv(nslip)
!     flag for cubic slip systems
      integer, intent(in) :: cubicslip 
!     c/a ratio for hcp crystals
      real(8), intent(in) :: caratio
!     slip model no.
      integer, intent(in) :: slipmodel
!     slip model parameters
      real(8), intent(in) :: slipparam(maxnparam)
!     creep model no.
      integer, intent(in) :: creepmodel
!     creep model parameters
      real(8), intent(in) :: creepparam(maxnparam)      
!     irrradiation model no.
      integer, intent(in) :: irradiationmodel
!     irradiation model parameters
      real(8), intent(in) :: irradiationparam(maxnparam)        
!
!     time increment
      real(8), intent(in) :: dt  
!     trial stress
      real(8), intent(in) :: sigmatr(6)
!     overall crss
      real(8), intent(in) :: tauceff(nslip)
!     total forest dislocation density per slip system at the former time step
      real(8), intent(in) :: rhofor(nslip)
!     crss at the current time step
      real(8), intent(in) :: X(nslip)
!     total slip per slip system accumulated over the time
!     at the current time step
      real(8), intent(in) :: gammasum(nslip)
!
!     INOUTS
!     Cauchy stress
      real(8), intent(inout) :: sigma(6)
!     overall crss
      real(8), intent(inout) :: tau(nslip)
!     convergence flag
      integer, intent(inout) :: cpconv
!
!     OUTPUTS
!     slip rates at the current time step
      real(8), intent(out) :: gammadot(nslip)  
!     plastic velocity gradient
      real(8), intent(out) :: Lp(3,3)
!     Total deformation rate
      real(8), intent(out) :: Dp(3,3)
!     Total deformation rate
      real(8), intent(out) :: dstranp33(3,3)
!     
!     Jacobian of the Newton-Raphson loop
!     and its inverse
      real(8), intent(out) :: invdpsi_dsigma(6,6)
!     number of iterations
      integer, intent(out) :: iter

!
!     Local variables used within this subroutine       
!        
!     plastic velocity gradient for slip
      real(8) Lp_s(3,3)
!     plastic velocity gradient for creep
      real(8) Lp_c(3,3)
!     Sym part of plastic velocity gradient for slip
      real(8) Dp_s(3,3)
!     Sym part of plastic velocity gradient for creep
      real(8) Dp_c(3,3)
!     plastic tangent stiffness for slip
      real(8) Pmat_s(6,6)
!     plastic tangent stiffness for creep
      real(8) Pmat_c(6,6)
!     tangent matrix for NR iteration
      real(8) Pmat(6,6)
!     slip rates for slip
      real(8) gammadot_s(nslip)
!     slip rates for creep
      real(8) gammadot_c(nslip)
!     derivative of slip rates wrto rss for slip
      real(8) dgammadot_dtau_s(nslip)
!     derivative of slip rates wrto rss for creep
      real(8) dgammadot_dtau_c(nslip)
!     derivative of slip rates wrto rss for slip
      real(8) dgammadot_dtauc_s(nslip)
!     derivative of slip rates wrto rss for creep
      real(8) dgammadot_dtauc_c(nslip)
!
!     plastic strain increment
      real(8) :: dstranp(6)
!
!
!     Jacobian of the Newton-Raphson loop
!     and its inverse
      real(8)  :: dpsi_dsigma(6,6)
!     residual of the Newton-Raphson loop
!     vector and scalar
      real(8) :: psinorm, psi(6)
!
!
!
!     stress increment
      real(8) :: dsigma(6)
!
!     error flag for svd inversion
      integer :: err
!
      integer :: is
!
!     Reset variables for the inner iteration
      psinorm = 1.
      iter = 0
!
!
!     Newton-Raphson (NR) iteration to find stress increment
      do while ((psinorm >= tolerance)
     +.and.(iter < maxniter).and.(cpconv == 1))
!
!
!         Slip models to find slip rates
!
!         none
          if (slipmodel == 0) then
!
              Lp_s = 0.
              Dp_s = 0.
              Pmat_s = 0.
              gammadot_s = 0.
!
!
!         sinh law
          elseif (slipmodel == 1) then
!
              call sinhslip(Schmid_0,Schmid,SchmidxSchmid,
     + tau,X,tauceff,rhofor,burgerv,dt,
     + nslip,phaid,mattemp,slipparam,
     + irradiationmodel,irradiationparam,
     + cubicslip,caratio,Lp_s,Dp_s,Pmat_s,
     + gammadot_s,dgammadot_dtau_s,
     + dgammadot_dtauc_s)
!
!
!         exponential law
          elseif (slipmodel == 2) then
!
!
              call doubleexpslip(Schmid_0,Schmid,SchmidxSchmid,
     + tau,X,tauceff,burgerv,dt,nslip,phaid,
     + mattemp,slipparam,irradiationmodel,
     + irradiationparam,cubicslip,caratio,
     + Lp_s,Dp_s,Pmat_s,gammadot_s,dgammadot_dtau_s,
     + dgammadot_dtauc_s)
!
!
!         power law
          elseif (slipmodel == 3) then
!
!
              call powerslip(Schmid_0,Schmid,SchmidxSchmid,
     + tau,X,tauceff,burgerv,dt,
     + nslip,phaid,mattemp,slipparam,
     + irradiationmodel,irradiationparam,
     + cubicslip,caratio,Lp_s,Dp_s,Pmat_s,
     + gammadot_s,dgammadot_dtau_s,
     + dgammadot_dtauc_s)
!
!
          end if
!
!
!
!         Slip due to creep
          if (creepmodel == 0) then
!
!
              Lp_c = 0.
              Dp_c = 0.
              Pmat_c = 0.
              gammadot_c = 0.
!
!
          elseif (creepmodel == 1) then
!
!
!
              call expcreep(Schmid_0,Schmid,SchmidxSchmid,
     + tau,X,tauceff,dt,nslip,phaid,
     + mattemp,creepparam,gammasum,Lp_c,Dp_c,Pmat_c,
     + gammadot_c,dgammadot_dtau_c,
     + dgammadot_dtauc_c)
!
!
!
!
          endif
!
!
!         Sum the effects of creep and slip rates
          Lp = Lp_s + Lp_c
          Dp = Dp_s + Dp_c
          Pmat = Pmat_s + Pmat_c
          gammadot = gammadot_s + gammadot_c
!
!
!
!         Check for NaN in the slip rates
          if(any(gammadot/=gammadot)) then
!             did not converge
              cpconv = 0
              return
          endif
!
!
!         Check for Inf in the slip rate vector
          if(any(gammadot*0./=gammadot*0.)) then
!             did not converge
              cpconv = 0
              return
          endif
!
!         Check for the Pmat
          if(any(Pmat /= Pmat)) then
!             did not converge
              cpconv = 0
              return
          endif
!
!
!
!         Plastic strain increment
          dstranp33 = Dp*dt
          call matvec6(dstranp33,dstranp)
          dstranp(4:6)=2.*dstranp(4:6)
!
!
!
!
!         Tangent-stiffness calculation
!         Jacobian of the Newton loop (see Dunne, Rugg, Walker, 2007)
          dpsi_dsigma = I6 + matmul(Cs, Pmat)
!

!

!
!         Invert (double precision version)
          call nolapinverse(dpsi_dsigma,invdpsi_dsigma,6)

!
!
!         If inversion is not successfull!
!         Check for the inverse
          if(any(invdpsi_dsigma /= invdpsi_dsigma)) then
!
!             Try using singular value decomposition
!             If singular value decomposition is ON
              if (SVDinversion==1) then
!
!                 Invert
                  call SVDinverse(dpsi_dsigma,6,invdpsi_dsigma,err)
!
!
!
              else
!
!
!                 did not converge
                  err = 1
!
!
!
              end if
!
!             Check again and if still not successfull
              if(err==1) then
!                 did not converge
                  cpconv = 0
                  return
              end if
!
!
          endif
!
!
!
!         residual (predictor - corrector scheme)
          psi = sigmatr - sigma - matmul(Cs,dstranp)
!
!         norm of the residual
          psinorm = sqrt(sum(psi*psi))
!
!
!         stress increment
          dsigma = matmul(invdpsi_dsigma,psi)
!
!
!         stress update
          sigma = sigma + dsigma
!
!
!
!         calculate resolved shear stress on slip systems
!         rss and its sign
          do is = 1, nslip
              tau(is) = dot_product(Schmidvec(is,:),sigma)
          end do
!
!
!         increment iteration no.
          iter = iter + 1
!
!     End of NR iteration (inner loop)
      end do
!
!
!
!
!         
      return
      end subroutine Dunne_innerloop
!
!
!
!
!
!
!
!
!     This subroutine is written by Chris Hardie (11/10/2023)      
!     Explicit state update rule
!     Solution using state variables at the former time step
      subroutine Hardie_innerloop(
     + Schmid_0,Schmid,
     + SchmidxSchmid,Schmidvec,
     + phaid,nslip,mattemp,Cs,
     + burgerv,cubicslip,caratio,
     + slipparam,slipmodel,
     + irradiationmodel,
     + irradiationparam,
     + dt,sigmatr,
     + abstautr,signtautr,
     + tauceff,rhofor,X,
     + sigma,iterinverse)
!
      use globalvariables, only : I3, I6, smallnum
!
      use userinputs, only : maxniter, maxnparam, maxnslip,
     + inversetolerance 
!
      use utilities, only : matvec6, gmatvec6
!
      use slipreverse, only : sinhslipreverse, powerslipreverse, 
     + doubleexpslipreverse
!
!
      implicit none
!
!
!     INPUTS
!
!     Initial Schmid tensor
      real(8), intent(in) :: Schmid_0(nslip,3,3) 
!     Schmid tensor
      real(8), intent(in) :: Schmid(nslip,3,3)
!     Schmid dyadic
      real(8), intent(in) :: SchmidxSchmid(nslip,6,6)
!     Vectorized Schmid tensor
      real(8), intent(in) :: Schmidvec(nslip,6)
!     phase-id
      integer, intent(in) :: phaid
!     number of slip sytems
      integer, intent(in) :: nslip
!     temperature
      real(8), intent(in) :: mattemp
!     elastic compliance
      real(8), intent(in) :: Cs(6,6)
!     Burgers vectors
      real(8), intent(in) :: burgerv(nslip)
!     flag for cubic slip systems
      integer, intent(in) :: cubicslip 
!     c/a ratio for hcp crystals
      real(8), intent(in) :: caratio
!     slip model no.
      integer, intent(in) :: slipmodel
!     slip model parameters
      real(8), intent(in) :: slipparam(maxnparam)    
!     irrradiation model no.
      integer, intent(in) :: irradiationmodel
!     irradiation model parameters
      real(8), intent(in) :: irradiationparam(maxnparam)        
!
!     time increment
      real(8), intent(in) :: dt  
!     trial stress
      real(8), intent(in) :: sigmatr(6)
!     absolute value of trial RSS
      real(8), intent(in) :: abstautr(nslip)
!     sign of trial RSS
      real(8), intent(in) :: signtautr(nslip)
!     overall crss
      real(8), intent(in) :: tauceff(nslip)
!     total forest dislocation density per slip system at the former time step
      real(8), intent(in) :: rhofor(nslip)
!     crss at the current time step
      real(8), intent(in) :: X(nslip)
!
!     OUTPUTS
!     Cauchy stress
      real(8), intent(out) :: sigma(6)
!     Number of iterations
      integer, intent(out) :: iterinverse
!
!
!     Local variables used within this subroutine       
!           
!     plastic velocity gradient for slip
      real(8) Lp_s(3,3)
!     slip rates for slip
      real(8) gammadot_s(nslip)
!     absolute slip rates      
      real(8) absgammadot_s(nslip)
!     sign of gammadot_s      
      real(8) signgammadot_s(nslip)
!     derivative of rss wrt slip rates for slip
      real(8) dtau_dgammadot_s(nslip,nslip)
!     derrivative of error function      
      real(8) dpsi_dgammadot(nslip,nslip)
!     inverse of derivative above      
      real(8) dpsi_dgammadotinv(nslip,nslip)
!     slip correction
      real(8) dgammadot(nslip)
!     Derivative of slip law in forward direction
      real(8) :: dgammadot_dtau(nslip)
!     rss at the former time step
      real(8) :: tau(nslip), tau_e(nslip)
!     absolute value of rss at the former time step
      real(8) :: abstau(nslip)
!     sign of rss at the former time step
      real(8) :: signtau(nslip)
!
!     residual of the Newton-Raphson loop
!     vector and scalar
      real(8) :: psinorm, psinorm2, psi(nslip)
!
!     Forward Jacobian
      real(8) :: Pmat(6,6)
      real(8) :: dpsi_dsigma(6,6)
!     LU decomposition matricies for calculation of determinant
!      real(8), allocatable :: l(:,:), u(:,:)
!      integer, allocatable :: p(:,:)
!
!     plastic strain increment
      real(8) :: plasstraininc33(3,3), plasstraininc(6)
      real(8) :: plasstrainrate(3,3)
!
!     Shear Stiffness Matrix
      real(8) :: G(nslip,nslip)
!
!     other variables
      real(8) :: dummy6(6), damping, iter_max, activesum
      integer :: is
!
!      allocate(dpsi_dsigma(6,6),l(6,6),u(6,6),p(6,6))
!
!
!     Reset variables for iteration
      iter_max=10000.0
      iterinverse = 0
      psinorm = 1.0d+200
      psinorm2 = 1.0d+200
!
!     Initial guess for NR scheme
!     Stress at the former time step
      sigma = sigmatr
      abstau = abstautr
      signtau =signtautr
      tau_e = abstau*signtau
      gammadot_s=0.
      absgammadot_s=0.
      Lp_s=0.
!
!     Build Shear stiffness matrix
!
      G = 0.
      do is = 1, nslip
          dummy6 = matmul(Cs, Schmidvec(is,:))
          G(is,is) = dot_product(Schmidvec(is,:), dummy6)
      end do
!
!     Newton-Raphson (NR) iteration to find stress increment
      do while ((psinorm >= inversetolerance))!  .OR. (psinorm2 >= tol2))!((psinorm >= tol)) !((psinorm >= tol)  .AND. (psinorm2 >= tol2))
!
!         increment iteration no.
          iterinverse = iterinverse + 1
!
!         Slip models to find slip rates
!
!         none
          if (slipmodel == 0) then
!
              gammadot_s = 0.
!
!         sinh law
          elseif (slipmodel == 1) then
!
              call sinhslipreverse(
     + Schmid,SchmidxSchmid,signtau,
     + abstau,tau,X,tauceff,rhofor,burgerv,dt,
     + nslip,phaid,mattemp,slipparam,
     + irradiationmodel,irradiationparam,
     + cubicslip,caratio,Lp_s,
     + absgammadot_s, gammadot_s,
     + dtau_dgammadot_s, dgammadot_dtau,Pmat)
!
!
!
!         double exponent law (exponential law)
          elseif (slipmodel == 2) then
!
!
              call doubleexpslipreverse(
     + Schmid,SchmidxSchmid,
     + tau,X,abstau,signtau,tauceff,
     + burgerv,dt,nslip,phaid,mattemp,slipparam,
     + irradiationmodel,irradiationparam,
     + cubicslip,caratio,
     + Lp_s,Pmat,absgammadot_s,gammadot_s,
     + dtau_dgammadot_s,dgammadot_dtau) 
!
!
!         power law
          elseif (slipmodel == 3) then
!
!
              call powerslipreverse(
     + Schmid,SchmidxSchmid,
     + tau,X,abstau,signtau,tauceff,
     + burgerv,dt,nslip,phaid,mattemp,slipparam,
     + irradiationmodel,irradiationparam,
     + cubicslip,caratio,
     + Lp_s,absgammadot_s, gammadot_s,dtau_dgammadot_s,
     + dgammadot_dtau, Pmat)  
!
!
          end if
!
!         calculate resolved shear stress on slip systems
!         rss and its sign
          activesum=0
          do is = 1, nslip
              tau(is) = signtau(is)*abstau(is)
              if (abstau(is) .GE. tauceff(is)) then
                  activesum=activesum+1.0
              end if
          end do             
!
!         residual (predictor - corrector) including backstress
          psi = tau - X - tau_e
!
!         norm of the residual
          psinorm2 = sqrt(sum(psi*psi))
!
!         Forward Jacobian
          dgammadot_dtau=abs(dgammadot_dtau)*dt
          psinorm = maxval(dgammadot_dtau)
!
!         ! CALL LU_DECOMP(dpsi_dsigma, p, l, u)
!
!
!
!         dpsi_dgammadot
!
          dpsi_dgammadot = dtau_dgammadot_s + G*dt
!
!         Invert diagonal matrix
          dpsi_dgammadotinv=0.
          do is = 1, nslip
              dpsi_dgammadotinv(is,is)=1/dpsi_dgammadot(is,is)
          end do
!
!
          damping=min(2.0/activesum, 1.0)!0.5)
!          damping = 0.5
!          if (psinorm > 200.0) then
!              damping=min(2.0/activesum, 0.5)
!          end if
!
!         slip increment
          dgammadot = matmul(dpsi_dgammadotinv,psi)
!
          gammadot_s = gammadot_s-damping*dgammadot
!
!
!         Calculate Plastic Strain
          Lp_s=0.; plasstrainrate=0.
          do is = 1, nslip
              Lp_s = Lp_s + gammadot_s(is)*Schmid_0(is,:,:)
              plasstrainrate = plasstrainrate +
     + gammadot_s(is)*Schmid(is,:,:)
              absgammadot_s(is) = abs(gammadot_s(is))
              signgammadot_s(is)= sign(1.0,gammadot_s(is))
          end do
!
!         plastic strain rate
          plasstrainrate = (plasstrainrate + 
     + transpose(plasstrainrate))/2.
!
!
!         Plastic strain increment
          plasstraininc33 = plasstrainrate*dt
          call matvec6(plasstraininc33,plasstraininc)
          plasstraininc(4:6)=2.*plasstraininc(4:6)
          call gmatvec6(plasstraininc33,plasstraininc)
!
!         Calculate rss following plastic relaxation
!
          sigma=sigmatr-matmul(Cs,plasstraininc)
!
!         calculate resolved shear stress on slip systems
!         rss and its sign
          do is = 1, nslip
              tau_e(is) = dot_product(Schmidvec(is,:),sigma)
              tau(is) = tau_e(is)
              signtau(is) = sign(1.0,tau(is))
              abstau(is) = abs(tau(is))
          end do          
!
!         check if slip is changing direction
          do is = 1, nslip
              if (signgammadot_s(is) /= signtau(is)) then
                  gammadot_s(is)=0.0
                  absgammadot_s(is)=0.0
              end if
          end do
!
          if (iterinverse > iter_max) then
              psinorm=1e-10
              psinorm2=1e-10
          end if
!
!     End of NR iteration for stress approximation
      end do
          !active_s=1
          !do is = 1, nslip
          !    if (absgammadot_s(is) > 0.0) then
          !        active_s(is)=1
          !    end if
          !end do  

      return
      end subroutine Hardie_innerloop
!
      end module innerloop