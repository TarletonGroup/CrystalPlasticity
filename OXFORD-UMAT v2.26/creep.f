!     Oct. 6th, 2022
!     Eralp Demir
!     Creep laws
!     1. exponential law (used for Ni-superalloy)
!
      module creep
      implicit none
      contains
!
!
      subroutine expcreep(Schmid_0,
     + Schmid,SchmidxSchmid,tau,X,tauc,
     + dtime,nslip,iphase,T,creepparam,slipsysplasstran,
     + Lp,Dp,Pmat,gammadot,dgammadot_dtau,dgammadot_dtauc)
!
      use globalvariables, only : Rgas
      use utilities, only : gmatvec6
      use userinputs, only: maxnparam
      implicit none
!     Schmid tensor - undeformed
      real(8), intent(in) :: Schmid_0(nslip,3,3)
!     Schmid tensor - deformed
      real(8), intent(in) :: Schmid(nslip,3,3)
!     Schmid dyadic
      real(8), intent(in) :: SchmidxSchmid(nslip,6,6)
!     RSS
      real(8), intent(in) :: tau(nslip)
!     Backstress
      real(8), intent(in) :: X(nslip)
!     CRSS
      real(8), intent(in) :: tauc(nslip)
!     time increment
      real(8), intent(in) :: dtime
!     number of slip systems
      integer, intent(in) :: nslip
!     phase id
      integer, intent(in) :: iphase
!     temperature
      real(8), intent(in) :: T
!     slip parameters
      real(8), intent(in) :: creepparam(maxnparam)
!     accumulated plastic strain on each slip system, signed
      real(8), intent(in) :: slipsysplasstran(nslip)
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
!     variables used within this subroutine
      real(8) :: gammadotc, gammadotd, bc, bd, Qc, Qd
      real(8) :: abstau, signtau
      integer :: is
!
!
!
!
!     Obtain creep parameters
!     reference creep rate (1/s)
      gammadotc = creepparam(1)
!     creep stress multiplier (1/MPa)
      bc = creepparam(5)
!     activation energy for creep (J/mol)
      Qc = creepparam(3)
!     reference damage rate (1/s)
      gammadotd = creepparam(4)
!     damage stress multiplier (1/MPa)
      bd = creepparam(5)
!     activation energy for damage (J/mol)
      Qd = creepparam(6)
!
!
!
      Pmat = 0.; Lp = 0.; Dp = 0.
!
      gammadot=0.
      dgammadot_dtau=0.
      dgammadot_dtauc=0.
!
!
!
!     contribution to Lp of all slip systems
      do is=1,nslip
!
!          
!         Absolute value of RSS
          abstau = abs(tau(is)-X(is))
!
!         Sign of RSS
          signtau = sign(1.0,tau(is)-X(is))
!
!          This statement is redundant so commented out!
!          if (abstau > 0.0) then
!
!
              gammadot(is) = 
     + signtau*gammadotc*exp(bc*abstau-Qc/Rgas/T) +
     + signtau*abs(slipsysplasstran(is))*gammadotd*
     + exp(bd*abstau-Qd/Rgas/T)
!
              dgammadot_dtau(is) = gammadotc*bc*
     + exp(bc*abstau-Qc/Rgas/T) +
     + abs(slipsysplasstran(is))*
     + gammadotd*bd*exp(bd*abstau-Qd/Rgas/T)
!
!
!
!
!             contribution to Jacobian
              Pmat = Pmat + dtime*dgammadot_dtau(is)
     + *SchmidxSchmid(is,:,:)
!
!             plastic velocity gradient contribution
              Lp = Lp + gammadot(is)*Schmid_0(is,:,:)
!
!
!             plastic velocity gradient contribution
              Dp = Dp + gammadot(is)*Schmid(is,:,:)
!
!
!          endif
!
!
!
      end do
!
!
!
!
!     Find plastic flow direction
      Dp = (Dp + transpose(Dp))/2.          
!
!     Symmetrize Pmat
      Pmat = (Pmat + transpose(Pmat))/2.
!
!      
!
!     
      end subroutine expcreep
!
!
!
!
!
      end module creep