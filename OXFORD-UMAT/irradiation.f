!     Specific subroutines for irradiation models
!
!     Strength and hardening interaction matrices for irradiation model-2
      module irradiation
      implicit none
!
!
      contains
!
!
!     Subroutine to calculate strength interaction matrix for
!     irradiation model-2
!     Ref: https://doi.org/10.1016/j.actamat.2022.118361
      subroutine calculateintmats4irradmodel2(nslip, dirc, norc,
     + irradiationparam, sintmat, hintmat)
      use userinputs, only: maxnslip, maxnparam
      implicit none
!     Inputs
!     Number of slip systems
      integer, intent(in) :: nslip
!     Normalize slip directions
      real(8), intent(in) :: dirc(maxnslip,3)
!     Normalized slip plane normals
      real(8), intent(in) :: norc(maxnslip,3)
!     Irradiation model parameters
      real(8), intent(in) :: irradiationparam(maxnparam)
!     Outputs
!     Strength interaction matrix
!     Strength interaction:   sintmat (nslip x nloop)
      real(8), intent(out) :: sintmat(maxnslip,maxnslip)
!     Hardening (Softening) interaction matrix
!     Hardening interaction:  hintmat (nloop x nslip)
      real(8), intent(out) :: hintmat(maxnslip,maxnslip)
!     Variables used in this subroutine
      integer ::  nloop
!     reaction vector
      real(8) ::  r(3)
!     Projected vector (reaction segment)
      real(8) ::  a
      integer ::  i, j, ind
!
!
!
!
!     Reset arrays
      sintmat = 0.
      hintmat = 0.
!
!
!     Number of types of defects
      nloop = int(irradiationparam(1))
!
!
      do i=1,nslip
!
          do j=1,nloop
!
!             Slip system index corresponding to the defect type
              ind = int(irradiationparam(7+j))
!
!             What is the reaction product for the gliding dislocation on slip
!             system 'i' and defect type 'j'?
              r = dirc(i,:) + dirc(ind,:)
!
!             Is this reaction in the slip plane of slip system 'i'?
              a = dot_product(norc(i,:), r)
!
              if (a==0.) then
                  sintmat(i,j)=irradiationparam(12)
                  hintmat(j,i)=irradiationparam(14)
              else
                  sintmat(i,j)=irradiationparam(11)
                  hintmat(j,i)=irradiationparam(13)
              end if
!
          end do
!
      end do
!
!
!
!
!
!
!
      end subroutine calculateintmats4irradmodel2
!
!
!
      end module irradiation