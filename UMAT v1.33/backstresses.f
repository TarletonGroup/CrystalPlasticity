!     May 03rd, 2022
!     Eralp Demir
!
      module backstresses
      implicit none
      contains
!
!     ********************************************************
!     ** crss calculates the crss of slip systems           **
!     ********************************************************
      subroutine backstressmodel1
      use userinputs, only: numel, maxnslip
      use globalvariables, only: numpt,
     + phaseid, numslip_all, numscrew_all, screw_all,
     + burgerv_all, gf_all, G12_all, v12_all,
     + dirc_0_all, norc_0_all, statev_gmatinv,
     + statev_backstress_t, statev_gnd
      use utilities, only: matvec6
      implicit none
      integer :: matid, nslip, nscrew
      real(8) :: burgerv(maxnslip)
      integer :: screw(maxnslip)
      real(8) :: dirc_0(maxnslip,3)
      real(8) :: norc_0(maxnslip,3)
      real(8) :: gmatinv(3,3)
      real(8) :: dirs(3)
      real(8) :: nors(3)
      real(8) :: X(maxnslip)
      real(8) :: dyad33(3,3), dyad6(6)
      real(8) :: sigmab(6)
      real(8) :: gf, G12, v12
      real(8) :: rhoGNDe, rhoGNDs
      integer :: ie, ip, is, i, k, l
!
!
!
!
!
!
      do ie=1, numel
!
!         Reset arrays
          burgerv=0.;screw=0
          dirc_0=0.; norc_0=0.
!
!
!
!         Assume the same material for al the Gaussian points of an element
          matid = phaseid(ie,1)
!
!         Number of slip systems
          nslip = numslip_all(matid)
!
!         Number of screw systems
          nscrew = numscrew_all(matid)
!
!         Screw systems
          screw = screw_all(matid,1:nscrew)
!
!         Geometric factor
          gf = gf_all(matid)
!
!         Shear Modulus
          G12 = G12_all(matid)
!
!         Poisson's ratio
          v12 = v12_all(matid)
!
!         Burgers vector
          burgerv(1:nslip) = burgerv_all(matid,1:nslip)
!
!
!         undeformed slip direction
          dirc_0(1:nslip,1:3) = dirc_0_all(matid,1:nslip,1:3)
!
!         undeformed line direction
          norc_0(1:nslip,1:3) = norc_0_all(matid,1:nslip,1:3)
!
!
          do ip=1,numpt
!
!
!
!
!
!
!             Reset backstress
              X = 0.
!
!
!             Edges
              do is = 1, nslip
!
!
!                 Edge dislocation density
                  rhoGNDe = statev_gnd(ie,ip,is)
!
!
!
!
!                 Backstress due to edges
                  X(is) = -gf*G12*burgerv(is)*
     + sqrt(abs(rhoGNDe))*sign(1.0,rhoGNDe)/(1.-v12)
!
!
!
              end do
!
!
!
!             Screws
              do i = 1, nscrew
!
                  is = screw(i)
!
!
!                 Screw dislocation density
                  rhoGNDs = statev_gnd(ie,ip,nslip+i)
!
!                 Backstress due to screws
                  X(is) = X(is) + gf*G12*burgerv(is)*
     + sqrt(abs(rhoGNDs))*sign(1.0,rhoGNDs)
!
!
!
              end do
!
!
!
!
!             Crystal to sample transformation matrix
              gmatinv = statev_gmatinv(ie,ip,:,:)              
!
!             Sum the backstress components in dyadic form
!             Reset backstress    
              sigmab = 0.
!
!             To obtain backstress tensor
              do is = 1, nslip
!                 Transform slip directions to sample reference
                  dirs = matmul(gmatinv, dirc_0(is,:))
!
!                 Transform line directions to sample reference
                  nors = matmul(gmatinv, norc_0(is,:))
!
!                 Dyadic product
                  do k=1,3
                      do l=1,3
                          dyad33(k,l) = dirs(k) * nors(l)
                      end do
                  end do
!
!                 Convert 3x3 matrix to 6x1 vector
                  call matvec6(dyad33,dyad6)
!
!
!
!                 Sum the backstress
                  sigmab = sigmab + dyad6 * X(is)
!
!
              end do
!
!
!
!
!             Assign to the global vector
              statev_backstress_t(ie,ip,1:6) = sigmab
!
!
!
!
!
!
          end do
!
!
!
!
!
!
!
      end do
!
!
!
!
!
!
!
      return
!
      end subroutine backstressmodel1
!
!
      end module backstresses