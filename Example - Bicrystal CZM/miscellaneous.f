!     Oct. 27th, 2025
!     Eralp Demir
!
      module miscellaneous
      implicit none
      contains
!
!     ********************************************************
!     **                 general functions                  **
!     ********************************************************
!
!
!
!
!     Subroutine to calculate slip and kink bands
!     Reference: https://doi.org/10.1016/j.actamat.2019.06.010
      subroutine SlipandKink(noel,npt,K,S)
      use userinputs, only: maxnslip
      use globalvariables, only:
     + statev_theta, statev_evmp,
     + numel, numpt 
      implicit none
      integer, intent(in) :: noel
      integer, intent(in) :: npt
      real(8), intent(out) :: K
      real(8), intent(out) :: S
      real(8) :: p_av, theta_av
      real(8) :: L, R
      
      
      p_av=sum(statev_evmp(:,:))/numpt/numel
      
      theta_av=sum(statev_theta(:,:))/numpt/numel
      
!     plastic slip
      if (statev_evmp(noel,npt).gt.(1.5*p_av)) then
          L = 1.
      else
          L =0.
      end if
      
!     rotation (changed from 3 to 2)
      if (statev_theta(noel,npt).gt.(3.*theta_av)) then
          R = 1.
      else
          R = 0.
      end if
      
!     Kink band
      K = L * R
      
!     Slip band
      S = L - K
      
      
      
      
      return
      end subroutine SlipandKink
!
!
!
!     Subroutine to calculate the active slip sytems
      subroutine SlipSystemActivity(noel,npt,SSA)
      use userinputs, only: maxnslip
      use globalvariables, only: statev_gammadot, smallnum,
     + materialid, numslip_all
      implicit none
      integer, intent(in) :: noel
      integer, intent(in) :: npt
      real(8), intent(out) :: SSA(maxnslip)
      integer :: is, nslip, matid
!
      SSA=0.
      matid = materialid(noel,npt)
      if (matid.ne.0) then
          nslip = numslip_all(matid)
          do is=1,nslip
              if (abs(statev_gammadot(noel,npt,is)).gt.smallnum) then
                  SSA(is)=1.
              end if
          end do
      end if
!
      return
      end subroutine SlipSystemActivity
!
!
!
!
!     Subroutine to calculate strain projections
      subroutine ProjectLatticeStrain(noel,npt,eps)
      use globalvariables, only: statev_Eec
      use utilities, only: vecmat6
      implicit none
      integer, intent(in) :: noel
      integer, intent(in) :: npt
      real(8), intent(out) :: eps
      real(8) :: hkl(3)
      real(8) :: Eec33(3,3), Eec(6)
      integer i, j
!
!     THIS IS DEFINED BY THE USER
!     Define hkl direction
      hkl(1) =  0.
      hkl(2) = -1.
      hkl(3) =  1.
!
!     Normalize
      hkl=hkl/norm2(hkl)
!
!     Read lattice strains
      Eec=statev_Eec(noel,npt,:)
!
!     Undo shears
      Eec(4:6)=Eec(4:6)/2.
!
!     Convert to a 3x3 matrix
      call vecmat6(Eec,Eec33)
!
!     Project
      eps=0.
      do i=1,3
          do j=1,3
              eps = eps +
     + hkl(i)*Eec33(i,j)*hkl(j)
          end do
      end do
!
!
      return
      end subroutine ProjectLatticeStrain
!
!
!
!
!
!
!     Subroutine to calculate Fatemi Socie parameter
      subroutine FatemiSocieParameter(noel,npt,FSp)
      use userinputs, only: maxnslip
      use globalvariables, only: statev_sigma, statev_gammasum,
     + statev_gmatinv, statev_tauceff, materialid, norc_0_all,
     + numslip_all 
      use utilities, only: vecmat6
      implicit none
      integer, intent(in) :: noel
      integer, intent(in) :: npt
      real(8), intent(out) :: FSp
!     Variables used in the subroutine
      real(8) :: sig6(6), sig3x3(3,3)
      real(8) :: gammasum(maxnslip)
      real(8) :: gmatinv(3,3)
      real(8) :: tauceff(maxnslip), tauceffmax
      integer :: matid
      real(8) :: norc(3), nors(3)
      real(8) :: shmax, sigmax, k
      integer :: ismax, is, i, j, nslip
!
!     Input parameter "k"
!     Reference: https://doi.org/10.1016/j.ijfatigue.2011.01.003
      k = 0.5
!
!
      FSp=0.
      sig6 = statev_sigma(noel,npt,:)
      gammasum = statev_gammasum(noel,npt,:)
      gmatinv = statev_gmatinv(noel,npt,:,:)
      tauceff = statev_tauceff(noel,npt,:)
      matid = materialid(noel,npt)
      if (matid.ne.0) then
          nslip = numslip_all(matid)
!
!
!         Find maximum slip
          shmax=0.; ismax=0
          do is=1,nslip
              if (abs(gammasum(is)).gt.shmax) then
                  shmax = abs(gammasum(is))
                  ismax = is
              end if
          end do
!
!         If there is no maximum or zero slip
          if (ismax.eq.0) then
!
              FSp=0.
!
!         If there some slip on any system
          else      
!             Maximum CRSS
              tauceffmax = tauceff(ismax)
!
!             Slip plane normal of maximum slip
              norc = norc_0_all(matid,ismax,:)
!
!             Transform to sample reference
              nors = matmul(gmatinv,norc)
!
!             Convert stress to 3x3 matrix
              call vecmat6(sig6,sig3x3)
!
!             Project stress to the slip plane of maximum slip
              sigmax = 0.
              do i = 1, 3
                  do j = 1, 3
                      sigmax = sigmax +
     + nors(i)*sig3x3(i,j)*nors(j)
                  end do
              end do
!
!             Check for zero tauceffmax
              if (tauceffmax.gt.0.) then
!                 Parameter calculation
                  FSp = shmax/2.*(1.+k*sigmax/tauceffmax)
              else
                  FSp = 0.
              end if
!             
          end if
      else
          FSp = 0.
      end if
!
      return
      end subroutine FatemiSocieParameter
!
!
!
!     Identify the neighboring points
      subroutine FindNeighbours
      use globalvariables, only: pi, 
     + numel, numpt, ipcoords, numdim, ipdomain,
     + numneigh, eleneigh, iptneigh, facneigh
      use userinputs, only: maxneigh, horizonR
      implicit none
!     Variables used in this subroutine
      integer :: iel, ipt, jel, jpt, k, ind
      real(8) :: d, sum2, Vi, Vj
!
!     Loop through each point
      do iel=1, numel
          
          do ipt = 1, numpt
              
              Vi = ipdomain(iel,ipt)
              
              
              ind=0
              do jel=1,numel
                  
                  if (jel.ne.iel) then
                  
                      do jpt=1,numpt
              
                         
                          Vj = ipdomain(jel,jpt)
                          sum2=0.
                          do k = 1, numdim
                      
                              sum2 = sum2 + 
     + (ipcoords(iel,ipt,k) - ipcoords(jel,jpt,k))**2 
                          
                          end do
                      
              
                          d = sqrt(sum2)
                      
                          if (d.le.horizonR) then
                      
                              ind=ind+1
                          
                              
!                             Store the first maxneigh - array size
                              if (ind.le.maxneigh) then
                                  numneigh(iel,ipt)=ind
                                  eleneigh(iel,ipt,ind)=jel
                                  iptneigh(iel,ipt,ind)=jpt
                                  facneigh(iel,ipt,ind)=
     + Vj/Vi*exp(-pi*sum2/horizonR**2)
                              end if
                          
                          end if
                      
              
                      end do
                  
                  end if

              end do
                  
          enddo
      enddo
      
      return
      end subroutine FindNeighbours
!
!
!
!
      end module miscellaneous