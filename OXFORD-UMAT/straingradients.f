!     Jan. 1st, 2023
!     Eralp Demir
!
      module straingradients
      implicit none
!
      contains
!
!
!
!
!
!     GND model-1
!     GND calculation using curl of (Fp) together with L2 approximation
!     Cumulative (total) calculation of GNDs
!     Shutting down non-active slip systems followed by singular value decompostion
!     Proposed by Chris Hardie
      subroutine gndmodel1
      use userinputs, only : maxnslip,
     + gndhomogenization, gndthreshold
      use globalvariables, only : numel, numdim, numpt, nnpel,
     + materialid, numslip_all, numscrew_all, screw_all,
     + slip2screw_all, burgerv_all, dirc_0_all, trac_0_all, gradip2ip,
     + eijk, I3, statev_curvature, statev_Lambda, statev_gammasum,
     + statev_Fp, statev_gmatinv_0, statev_gnd, statev_gnd_t
      use utilities, only: matvec9
      implicit none
!     Local variables
      integer matid, nslip, nscrew
!     Some variables are allotable since material type can vary hence
!     the NUMBER OF SLIP SYSTEMS can alter within the mesh
      real(8) :: burgerv(maxnslip)
      real(8) :: gammasum(maxnslip)
      integer :: screw(maxnslip)
      real(8) :: slip2screw(maxnslip,maxnslip)
      real(8) :: dirc_0(maxnslip,3)
      real(8) :: trac_0(maxnslip,3)
!
      real(8) :: dirs(maxnslip,3)
      real(8) :: tras(maxnslip*2,3)
      real(8) :: Bmat(maxnslip*2,9)
      real(8) :: rhoGND(maxnslip*2)
      real(8) :: drhoGND(maxnslip*2)
!
      real(8) :: grad_invN(3,numpt), grad(3,3,3)
      real(8) :: Lambda(3,3), sum, Lambda_vec(9)
      real(8) :: kappa(3,3), kappa_vec(9), trace
!
      real(8) :: Fp_ip(numpt,3,3)
      real(8) :: gmatinv(3,3)
!
      integer :: i, j, k, l
      integer :: ie, ip, is
!
!
!
!
!     Loop through the elements
      do ie=1,numel
!
!         Reset arrays
          burgerv=0.; screw=0
          dirc_0=0.; trac_0=0.
          dirs=0.; tras=0.
          slip2screw=0.
!
!
!
!
!         Assume the same material for al the Gaussian points of an element
          matid = materialid(ie,1)
!
!         Number of slip systems
          nslip = numslip_all(matid)
!
!         Number of screw systems
          nscrew = numscrew_all(matid)
!
!         Screw systems
          screw(1:nscrew) = screw_all(matid,1:nscrew)
!
!         Slip to screw mapping
          slip2screw(1:nscrew,1:nslip) = 
     + slip2screw_all(matid,1:nscrew,1:nslip)
!
!         Burgers vector
          burgerv(1:nslip) =  burgerv_all(matid,1:nslip)
!
!         undeformed slip direction
          dirc_0(1:nslip,1:3) = dirc_0_all(matid,1:nslip,1:3)
!
!         undeformed line direction
          trac_0(1:nslip,1:3) = trac_0_all(matid,1:nslip,1:3)
!
!
!         Store Fp for each IP
!         Plastic part of the deformation gradient
          Fp_ip = statev_Fp(ie,1:numpt,1:3,1:3)
!
!
!         Calculate the gradient of Fp
!         Calculate the gradients using gradient operator
          do ip = 1, numpt
!
!
!             Reset arrays
              Bmat=0.; rhoGND=0.; gammasum=0.
!
!             Crystal to sample transformation matrix
              gmatinv = statev_gmatinv_0(ie,ip,:,:)
!
!
!             Slip rate per slip system
              gammasum(1:nslip) = statev_gammasum(ie,ip,1:nslip)
!
!
              do is = 1, nslip
!                 Transform slip directions to sample reference
                  dirs(is,:) = matmul(gmatinv, dirc_0(is,:))
!
!                 Transform line directions to sample reference
                  tras(is,:) = matmul(gmatinv, trac_0(is,:))
!
              end do
!
!
!             Calculate Bmatrix - using singular value decomposition
!             Arsenlis, A. and Parks, D.M., 1999. Acta materialia, 47(5), pp.1597-1611.
              call calculateBmatPINV(nslip,nscrew,
     + screw(1:nscrew), slip2screw(1:nscrew,1:nslip),
     + dirs(1:nslip,:),tras(1:nslip,:),burgerv(1:nslip),
     + gammasum(1:nslip), Bmat(1:nslip+nscrew,1:9))
!
!
!             Use gradients per integration point
              if (gndhomogenization == 0) then
!
                  grad_invN = gradip2ip(ie,ip,:,:)
!
!             Use the gradient at the element center
              else
!
                  grad_invN = gradip2ip(ie,numpt+1,:,:)
!
              end if
!
!
!
!
!             grad(k,i,j) = Fp(i,j,k)
              do k = 1,3
                  do j = 1,3
                      do i = 1,3
                          grad(k,i,j) =
     + dot_product(grad_invN(k,:),Fp_ip(:,i,j))
                      end do
                  end do
              end do
!
!
!             calculate curl
!             NOTE THE NEGATIVE SIGN AND TRANSPOSE ARE MISSING IN THE ORIGINAL REFERENCE
!             lambda(l,k) = -eijk(i,j,k) * Fp(l,j,i)
!             index "i" refers to the gradient direction
              do k = 1,3
                  do l = 1,3
                      sum = 0.
                      do i = 1, 3
                          do j = 1, 3
                              sum = sum - eijk(i,j,k)*grad(i,l,j)
!!                             earlier version (no "-" sign)
!                              sum = sum + eijk(i,j,k)*grad(i,l,j)
                          end do
                      end do
                      Lambda(l,k) = sum
                  end do
              end do
!
!             Vectorize the incompatibility
              call matvec9(Lambda,Lambda_vec)
!
!
!             Assign the Incompatibility
              statev_Lambda(ie,ip,:) = Lambda_vec
!
!
!             Trace
              trace = Lambda(1,1) + Lambda(2,2) + Lambda(3,3)
!
!             Calculate curvature (negative transpose + trace term)
              kappa = -transpose(Lambda) + I3 * trace / 2.
!
!             Convert curvature to vector
              call matvec9(kappa,kappa_vec)
!
!
!             Assign curvature
              statev_curvature(ie,ip,1:9) = kappa_vec(1:9)
!
!
!
!
!             Reset arrays
              rhoGND=0.; drhoGND=0.
!
!             Compute the dislocation densities using L2 minimization
              rhoGND(1:nslip+nscrew) =
     + matmul(Bmat(1:nslip+nscrew,1:9),Lambda_vec)
!
!
!             Calculate the increment of GNDs
              drhoGND(1:nslip+nscrew) =
     + rhoGND(1:nslip+nscrew) - statev_gnd_t(ie,ip,1:nslip+nscrew)
!
!
!             Check for a threshold
              do is = 1, nslip+nscrew
                  if (abs(drhoGND(is))<gndthreshold) then
                      drhoGND(is) = 0.
                  end if
              end do
!
!
!             Assign the GND value
              statev_gnd(ie,ip,1:nslip+nscrew) =
     + statev_gnd_t(ie,ip,1:nslip+nscrew) + drhoGND(1:nslip+nscrew)
!
!
!
          end do
!
!
!
      end do
!
!
!
      return
      end subroutine gndmodel1
!
!
!
!
!
!     GND model-2
!     Gurtin's measure for incompatibility: (curl x Fp)^T * Fp^T 
!     GND calculation using curl of (Fp) together with L2 approximation
!     Cumulative (total) calculation of GNDs
!     Shutting down non-active slip systems followed by singular value decompostion
!     Proposed by Chris Hardie 
      subroutine gndmodel2
      use userinputs, only : maxnslip,
     + gndhomogenization, gndthreshold  
      use globalvariables, only : numel, numdim, numpt, nnpel,
     + materialid, numslip_all, numscrew_all, screw_all,
     + slip2screw_all, burgerv_all, dirc_0_all, trac_0_all, gradip2ip, 
     + eijk, I3, statev_curvature, statev_Lambda, statev_gammasum,
     + statev_Fp, statev_gmatinv_0, statev_gnd, statev_gnd_t
      use utilities, only: matvec9
      implicit none
!     Local variables
      integer matid, nslip, nscrew
!     Some variables are allotable since material type can vary hence
!     the NUMBER OF SLIP SYSTEMS can alter within the mesh
      real(8) :: burgerv(maxnslip)
      real(8) :: gammasum(maxnslip)
      integer :: screw(maxnslip)
      real(8) :: slip2screw(maxnslip,maxnslip)
      real(8) :: dirc_0(maxnslip,3)
      real(8) :: trac_0(maxnslip,3)
!
      real(8) :: dirs(maxnslip,3)
      real(8) :: tras(maxnslip*2,3)
      real(8) :: Bmat(maxnslip*2,9)
      real(8) :: rhoGND(maxnslip*2)
      real(8) :: drhoGND(maxnslip*2)
!
      real(8) :: grad_invN(3,numpt), grad(3,3,3)
      real(8) :: Lambda(3,3), sum, Lambda_vec(9)
      real(8) :: kappa(3,3), kappa_vec(9), trace
      
      real(8) :: Fp_ip(numpt,3,3), Fp(3,3)
      real(8) :: gmatinv(3,3)

      integer :: i, j, k, l
      integer :: ie, ip, is
!
!
!
!
!     Loop through the elements
      do ie=1,numel
!
!         Reset arrays
          burgerv=0.; screw=0
          dirc_0=0.; trac_0=0.
          dirs=0.; tras=0.
          slip2screw=0.
!
!
!
!
!         Assume the same material for al the Gaussian points of an element
          matid = materialid(ie,1)              
!
!         Number of slip systems
          nslip = numslip_all(matid)
!
!         Number of screw systems
          nscrew = numscrew_all(matid)
!
!         Screw systems
          screw(1:nscrew) = screw_all(matid,1:nscrew)
!
!         Slip to screw mapping
          slip2screw(1:nscrew,1:nslip) = 
     + slip2screw_all(matid,1:nscrew,1:nslip)
!
!         Burgers vector
          burgerv(1:nslip) =  burgerv_all(matid,1:nslip)
!
!         undeformed slip direction
          dirc_0(1:nslip,1:3) = dirc_0_all(matid,1:nslip,1:3)
!
!         undeformed line direction
          trac_0(1:nslip,1:3) = trac_0_all(matid,1:nslip,1:3)
!
!
!         Store Fp for each IP
!         Plastic part of the deformation gradient
          Fp_ip = statev_Fp(ie,1:numpt,1:3,1:3)
!
!
!         Calculate the gradient of Fp
!         Calculate the gradients using gradient operator
          do ip = 1, numpt
!
!
!             Reset arrays
              Bmat=0.; rhoGND=0.; gammasum=0.              
!
!             Crystal to sample transformation matrix
              gmatinv = statev_gmatinv_0(ie,ip,:,:)
!
!
!             Cumulative slip per slip system
              gammasum(1:nslip) = statev_gammasum(ie,ip,1:nslip)
!
!
              do is = 1, nslip
!                 Transform slip directions to sample reference
                  dirs(is,:) = matmul(gmatinv, dirc_0(is,:))
!
!                 Transform line directions to sample reference
                  tras(is,:) = matmul(gmatinv, trac_0(is,:))
!
              end do
!
!
!             Calculate Bmatrix - using singular value decomposition
!             Arsenlis, A. and Parks, D.M., 1999. Acta materialia, 47(5), pp.1597-1611.
              call calculateBmatPINV(nslip, nscrew,
     + screw(1:nscrew), slip2screw(1:nscrew,1:nslip),
     + dirs(1:nslip,:), tras(1:nslip,:), burgerv(1:nslip),
     + gammasum(1:nslip), Bmat(1:nslip+nscrew,1:9))
!
!
!             Use gradients per integration point
              if (gndhomogenization == 0) then
!
                  grad_invN = gradip2ip(ie,ip,:,:)
! 
!             Use the gradient at the element center
              else
!
                  grad_invN = gradip2ip(ie,numpt+1,:,:)
!
              end if
!
!
!
!
!             grad(k,i,j) = Fp(i,j,k)
              do k = 1,3
                  do j = 1,3
                      do i = 1,3
                          grad(k,i,j) = 
     + dot_product(grad_invN(k,:),Fp_ip(:,i,j))
                      end do
                  end do
              end do
!
!
!             calculate curl
!             NOTE THE NEGATIVE SIGN AND TRANSPOSE ARE MISSING IN THE ORIGINAL REFERENCE
!             lambda(l,k) = -eijk(i,j,k) * Fp(l,j,i)
!             index "i" refers to the gradient direction
              do k = 1,3
                  do l = 1,3
                      sum = 0.
                      do i = 1, 3
                          do j = 1, 3
                              sum = sum - eijk(i,j,k)*grad(i,l,j)
!!                             earlier version (no "-" sign)
!                              sum = sum + eijk(i,j,k)*grad(i,l,j)
                          end do
                      end do
                      Lambda(l,k) = sum
                  end do
              end do
!
!
!             Plastic part of the deformation gradient
              Fp = Fp_ip(ip,:,:)
!
!
!             (curlxFp)^T * Fp^T
!             Curl post-multiplies with Fp^T
              Lambda =  matmul(Lambda, transpose(Fp))
!
!
!
!             Vectorize the incompatibility
              call matvec9(Lambda,Lambda_vec)
!
!
!             Assign the Incompatibility
              statev_Lambda(ie,ip,:) = Lambda_vec
!
!
!             Trace
              trace = Lambda(1,1) + Lambda(2,2) + Lambda(3,3)
!
!             Calculate curvature (negative transpose + trace term)
              kappa = -transpose(Lambda) + I3 * trace / 2.
!
!             Convert curvature to vector
              call matvec9(kappa,kappa_vec)
!
!
!             Assign curvature
              statev_curvature(ie,ip,1:9) = kappa_vec(1:9)
!
!
!
!
!             Reset arrays
              rhoGND=0.; drhoGND=0.
!
!             Compute the dislocation densities using L2 minimization
              rhoGND(1:nslip+nscrew) = 
     + matmul(Bmat(1:nslip+nscrew,1:9),Lambda_vec)
!
!
!             Calculate the increment of GNDs
              drhoGND(1:nslip+nscrew) = 
     + rhoGND(1:nslip+nscrew) - statev_gnd_t(ie,ip,1:nslip+nscrew)
!
!
!             Check for a threshold
              do is = 1, nslip+nscrew
                  if (abs(drhoGND(is))<gndthreshold) then
                      drhoGND(is) = 0.
                  end if
              end do    
!
!
!             Assign the GND value
              statev_gnd(ie,ip,1:nslip+nscrew) = 
     + statev_gnd_t(ie,ip,1:nslip+nscrew) + drhoGND(1:nslip+nscrew)
!
!
!
          end do
!
!
!
!
      end do
!
!
!
      return
      end subroutine gndmodel2
!
!
!
!
!     GND model-3
!     Based on original rate formulation of GNDs
!     GND calculation using curl of (n^a*Fp*gammadot^a)
!     Incremental or rate form
!     followed by direct projections
!     Dai, H., 1997. Doctoral dissertation, Massachusetts Institute of Technology.
      subroutine gndmodel3(dt)
      use userinputs, only : maxnslip,
     + gndhomogenization, gndthreshold
      use globalvariables, only : numel, numdim, numpt, nnpel,
     + materialid, numslip_all, numscrew_all, eijk, I3,
     + dirc_0_all, trac_0_all, norc_0_all, burgerv_all,
     + slip2screw_all, gradip2ip, statev_curvature,
     + statev_gammadot, statev_gmatinv_0, statev_Fp,
     + statev_gnd_t, statev_gnd,
     + statev_Lambda, statev_Lambda_t
      use utilities, only: matvec9, vecmat9
      implicit none
!     time increment
      real(8), intent(in) :: dt
!     Local variables
      real(8) :: gdot, Fp(3,3), gmatinv(3,3)
      integer matid, nslip, nscrew
!     Some variables are allotable since material type can vary hence
!     the NUMBER OF SLIP SYSTEMS can alter within the mesh
      real(8) :: dirc_0(maxnslip,3)
      real(8) :: norc_0(maxnslip,3)
      real(8) :: trac_0(maxnslip,3)
      real(8) :: slip2screw(maxnslip,maxnslip)
      real(8) :: drhoGNDe(numpt,maxnslip)
      real(8) :: drhoGNDs(numpt,maxnslip)
      real(8) :: drhoGND(maxnslip*2)
      real(8) :: burgerv(maxnslip)
      real(8) :: vec_ip(numpt,3)
      real(8) :: n_a(3)
      real(8) :: tras_ip(numpt,3), dirs_ip(numpt,3)
!     Incompatibility
      real(8) :: dLambda_ip(numpt,3,3)
      real(8) :: dLambda(3,3), dLambda_vec(9)
      real(8) :: Lambda_vec(9), Lambda(3,3)
      real(8) :: kappa(3,3), kappa_vec(9)
!     Overall gradient mapping
      real(8) :: grad_invN(3,numpt), grad(3,3)
      real(8) :: lambdadot_a(3), sum, trace
      integer :: i, j, k, l
      integer :: ie, ip, is
!
!
!
!     Elemental calculation
!     For each element
      do ie = 1, numel
!
!         Reset arrays
          dirc_0=0.; norc_0=0.; trac_0=0.
          burgerv=0.; slip2screw=0.
!
!         Assume the same material for al the Gaussian points of an element
          matid = materialid(ie,1)
!
!         Number of slip systems
          nslip = numslip_all(matid)
!
!         Number of screw systems
          nscrew = numscrew_all(matid)
!
!         Burgers vector
          burgerv(1:nslip) =  burgerv_all(matid,1:nslip)
!
!         undeformed slip direction
          dirc_0(1:nslip,:) = dirc_0_all(matid,1:nslip,:)
!
!         undeformed slip plane normal
          norc_0(1:nslip,:) = norc_0_all(matid,1:nslip,:)
!
!         undeformed transverse direction
          trac_0(1:nslip,:) = trac_0_all(matid,1:nslip,:)
!
!         slip to screw system mapping
          slip2screw(1:nscrew,1:nslip) =
     + slip2screw_all(matid,1:nscrew,1:nslip)
!
!
!         Set incompatibility increment to zero
          dLambda_ip = 0.
!
!         Reset arrays
          drhoGNDe=0.; drhoGNDs=0.; drhoGND=0.
!
!
!         For each slip system
          do is = 1, nslip
!
!
!             Calculate the vector of known quantities
!
!             For each integration point
              do ip = 1, numpt
!
!
!                 Slip rate
                  gdot = statev_gammadot(ie,ip,is)
!
!                 Plastic part of the deformation gradient
                  Fp = statev_Fp(ie,ip,:,:)
!
!
!                 Crystal to sample transformation matrix
                  gmatinv = statev_gmatinv_0(ie,ip,:,:)
!
!                 Slip plane normal in sample reference
                  n_a = matmul(gmatinv, norc_0(is,:))
!
!                 transverse directions (store) in sample reference
                  tras_ip(ip,:) = matmul(gmatinv, trac_0(is,:))
!                 
!                 slip directions (store) in sample reference
                  dirs_ip(ip,:) = matmul(gmatinv, dirc_0(is,:))
!                 
!                 Result = vector
                  vec_ip(ip,:) = -matmul(transpose(Fp),n_a) * gdot
!
!
!             
!
              end do
!
!
!
!
!             Calculate the gradients using gradient operator
              do ip = 1, numpt
!
!                 Use gradients per integration point
                  if (gndhomogenization == 0) then
!
                      grad_invN = gradip2ip(ie,ip,:,:)
!
!                 Use the gradient at the element center
                  else
!
                      grad_invN = gradip2ip(ie,numpt+1,:,:)
!
                  end if
!
!
!
!                 grad(k,i) = - n^a(j) * Fp(ji,k) * gdot^a
!                 lambdadot^a(l) = eijk(k,i,l) * grad(k,i)
!                 index "k" refers to the gradient direction
                  grad = matmul(grad_invN,vec_ip)
!
!                 calculate curl
                  do l = 1, 3
                      sum = 0.
                      do k = 1, 3
                          do i = 1, 3
                              sum = sum + eijk(k,i,l)*grad(k,i)
                          end do
                      end do
                      lambdadot_a(l) = sum
                  end do                  
!
!
!                 calculate incompatibility increment
                  do i = 1, 3
                      do j = 1, 3
                          dLambda_ip(ip,i,j) =
     + dLambda_ip(ip,i,j) + dirs_ip(ip,i)*lambdadot_a(j)*dt
                      end do
                  end do                  
!
!
!
!                 gnd increment - edge dislocations                 
                  drhoGNDe(ip,is) = 1. / burgerv(is) *
     + dot_product(lambdadot_a, tras_ip(ip,:)) * dt
!
!
!
!                 gnd increment - screw dislocations                 
                  drhoGNDs(ip,is) = 1. / burgerv(is) *
     + dot_product(lambdadot_a, dirs_ip(ip,:)) * dt
!
!
!
!
!
!             end of IP loop                  
              end do
!
!
!
!
!
!         end of slip system loop          
          end do
!
!
!         Calculate overall GNDs
          do ip = 1, numpt
!
!             Use the projection to find the screw dislocation density
!
!             Edge dislocations
              drhoGND(1:nslip) = drhoGNDe(ip,1:nslip)
!
!             Screw dislocations
              drhoGND(nslip+1:nslip+nscrew) = matmul(
     + slip2screw(1:nscrew,1:nslip),drhoGNDs(ip,1:nslip))
!
!
!             Check for a threshold
              do is =1,nslip+nscrew
                  if (abs(drhoGND(is))<gndthreshold) then
                      drhoGND(is) = 0.
                  end if
              end do
!
!             Assign the overall GND density
              statev_gnd(ie,ip,1:nslip+nscrew) =
     + statev_gnd_t(ie,ip,1:nslip+nscrew) + 
     + drhoGND(1:nslip+nscrew)
!
!
!
!             Incompatibility increment dyadic
              dLambda = dLambda_ip(ip,:,:)
!
!             Vectorize
              call matvec9(dLambda,dLambda_vec)
!
!             Add the result to the incompatibility
              Lambda_vec=statev_Lambda_t(ie,ip,:) + dLambda_vec(:)
!
!             Assign it to state variable
              statev_Lambda(ie,ip,:) = Lambda_vec
!
!
!             Convert to 3x3 matrix
              call vecmat9(Lambda_vec,Lambda)
!
!             Trace
              trace = Lambda(1,1) + Lambda(2,2) + Lambda(3,3)
!
!             Calculate curvature (negative transpose + trace term)
              kappa = -transpose(Lambda) + I3 * trace / 2.
!
!             Convert curvature to vector
              call matvec9(kappa,kappa_vec)          
!
!
              statev_curvature(ie,ip,1:9) = kappa_vec(1:9)
!
!
!
!
!         end of IP loop                  
          end do
!
!
!
!     end of element loop
      end do
!
!
!
      return
      end subroutine gndmodel3
!
!
!
!
!
!
!
!
!
!     GND model-4
!     GND calculation using slip gradients
!     Incremental or rate form
!     followed by direct projections
!     Gerken, J.M. and Dawson, P.R., 2008. Journal of the Mechanics and Physics of Solids, 56(4), pp.1651-1672.
      subroutine gndmodel4(dt)
      use userinputs, only : maxnslip,
     + gndhomogenization, gndthreshold
      use globalvariables, only : I3, 
     + numdim, numel, numpt, nnpel,
     + materialid, numslip_all, numscrew_all,
     + dirc_0_all, trac_0_all, burgerv_all,
     + slip2screw_all, gradip2ip, statev_curvature,
     + statev_gammadot, statev_gmatinv_0,
     + statev_gnd_t, statev_gnd,
     + statev_Lambda, statev_Lambda_t
      use utilities, only: matvec9, vecmat9
      implicit none
!     time increment
      real(8), intent(in) :: dt
!     Local variables
      real(8) :: gdot(numpt), gmatinv(3,3)
      integer matid, nslip, nscrew
!     Some variables are allotable since material type can vary hence
!     the NUMBER OF SLIP SYSTEMS can alter within the mesh
      real(8) :: dirc_0(maxnslip,3)
      real(8) :: trac_0(maxnslip,3)
      real(8) :: slip2screw(maxnslip,maxnslip)
      real(8) :: drhoGNDe(numpt,maxnslip)
      real(8) :: drhoGNDs(numpt,maxnslip)
      real(8) :: drhoGND(maxnslip*2)
      real(8) :: burgerv(maxnslip)
      real(8) :: tras_ip(numpt,3), dirs_ip(numpt,3)
      real(8) :: dot_t, dot_s
!     Incompatibility
      real(8) :: dLambda(3,3), dLambda_vec(9)
      real(8) :: Lambda_vec(9), Lambda(3,3), trace
      real(8) :: kappa(3,3), kappa_vec(9)
      real(8) :: dLambda_ip(numpt,3,3)
!     Overall gradient mapping
      real(8) :: grad_invN(3,numpt), grad(3)
      integer :: i, j, k, q
      integer :: ie, ip, is   
!
!
!
!
!     Elemental calculation
!     For each element
      do ie = 1, numel
!
!
!         Reset arrays
          dirc_0=0.; trac_0=0.
          burgerv=0.; slip2screw=0.
!
!         Assume the same material for al the Gaussian points of an element
          matid = materialid(ie,1)
!
!         Number of slip systems
          nslip = numslip_all(matid)
!
!         Number of screw systems
          nscrew = numscrew_all(matid) 
!
!         Burgers vector
          burgerv(1:nslip) =  burgerv_all(matid,1:nslip)
!
!         undeformed slip direction
          dirc_0(1:nslip,:) = dirc_0_all(matid,1:nslip,:)
!
!         undeformed transverse direction
          trac_0(1:nslip,:) = trac_0_all(matid,1:nslip,:)
!
!         slip to screw system mapping
          slip2screw(1:nscrew,1:nslip) =
     + slip2screw_all(matid,1:nscrew,1:nslip)
!
!
!
!         Set incompatibility increment to zero
          dLambda_ip = 0.
!
!         Reset arrays
          drhoGNDe=0.; drhoGNDs=0.
!
!
!         For each slip system
          do is = 1, nslip
!
!
!             Calculate the vector of known quantities
!
!             For each integration point
              do ip = 1, numpt
!
!
!                 Slip rate
                  gdot(ip) = statev_gammadot(ie,ip,is)
!
!
!
!
!                 Crystal to sample transformation matrix
                  gmatinv = statev_gmatinv_0(ie,ip,:,:)
!
!
!
!                 transverse directions (store)
                  tras_ip(ip,:) = matmul(gmatinv, trac_0(is,:))
!
!
!                 slip directions (store)
                  dirs_ip(ip,:) = matmul(gmatinv, dirc_0(is,:))
!
!
              end do
!
!
!             Calculate the gradients using gradient operator
              do ip = 1, numpt
!
!                 Use gradients per integration point
                  if (gndhomogenization == 0) then
!
                      grad_invN = gradip2ip(ie,ip,:,:)
!
!                 Use the gradient at the element center
                  else
!
                      grad_invN = gradip2ip(ie,numpt+1,:,:)
!
                  end if
!
!
!                 grad(k) = gdot^a
!                 lambdadot^a(k) = grad(k,1:numpt) * gdot^a(1:numpt)
!                 index "k" refers to the gradient direction
                  grad = matmul(grad_invN,gdot)
!
!
!
!
!                 gnd increment - screw dislocations
                  dot_t = dot_product(grad, tras_ip(ip,:))
!
!                 
                  drhoGNDs(ip,is) = 1. / burgerv(is) * dot_t * dt
!
!
!                 gnd increment - edge dislocations
                  dot_s = dot_product(grad, dirs_ip(ip,:))
!
!                 Note the negative sign for screws
                  drhoGNDe(ip,is) = -1. / burgerv(is) * dot_s * dt
!
!
!
!
!                 Compute incompatibility dyadic
!
                  do i = 1, 3
                      do j = 1, 3
!
                          dLambda_ip(ip,i,j) = dLambda_ip(ip,i,j)
     + - dot_s * dirs_ip(ip,i) * tras_ip(ip,j) * dt
     + + dot_t * dirs_ip(ip,i) * dirs_ip(ip,j) * dt
!
!                         
                      end do
                  end do
!
!
!
!
!
!             end of IP loop
              end do
!
!
!
!
!
!         end of slip system loop
          end do
!
!
!         Calculate overall GNDs
          do ip = 1, numpt
!
              drhoGND=0.
!
!             Use the projection to find the screw dislocation density
!
!             Edge dislocations
              drhoGND(1:nslip) = drhoGNDe(ip,1:nslip)
!
!             Screw dislocations
              drhoGND(nslip+1:nslip+nscrew) = matmul(
     + slip2screw(1:nscrew,1:nslip),drhoGNDs(ip,1:nslip))
!
!             Check for a threshold
              do is =1,nslip+nscrew
                  if (abs(drhoGND(is))<gndthreshold) then
                      drhoGND(is) = 0.
                  end if
              end do              
!
!
!             Assign the overall GND density
              statev_gnd(ie,ip,1:nslip+nscrew) =
     + statev_gnd_t(ie,ip,1:nslip+nscrew) + 
     + drhoGND(1:nslip+nscrew)
!
!
!
!
!
!
!
!
!             Incompatibility increment dyadic
              dLambda = dLambda_ip(ip,:,:)
!
!             Vectorize
              call matvec9(dLambda,dLambda_vec)
!
!             Add the result to the incompatibility
              Lambda_vec=statev_Lambda_t(ie,ip,:) + dLambda_vec(:)
!
!             Assign it to state variable
              statev_Lambda(ie,ip,:) = Lambda_vec
!
!             Convert to 3x3 matrix
              call vecmat9(Lambda_vec,Lambda)
!
!             Trace
              trace = Lambda(1,1) + Lambda(2,2) + Lambda(3,3)
!
!             Calculate curvature (negative transpose + trace term)
              kappa = -transpose(Lambda) + I3 * trace / 2.
!
!             Convert curvature to vector
              call matvec9(kappa,kappa_vec)
!
!
              statev_curvature(ie,ip,1:9) = kappa_vec(1:9)
!
!
!
!             
!         end of IP loop
          end do
!
!
!
!     end of element loop
      end do
!
!
!
      return
      end subroutine gndmodel4
!
!
!
!!     GND model-1_old
!!     GND calculation using curl of (Fp) together with L2 approximation
!!     Cumulative (total) calculation of GNDs
!!     Kuksenko, V., Roberts, S. and Tarleton, E., 2019. International Journal of Plasticity, 116, pp.62-80.
!      subroutine gndmodel1_old
!      use userinputs, only : numel, maxnslip,
!     + gndhomogenization, gndthreshold
!      use globalvariables, only : numdim, numpt, nnpel,
!     + phaseid, numslip_all, numscrew_all, screw_all,
!     + burgerv_all, dirc_0_all, trac_0_all, gradip2ip,
!     + I3, eijk, statev_curvature, statev_Lambda,
!     + statev_Fp, statev_gmatinv_0, statev_gnd, statev_gnd_t
!      use utilities, only: matvec9
!      implicit none
!!     Local variables
!      integer matid, nslip, nscrew
!!     Some variables are allotable since material type can vary hence
!!     the NUMBER OF SLIP SYSTEMS can alter within the mesh
!      real(8) :: burgerv(maxnslip)
!      integer :: screw(maxnslip)
!      real(8) :: dirc_0(maxnslip,3)
!      real(8) :: trac_0(maxnslip,3)
!!
!      real(8) :: dirs(maxnslip,3)
!      real(8) :: tras(maxnslip*2,3)
!      real(8) :: Bmat(maxnslip*2,9)
!      real(8) :: rhoGND(maxnslip*2)
!      real(8) :: drhoGND(maxnslip*2)
!!
!      real(8) :: grad_invN(3,numpt), grad(3,3,3)
!      real(8) :: Lambda(3,3), sum, Lambda_vec(9)
!      real(8) :: kappa(3,3), kappa_vec(9), trace
!      real(8) :: Fp_ip(numpt,3,3)
!      real(8) :: gmatinv(3,3)
!
!      integer :: i, j, k, l
!      integer :: ie, ip, is
!!     
!!
!!
!!
!!     Loop through the elements
!      do ie=1,numel
!!
!!         Reset arrays
!          burgerv=0.;screw=0
!          dirc_0=0.; trac_0=0.
!          dirs=0.; tras=0.
!!
!!
!!         Assume the same material for al the Gaussian points of an element
!          matid = phaseid(ie,1)
!!
!!         Number of slip systems
!          nslip = numslip_all(matid)
!!
!!         Number of screw systems
!          nscrew = numscrew_all(matid)
!!
!!         Screw systems
!          screw = screw_all(matid,1:nscrew)
!!
!!         Burgers vector
!          burgerv(1:nslip) =  burgerv_all(matid,1:nslip)
!!
!!         undeformed slip direction
!          dirc_0(1:nslip,1:3) = dirc_0_all(matid,1:nslip,1:3)
!!
!!         undeformed line direction
!          trac_0(1:nslip,1:3) = trac_0_all(matid,1:nslip,1:3)
!!
!!
!!
!!         Store Fp for each IP
!!         Plastic part of the deformation gradient
!          Fp_ip = statev_Fp(ie,1:numpt,1:3,1:3)
!!
!!
!!         Calculate the gradient of Fp
!!         Calculate the gradients using gradient operator
!          do ip = 1, numpt
!!
!!             Reset arrays
!              Bmat=0.; rhoGND=0.
!!
!!             Crystal to sample transformation matrix
!              gmatinv = statev_gmatinv_0(ie,ip,:,:)              
!!
!!
!              do is = 1, nslip
!!                 Transform slip directions to sample reference
!                  dirs(is,:) = matmul(gmatinv, dirc_0(is,:))
!!
!!                 Transform line directions to sample reference
!                  tras(is,:) = matmul(gmatinv, trac_0(is,:))
!!
!!
!              end do
!!
!!             Calculate Bmatrix - using singular value decomposition
!!             Arsenlis, A. and Parks, D.M., 1999. Acta materialia, 47(5), pp.1597-1611.
!              call calculateBmat(nslip,nscrew,screw,
!     + dirs(1:nslip,:),tras(1:nslip,:),burgerv(1:nslip),
!     + Bmat(1:nslip+nscrew,1:9))
!!
!!
!!             Use gradients per integration point
!              if (gndhomogenization == 0) then
!!
!                  grad_invN = gradip2ip(ie,ip,:,:)
!!
!!             Use the gradient at the element center
!              else
!!
!                  grad_invN = gradip2ip(ie,numpt+1,:,:)
!!
!              end if
!!
!!
!!
!!
!!             Calculate gradient
!              do k = 1,3
!                  do j = 1,3
!                      do i = 1,3
!                          sum = 0.
!                          do l = 1, numpt
!                              sum = sum + grad_invN(k,l)*Fp_ip(l,i,j)
!                          end do
!                          grad(i,j,k) = sum
!                      end do
!                  end do
!              end do
!!
!!
!!             calculate curl
!!             NOTE THE NEGATIVE SIGN AND TRANSPOSE ARE MISSING IN THE ORIGINAL REFERENCE
!!             lambda(l,k) = -eijk(i,j,k) * Fp(l,j,i)
!!             index "i" refers to the gradient direction
!              do k = 1,3
!                  do l = 1,3
!                      sum = 0.
!                      do i = 1, 3
!                          do j = 1, 3
!                              sum = sum - eijk(i,j,k)*grad(l,j,i)
!!!                             earlier version (no "-" sign)
!!                              sum = sum + eijk(i,j,k)*grad(i,l,j)
!                          end do
!                      end do
!                      Lambda(l,k) = sum
!                  end do
!              end do
!!
!!             Vectorize the incompatibility
!              call matvec9(Lambda,Lambda_vec)
!!
!!
!!             Assign the Incompatibility
!              statev_Lambda(ie,ip,:) = Lambda_vec
!!
!!             Trace
!              trace = Lambda_vec(1) + Lambda_vec(2) + Lambda_vec(3)
!!
!!             Calculate curvature (negative transpose + trace term)
!              kappa = -transpose(Lambda) + I3 * trace / 2.
!!
!!             Convert curvature to vector
!              call matvec9(kappa,kappa_vec)
!!
!!             Assign curvature
!              statev_curvature(ie,ip,1:9) = kappa_vec(1:9)
!!
!!             Reset arrays
!              rhoGND=0.; drhoGND=0.
!!             Compute the dislocation densities using L2 minimization
!              rhoGND(1:nslip+nscrew) =
!     + matmul(Bmat(1:nslip+nscrew,1:9),Lambda_vec)
!!
!!
!!             Calculate the increment of GNDs
!              drhoGND(1:nslip+nscrew) =
!     + rhoGND(1:nslip+nscrew) - statev_gnd_t(ie,ip,1:nslip+nscrew)
!!
!!
!!             Check for a threshold
!              do is = 1, nslip+nscrew
!                  if (abs(drhoGND(is))<gndthreshold) then
!                      drhoGND(is) = 0.
!                  end if
!              end do
!!
!!
!!             Assign the GND value
!              statev_gnd(ie,ip,1:nslip+nscrew) =
!     + statev_gnd_t(ie,ip,1:nslip+nscrew) + drhoGND(1:nslip+nscrew)
!!
!!
!!
!          end do
!!      
!!
!!
!!
!!
!!
!      end do
!!
!!
!!
!      return
!      end subroutine gndmodel1_old
!
!
!
!
!
!     Calculate Bmatrix using KKT condition
      subroutine calculateBmatKKT(nslip,nscrew,screw,dirs,tras,burgerv,
     + BmatKKT)
      use utilities, only : nolapinverse
      use errors, only: error
      implicit none
!     number of slip systems
      integer, intent(in) :: nslip
!     number of screw systems
      integer, intent(in) :: nscrew
!     number of screw systems
      integer, dimension(nscrew), intent(in) :: screw
!     slip direction in sample reference
      real(8), dimension(nslip,3), intent(in) :: dirs
!     transverse (to slip) direction in sample reference
      real(8), dimension(nslip,3), intent(in) :: tras
!     Burgers vector
      real(8), dimension(nslip), intent(in) :: burgerv
!     L2 mapping at the deformed configuration
      real(8), dimension(nslip+nscrew+9,nslip+nscrew+9),
     + intent(out) :: BmatKKT
!
!     Local variables used within this subroutine
      real(8) :: Amat(9,nslip+nscrew)
!     KKT matrix in L2 method
      real(8) :: KKTmat(nslip+nscrew+9,nslip+nscrew+9)
      real(8) :: s(3), l(3), b
      integer i, is
!
!     Construct Nye's tensor
      Amat=0.
!     Loop through dislocation configurations
      do i = 1, nslip+nscrew
!
!         Slip direction
!         For screws
          if (i>nslip) then
              is = screw(i-nslip)
              s = dirs(is,:)
              l = dirs(is,:)
              b = burgerv(is)
!         For edges
          else
              s = dirs(i,:)
              l = tras(i,:)
              b = burgerv(i)
          end if
!
!         The ordering of Amat is as follows:
!         11-12-13-21-22-23-31-32-33
          Amat(1,i)=s(1)*l(1)*b; Amat(2,i)=s(1)*l(2)*b
          Amat(3,i)=s(1)*l(3)*b; Amat(4,i)=s(2)*l(1)*b
          Amat(5,i)=s(2)*l(2)*b; Amat(6,i)=s(2)*l(3)*b
          Amat(7,i)=s(3)*l(1)*b; Amat(8,i)=s(3)*l(2)*b
          Amat(9,i)=s(3)*l(3)*b
!
      end do
!
!     L2 solution by KKT condition
!     Assign zero initially
      KKTmat = 0.
!     Assign the identity part
      do i = 1, nslip+nscrew
          KKTmat(i,i)=2.
      end do
!
!     Assign A^T - upper right part
      KKTmat(1:nslip+nscrew,nslip+nscrew+1:nslip+nscrew+9) =
     + transpose(Amat)
!
!     Assign A - lower left part
      KKTmat(nslip+nscrew+1:nslip+nscrew+9,1:nslip+nscrew) =
     + Amat
!
!     Take the inverse
      call nolapinverse(KKTmat,BmatKKT,nslip+nscrew+9)
!     
!     Set to zero if not invertible
      if(any(BmatKKT /= BmatKKT)) then
!         Set the outputs to zero initially
          BmatKKT = 0.
!         error message in .dat file
          call error(10)
      end if            
!
!
!
!
!
!
      return
      end subroutine calculateBmatKKT
!
!
!
!
!     Calculate Bmatrix using singular value decomposition
      subroutine calculateBmat(nslip,nscrew,screw,dirs,tras,burgerv,
     + Bmat)
      use utilities, only : nolapinverse
      use errors, only: error
      implicit none
!     number of slip systems
      integer, intent(in) :: nslip
!     number of screw systems
      integer, intent(in) :: nscrew
!     screw systems
      integer, dimension(nscrew), intent(in) :: screw
!     slip direction in sample reference
      real(8), dimension(nslip,3), intent(in) :: dirs
!     transverse (to slip) direction in sample reference
      real(8), dimension(nslip,3), intent(in) :: tras
!     Burgers vector
      real(8), dimension(nslip), intent(in) :: burgerv
!     Rank Deficit B-matrix
      real(8), dimension(nslip+nscrew,9), intent(out) :: Bmat
!
!     Local variables used within this subroutine
!     Note A^T*A is rank deficit
      real(8) :: Amat(9,nslip+nscrew)
      real(8) :: AmatAmatT(9,9)
      real(8) :: invAmatAmatT(9,9)
      real(8) :: s(3), l(3), b
      integer :: i, j, is
!
!     Construct Nye's tensor
      Amat=0.
!     Loop through dislocation configurations
      do i = 1, nslip+nscrew
!
!
!         Slip direction
!         For screws
          if (i>nslip) then
              is = screw(i-nslip)
              s = dirs(is,:)
              l = dirs(is,:)
              b = burgerv(is)
!         For edges
          else
              s = dirs(i,:)
              l = tras(i,:)
              b = burgerv(i)
          end if
!
!         The ordering of Amat is as follows:
!         11-12-13-21-22-23-31-32-33
          Amat(1,i)=s(1)*l(1)*b; Amat(2,i)=s(1)*l(2)*b
          Amat(3,i)=s(1)*l(3)*b; Amat(4,i)=s(2)*l(1)*b
          Amat(5,i)=s(2)*l(2)*b; Amat(6,i)=s(2)*l(3)*b
          Amat(7,i)=s(3)*l(1)*b; Amat(8,i)=s(3)*l(2)*b
          Amat(9,i)=s(3)*l(3)*b
!
      end do
!
!
!
!
!
!     Other solution (former method)
!     -----------------------------------------------
!     This is preserved to check its validity
!     Right pseudo inverse is used
!     This is exactly the same as the inversion by singular value decomposition
!     Compute the L2 coefficient
      AmatAmatT = matmul(Amat,transpose(Amat))     
!
!     Take the inverse
      call nolapinverse(AmatAmatT,invAmatAmatT,9)
!
!      
!     Compute Bmat
      Bmat = matmul(transpose(Amat),invAmatAmatT)    
!
!     Set to zero if not invertible
      if(any(Bmat /= Bmat)) then
!         Set the outputs to zero initially
          Bmat = 0.
!         error message in .dat file
          call error(10)
      end if        
!
!
!
      return
      end subroutine calculateBmat
!
!
!     Calculate Bmatrix using singular value decomposition
      subroutine calculateBmatPINV(nslip,nscrew,screw,slip2screw,
     + dirs,tras,burgerv,gammadot,BmatPINV)
      use userinputs, only: slipthreshold
      use utilities, only: svdgeninverse
      use errors, only: error
      implicit none
!     number of slip systems
      integer, intent(in) :: nslip
!     number of screw systems
      integer, intent(in) :: nscrew
!     screw systems
      integer, dimension(nscrew), intent(in) :: screw
!     slip to screw mapping
      real(8), dimension(nscrew,nslip), intent(in) :: slip2screw
!     slip direction in sample reference
      real(8), dimension(nslip,3), intent(in) :: dirs
!     transverse (to slip) direction in sample reference
      real(8), dimension(nslip,3), intent(in) :: tras
!     Burgers vector
      real(8), dimension(nslip), intent(in) :: burgerv
!     Cumulative slip per slip system
      real(8), dimension(nslip), intent(in) :: gammadot
!     Rank Deficit B-matrix
      real(8), dimension(nslip+nscrew,9), intent(out) :: BmatPINV
!
!     Local variables used within this subroutine
!     Note A^T*A is rank deficit
      real(8) :: Amat(9,nslip+nscrew)
      real(8) :: gdot_screw(nscrew)
!      real(8) :: AmatTAmat(nslip+nscrew,nslip+nscrew)
!      real(8) :: invAmatTAmat(nslip+nscrew,nslip+nscrew)
      real(8) :: s(3), l(3), b
      integer :: i, j, is, err
!
!     Construct Nye's tensor
      Amat=0.
!     Loop through dislocation configurations
      do i = 1, nslip+nscrew
!
!
!         Slip direction
!         For screws
          if (i>nslip) then
              is = screw(i-nslip)
              s = dirs(is,:)
              l = dirs(is,:)
              b = burgerv(is)
!         For edges
          else
              s = dirs(i,:)
              l = tras(i,:)
              b = burgerv(i)
          end if
!
!
!         The ordering of Amat is as follows:
!         11-12-13-21-22-23-31-32-33
          Amat(1,i)=s(1)*l(1)*b; Amat(2,i)=s(1)*l(2)*b
          Amat(3,i)=s(1)*l(3)*b; Amat(4,i)=s(2)*l(1)*b
          Amat(5,i)=s(2)*l(2)*b; Amat(6,i)=s(2)*l(3)*b
          Amat(7,i)=s(3)*l(1)*b; Amat(8,i)=s(3)*l(2)*b
          Amat(9,i)=s(3)*l(3)*b
      end do
!
!
!     Sum slip on each system sharing common screw types
      gdot_screw= matmul(slip2screw,gammadot)
!
!
!     Shut down the slip systems 
!     that have a cumulative slip
!     less than 10^-10

!     For edge type of dislocations
      do is = 1, nslip
!
          if (abs(gammadot(is))<slipthreshold) then
!
!             Set columns to zero - for edges
              Amat(:,is) = 0.
!
          end if
!
      end do
!
!     For screw type of dislocations
      do i = 1, nscrew
!
!          is = screw(i)
!
          if (abs(gdot_screw(i))<slipthreshold) then
!
!             Set columns to zero - for edges
              Amat(:,i+nslip) = 0.
!
          end if
!
      end do      
!
!
!
!!     Other solution (former method)
!!     -----------------------------------------------
!!     This is preserved to check its validity
!!     BmatPINV: B-mat calculation using pseudo inverse
!!     Former method gives large number due to rank deficit inversion
!!     Compute the L2 coefficient
!      AmatTAmat = matmul(transpose(Amat),Amat)
!!
!!      if (any(isnan(AmatAmatT))) write(*,*) 'NaN before'
!!     Take the inverse
!
!!     Using singular value decomposition      
!      call svdinverse(AmatTAmat,nslip+nscrew,invAmatTAmat,err)
!!      if (any(isnan(invAmatAmatT))) write(*,*) 'NaN after'
!!
!!      
!!      Compute Bmat
!      BmatPINV = matmul(invAmatTAmat,transpose(Amat))
!
!     Generalized inverse
!     Subroutine implemented by Alvaro 8-3-2023
!
      call svdgeninverse(Amat,9,nslip+nscrew,BmatPINV,err)
!
!     Check if there is an error
      if (err/=0) then
          BmatPINV = 0.
!         error message
          call error(10)
      end if
!
!
!     Set to zero if not invertible
      if(any(BmatPINV /= BmatPINV)) then
!         Set the outputs to zero initially
          BmatPINV = 0.
!         error message
          call error(10)
      end if
!
!
!
      return
      end subroutine calculateBmatPINV
!
!     
!
      end module straingradients