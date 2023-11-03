!         Else Forward Gradient Predictor scheme computes sigma0
          elseif (predictor == 1) then
!
!
              call CP_ForwardGradientPredictor(
     + matid, phaid, nslip, nscrew,
     + mattemp, Cs, gf, G12,
     + burgerv, cubicslip, caratio,
     + Fp_t, gmatinv_t, Eec_t,
     + gammadot_t, gammasum_t,
     + totgammasum_t, evmp_t, sigmarot_t,
     + forestproj, slip2screw, dirs_t, nors_t,
     + Schmidvec, Schmid_0,
     + Schmid, SchmidxSchmid,
     + smodel, sparam, cmodel, cparam,
     + imodel, iparam, hmodel, hparam,
     + bparam, sintmat1, sintmat2,
     + hintmat1, hintmat2,
     + tauceff_t, tauc_t, rhotot_t,
     + sumrhotot_t, ssdtot_t,
     + rhofor_t, forest_t, substructure_t,
     + gnd_t, ssd_t, loop_t, X_t,
     + dt, dstran, domega,
     + Fp0, gmatinv0, Eec0,
     + gammadot0, gammasum0,
     + totgammasum0, evmp0,
     + tauc0, tausolute0,
     + ssdtot0, ssd0, loop0, X0,
     + forest0, substructure0,
     + sigma0, jacobi0, cpconv0)
!
!
!
              if (cpconv0==0) then
                  sigma0 = (1.-phi)*sigma_t + phi*sigmatr
                  endif
!
!
                  
!     Fully implicit solution
      subroutine CP_Huang_fullyimplicit(matid, phaid,
     + nslip, nscrew, mattemp, Cs, gf, G12,
     + burgerv, cubicslip, caratio,
     + Fp_t, gmatinv_t, Eec_t,
     + gammadot_t, gammasum_t,
     + totgammasum_t, evmp_t, sigma_t,
     + forestproj, slip2screw, dirs_t, nors_t,
     + slipmodel, slipparam,
     + creepmodel, creepparam,
     + irradiationmodel, irradiationparam,
     + hardeningmodel, hardeningparam,
     + backstressparam, 
     + sintmat1, sintmat2,
     + hintmat1, hintmat2,
     + tauceff_t, tauc_t, rhotot_t,
     + sumrhotot_t, ssdtot_t,
     + rhofor_t, forest_t, substructure_t,
     + gnd_t, ssd_t, loop_t, X_t,
     + dt, dstran, domega,
     + Fp, gmatinv, Eec,
     + gammadot, gammasum,
     + totgammasum, evmp,
     + tauc, tausolute,
     + ssdtot, ssd, loop, X,
     + forest, substructure,
     + sigma, jacobi, cpconv)
!
      use globalvariables, only : I3, I6, smallnum
!
      use userinputs, only : theta, maxnparam, maxnloop,
     + maxniter, tolerance, backstressmodel
!
      use utilities, only : vecmat6, matvec6,
     + nolapinverse, deter3x3, inv3x3, trace3x3,
     + vecmat9, matvec9, gmatvec6
!
      use slip, only : sinhslip, doubleexpslip,
     + powerslip
!
      use creep, only : expcreep
!
      use hardening, only: hardeningrules
!
      use backstress, only: backstressmodel1
!
      use crss, only: slipresistance
!
      use errors, only : error
!
      implicit none
!
!     INPUTS
!
!     material-id
      integer, intent(in) :: matid
!     phase-id
      integer, intent(in) :: phaid
!     number of slip sytems
      integer, intent(in) :: nslip
!     number of screw sytems
      integer, intent(in) :: nscrew
!     temperature
      real(8), intent(in) :: mattemp
!     elastic compliance
      real(8), intent(in) :: Cs(6,6)
!     geometric factor
      real(8), intent(in) :: gf
!     elastic shear modulus
      real(8), intent(in) :: G12
!     Burgers vectors
      real(8), intent(in) :: burgerv(nslip)
!     flag for cubic slip systems
      integer, intent(in) :: cubicslip
!     c/a ratio for hcp crystals
      real(8), intent(in) :: caratio
!     plastic part of the deformation gradient at former time step
      real(8), intent(in) :: Fp_t(3,3)
!     Crystal to sample transformation martrix at former time step
      real(8), intent(in) :: gmatinv_t(3,3)
!     Lattice strains
      real(8), intent(in) :: Eec_t(6)
!     slip rates at the former time step
      real(8), intent(in) :: gammadot_t(nslip)
!     total slip per slip system accumulated over the time
!     at the former time step
      real(8), intent(in) :: gammasum_t(nslip)
!     overall total slip at the former time step
      real(8), intent(in) :: totgammasum_t
!     Von-Mises equivalent total plastic strain at the former time step
      real(8), intent(in) :: evmp_t
!     Cauchy stress at the former time step
      real(8), intent(in) :: sigma_t(6)
!     Forest projection for GND
      real(8), intent(in) :: forestproj(nslip,nslip+nscrew)
!     Slip to screw system mapping
      real(8), intent(in) :: slip2screw(nscrew,nslip)
!     deformed slip direction in sample reference frame
      real(8), intent(in) :: dirs_t(nslip,3)
!     deformed slip plane normal in sample reference frame
      real(8), intent(in) :: nors_t(nslip,3)
!     Schmid Dyadic
      real(8), intent(in) :: Schmid_0(nslip,3,3)
!
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
!     hardening model no.
      integer, intent(in) :: hardeningmodel
!     hardening model parameters
      real(8), intent(in) :: hardeningparam(maxnparam)
!     backstress model parameters
      real(8), intent(in) :: backstressparam(maxnparam)
!
!     Interaction matrices
!     Strength interaction between dislocations
      real(8), intent(in) :: sintmat1(nslip,nslip)
!     Strength interaction dislocation loops related with irradiation
      real(8), intent(in) :: sintmat2(nslip,nslip)
!     Latent hardening
      real(8), intent(in) :: hintmat1(nslip,nslip)
!     Hardening interaction matrix between dislocations
      real(8), intent(in) :: hintmat2(nslip,nslip)
!
!
!     overall crss
      real(8), intent(in) :: tauceff_t(nslip)
!     crss at former time step
      real(8), intent(in) :: tauc_t(nslip)
!     total dislocation density over all slip systems at the former time step
      real(8), intent(in) :: rhotot_t(nslip)
!     total scalar dislocation density over all slip systems at the former time step
      real(8), intent(in) :: sumrhotot_t
!     total dislocation density over all slip systems at the former time step
      real(8), intent(in) :: ssdtot_t
!     total forest dislocation density per slip system at the former time step
      real(8), intent(in) :: rhofor_t(nslip)
!     forest dislocation density per slip system at the former time step (hardening model = 4)
      real(8), intent(in) :: forest_t(nslip)
!     substructure dislocation density at the former time step
      real(8), intent(in) :: substructure_t
!     statistically-stored dislocation density per slip system at the former time step
      real(8), intent(in) :: gnd_t(nslip+nscrew)
!     statistically-stored dislocation density per slip system at the former time step
      real(8), intent(in) :: ssd_t(nslip)
!     defect loop density per slip system at the former time step
      real(8), intent(in) :: loop_t(maxnloop)
!     backstress at former time step
      real(8), intent(in) :: X_t(nslip)
!     time increment
      real(8), intent(in) :: dt
!     mechanical strain increment
      real(8), intent(in) :: dstran(6)
!     mechanical spin increment
      real(8), intent(in) :: domega(3)
!
!
!     OUTPUTS
!
!     plastic part of the deformation gradient
      real(8), intent(out) :: Fp(3,3)
!     Crystal to sample transformation martrix at current time step
      real(8), intent(out) :: gmatinv(3,3)
!     Green-Lagrange strains in the crystal reference
      real(8), intent(out) :: Eec(6)
!     slip rates at the current time step
      real(8), intent(out) :: gammadot(nslip)
!     total slip per slip system accumulated over the time
!     at the current time step
      real(8), intent(out) :: gammasum(nslip)
!     overall total slip at the current time step
      real(8), intent(out) :: totgammasum
!     Von-Mises equivalent total plastic strain at the current time step
      real(8), intent(out) :: evmp
!     crss at the current time step
      real(8), intent(out) :: tauc(nslip)
!     solute strength due to irradiation hardening
      real(8), intent(out) :: tausolute
!     total dislocation density over all slip systems at the current time step
      real(8), intent(out) :: ssdtot
!     forest dislocation density per slip system at the current time step
      real(8), intent(out) :: forest(nslip)
!     substructure dislocation density at the current time step
      real(8), intent(out) :: substructure
!     statistically-stored dislocation density per slip system at the current time step
      real(8), intent(out) :: ssd(nslip)
!     defect loop density per slip system at the current time step
      real(8), intent(out) :: loop(maxnloop)
!     backstress at the current time step
      real(8), intent(out) :: X(nslip)
!     Cauchy stress
      real(8), intent(out) :: sigma(6)
!     material tangent
      real(8), intent(out) :: jacobi(6,6)
!     convergence flag
      integer, intent(out) :: cpconv
!
!
!     Variables used within
!     deviatoric strains
      real(8) :: dev
!     Plastic spin dyadic (without shears)
      real(8) :: W(3,nslip), W33(3,3)
!     Resolved shear stress
      real(8) :: tau_t(nslip)
!     Two coefficients used in the solution
!     ddemsd = D * P + beta
      real(8) :: ddemsd(6,nslip)
!     Stress correction for large rotations
!     beta = sigma * W - W * sigma
      real(8) :: beta(6,nslip)
!     Results from the constitutive laws
      real(8) :: Lp_s(3,3), Lp_c(3,3)
      real(8) :: Dp_s(3,3), Dp_c(3,3)
!     Quad-precision variables
      real(16) :: Pmat_s(6,6), Pmat_c(6,6)
      real(16) :: gammadot_s(nslip), gammadot_c(nslip)    
!     derivative of slip rates wrto rss for slip
      real(8) :: dgammadot_dtau_s(nslip)
!     derivative of slip rates wrto rss for creep
      real(8) :: dgammadot_dtau_c(nslip)
!     derivative of slip rates wrto crss for slip
      real(8) :: dgammadot_dtauc_s(nslip)
!     derivative of slip rates wrto crss for creep
      real(8) :: dgammadot_dtauc_c(nslip)
!     total derivative of slip rates wrto rss
      real(8) :: dgammadot_dtau(nslip)
!     total derivative of slip rates wrto crss
      real(8) :: dgammadot_dtauc(nslip)
!     Plastic strain increment related quantities
      real(8) :: Lp(3,3), pdot, Dp(3,3)
!     Hardening increment mapping (numerically calculated)
      real(8) :: Hab(nslip,nslip)
!
!
!
!
!
!
!     Variables used in the solution for slip increments
      real(8) :: Nab(nslip,nslip)
      real(8) :: Mab(nslip,nslip)
!     The derivative of shear rates with respecto to rss
      real(8) :: ddgdde(nslip,6)
!     Slip increments
      real(8) :: dgamma(nslip), dgamma0(nslip)
!     RSS increments
      real(8) :: dtau(nslip), dtau0(nslip)
!     Total spin increment
      real(8) :: domega33(3,3)
!     Plastic spin increment
      real(8) :: domega33_p(3,3)
!     Elastic spin increment
      real(8) :: domega33_e(3,3)
!     Rotation matrix increment
      real(8) :: dgmatinv(3,3)
!
!     Increment in Cauchy stress
      real(8) :: dsigma(6), dsigma0(6)
!
!
!     Determinant of plastic deformation gradient
      real(8) :: detFp
!     Strain increment and related variables
      real(8) :: dstranp(6), dstrane(6)
      real(8) :: dstrane33(3,3), dstranp33(3,3)
!     rss increment
      real(8) :: dtau(nslip), dtau0(nslip)
!     crss increment
      real(8) :: dtauc(nslip), dtauc0(nslip)
!     ssd density increment
      real(8) :: dssd(nslip), dssd0(nslip)
!     loop density increment
      real(8) :: dloop(maxnloop), dloop0(maxnloop)
!     backstress increment
      real(8) :: dX(nslip), dX0(nslip)
!     total ssd density increment
      real(8) :: dssdtot, dssdtot0
!     forest dislocation density increment
      real(8) :: dforest(nslip), dforest0(nslip)
!     substructure dislocation density increment
      real(8) :: dsubstructure, dsubstructure0
!
!     deformed slip direction in sample reference frame
      real(8) :: dirs(nslip,3)
!     deformed slip plane normal in sample reference frame
      real(8) :: nors(nslip,3)
!
!     deformed slip direction increment at sample reference frame
      real(8) :: ddirs(nslip,3), ddirs0(nslip,3)
!     deformed slip plane normal increment at sample reference frame
      real(8) :: dnors(nslip,3), dnors0(nslip,3)
!
!     Schmid related tensors
      real(8) :: Schmid(nslip,3,3)
      real(8) :: Schmidvec(nslip,6)
      real(8) :: SchmidxSchmid(nslip,6,6)
      real(8) :: sdir(3), ndir(3), SNij(3,3), NSij(3,3)
!
!     overall crss
      real(8) :: tauceff(nslip)
!     forest density
      real(8) :: rhofor(nslip)
!     total density
      real(8) :: rhotot(nslip)
!     total scalar density
      real(8) :: sumrhotot
!     residual
      real(8) :: res(nslip)
!
!     Dummy variables
      real(8) :: dummy33(3,3), dummy33_(3,3)
      real(8) :: dummy0, dummy6(6), dummy66(6,6)
!
!     Variables uses within the subroutine
      integer :: is, js, il, i, j, k, a, b
      integer :: it
!
!
!
!     Initiate convergence flag
      cpconv = 1
!
!     Assign initial variables
!
!     Orientations
      gmatinv = gmatinv_t
!
!     Slip vectors
      dirs = dirs_t
      dirn = dirn_t
!      
!     Stress
      sigma = sigma_t
!
!
!     shear
      gamma = gamma_t
!
!     Set the old increments to zero
      dsigma0=0.; dgamma0=0.
      dtau0=0.;  dtauc0=0.
      dssd0=0.; dX0=0.
      dloop0=0.; dssdtot0=0.
      dforest0=0.; dsubstructure0=0.
      ddirs0=0.; ddirn0=0.
!
!     Assign former state variables
      tauc = tauc_t
      X = X_t
      ssd = ssd_t
      ssdtot = ssdtot_t
      forest = forest_t
      loop = loop_t
      substructure = substructure_t
!
!
!
!     Schmid tensor
      do is = 1, nslip
          Schmid_vec(is,1) = dirs(is,1)*dirn(is,1)
          Schmid_vec(is,2) = dirs(is,2)*dirn(is,2)
          Schmid_vec(is,3) = dirs(is,3)*dirn(is,3)
          Schmid_vec(is,4) = dirs(is,1)*dirn(is,2)+dirs(is,2)*dirn(is,1)
          Schmid_vec(is,5) = dirs(is,1)*dirn(is,3)+dirs(is,3)*dirn(is,1)
          Schmid_vec(is,6) = dirs(is,2)*dirn(is,3)+dirs(is,3)*dirn(is,2)
      end do
      
!     Calculate RSSS
      do is = 1, nslip
          tau_t(is) = dot_product(Schmidvec(is,:),sigma_t)
      end do
!
!
!
!     RSS
      tau = tau_t
      
      
!
!     Volumetric change in the strain
      dev = dstran(1) + dstran(2) + dstran(3)
!
!
      
      
      
!     Iterations start here   
      do it = 1, maxniter
      
!         Calculate Schmid tensors and Schmid dyadic
          Schmid=0.; SchmidxSchmid=0.
          Schmidvec=0.; 
          do is=1, nslip
!
!             Slip direction
              sdir = dirs_t(is,:)
!             Slip plane normal
              ndir = nors_t(is,:)
!
              do i=1,3
                  do j=1,3
                      SNij(i,j) = sdir(i)*ndir(j)
                      NSij(i,j) = ndir(j)*sdir(i)
                      Schmid(is,i,j) = SNij(i,j)
                  enddo
              enddo
!
!
!
!
              call gmatvec6(SNij,sni)
!
              call gmatvec6(NSij,nsi)
!
!             Vectorized Schmid tensor
              Schmidvec(is,1:6) = sni
!
              do i=1,6
                  do j=1,6
                      SchmidxSchmid(is,i,j)=sni(i)*nsi(j)
                  enddo
              enddo
!
!
              Schmid_vec(is,1) = dirs(is,1)*dirn(is,1)
              Schmid_vec(is,2) = dirs(is,2)*dirn(is,2)
              Schmid_vec(is,3) = dirs(is,3)*dirn(is,3)
              Schmid_vec(is,4) = dirs(is,1)*dirn(is,2)+
     + dirs(is,2)*dirn(is,1)
              Schmid_vec(is,5) = dirs(is,1)*dirn(is,3)+
     + dirs(is,3)*dirn(is,1)
              Schmid_vec(is,6) = dirs(is,2)*dirn(is,3)+
     + dirs(is,3)*dirn(is,2)
!
          end do          
          
      
!         Plastic spin dyadic
          W = 0.
          do is = 1, nslip
!
              W(1,is) = 0.5*(dirs_t(is,1)*nors_t(is,2)-
     + dirs_t(is,2)*nors_t(is,1))
!
              W(2,is) = 0.5*(dirs_t(is,3)*nors_t(is,1)-
     + dirs_t(is,1)*nors_t(is,3))
!
              W(3,is) = 0.5*(dirs_t(is,2)*nors_t(is,3)-
     + dirs_t(is,3)*nors_t(is,2))
!
          end do
!
!
!

!
!
!         Calculate beta and ddemsd
          beta=0.; ddemsd=0.
          do is = 1, nslip
!
!             Symbolic math result
              beta(1,is) = -2.*W(2,is)*sigma_t(5) + 
     + 2.*W(1,is)*sigma_t(4)
              beta(2,is) = -2.*W(1,is)*sigma_t(4) + 
     + 2.*W(3,is)*sigma_t(6)
              beta(3,is) = -2.*W(3,is)*sigma_t(6) +
     + 2.*W(2,is)*sigma_t(5)
              beta(4,is) = -W(1,is)*sigma_t(1) + W(1,is)*sigma_t(2) - 
     + W(2,is)*sigma_t(6) + W(3,is)*sigma_t(5)
          beta(5,is) = -W(2,is)*sigma_t(3) + W(2,is)*sigma_t(1) + 
     + W(1,is)*sigma_t(6) - W(3,is)*sigma_t(4)
          beta(6,is) = -W(3,is)*sigma_t(2) - W(1,is)*sigma_t(5) + 
     + W(2,is)*sigma_t(4) + W(3,is)*sigma_t(3)
!

              ddemsd(:,is) = matmul(Cs,Schmidvec(is,:)) +
     + beta(:,is)
!
          end do
!
!
!
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
              dgammadot_dtau_s = 0.
              dgammadot_dtauc_s = 0.
!
!         sinh law
          elseif (slipmodel == 1) then
!
              call sinhslip(Schmid_0,Schmid,SchmidxSchmid,
     + tau_t,X_t,tauceff_t,rhofor_t,burgerv,dt,
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
     + tau_t,X_t,tauceff_t,burgerv,dt,nslip,phaid,
     + mattemp,slipparam,irradiationmodel,
     + irradiationparam,cubicslip,caratio,
     + Lp_s,Dp_s,Pmat_s,gammadot_s,
     + dgammadot_dtau_s,dgammadot_dtauc_s) 
!
!
!         power law
          elseif (slipmodel == 3) then
!
!
              call powerslip(Schmid_0,Schmid,SchmidxSchmid,
     + tau_t,X_t,tauceff_t,burgerv,dt,
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
              dgammadot_dtau_c = 0.
              dgammadot_dtauc_c = 0.
!
          elseif (creepmodel == 1) then
!
!
!
              call expcreep(Schmid_0,Schmid,SchmidxSchmid,
     + tau_t,X_t,tauceff_t,dt,nslip,phaid,
     + mattemp,creepparam,gammasum_t,
     + Lp_c,Dp_c,Pmat_c,gammadot_c,
     + dgammadot_dtau_c,dgammadot_dtauc_c)
!
!
!
!
          endif
!
!
!         Sum the effects of creep and slip rates
          gammadot = gammadot_s + gammadot_c
!
          dgammadot_dtau = dgammadot_dtau_s +
     + dgammadot_dtau_c
!
          dgammadot_dtauc = dgammadot_dtauc_s +
     + dgammadot_dtauc_c
!
          Lp = Lp_s + Lp_c
          Dp = Dp_s + Dp_c
!
!
!         Check for the slip rates
          if(any(gammadot /= gammadot)) then
!             did not converge
              cpconv = 0
!             enter dummy stress and jacobian
              sigma = 0.
              jacobi = I6
!             return to end of the subroutine
!             warning message
              call error(14)
              return
          endif      
!
!
!         Plastic strain-related quantities for hardening calculations
!
!
!
!         calculate von mises invariant plastic strain rate
          pdot=sqrt(2./3.*sum(Dp*Dp))
!
!
!
!         Total slip over time per slip system
          gammasum = 0.
          do is =1, nslip
!
              gammasum(is) = gammasum_t(is) +
     + abs(gammadot(is))*dt
!
          enddo
!
!
!         Total slip
          totgammasum = totgammasum_t +
     + sum(abs(gammadot))*dt
!
!
!
!
!         Update backstress
          if (backstressmodel==1) then
!
              call backstressmodel1(backstressparam,
     + nslip,X_t,gammadot,dt,dX)
!
          end if
!

!
!
!         Update the states using hardening laws
          call hardeningrules(phaid,nslip,
     + mattemp,dt,G12,burgerv,
     + totgammasum,gammadot,pdot,
     + irradiationmodel,irradiationparam,
     + hardeningmodel,hardeningparam,
     + hintmat1,hintmat2,
     + tauc_t,ssd_t,loop_t,
     + rhofor_t,substructure_t,
     + tausolute,dtauc,dssdtot,dforest,
     + dsubstructure,dssd,dloop)
!
!
!
!
!     Update the hardening states
!
      tauc = tauc_t + dtauc
!
      ssd = ssd_t + dssd
!
      loop = loop_t + dloop
!
      ssdtot = ssdtot_t + dssdtot
!
      forest = forest_t + dforest
!
      substructure = substructure_t + dsubstructure
!
!
!      
!     Check if the statevariables going negative due to softening
!     This may happen at high temperature and strain rates constants going bad
      if(any(tauc < 0.)) then
!         did not converge
          cpconv = 0
!         enter dummy stress and jacobian
          sigma = 0.
          jacobi = I6
!         return to end of the subroutine
!         warning message
          call error(21)
          return
      endif
!
      if(any(ssd < 0.)) then
!         did not converge
          cpconv = 0
!         enter dummy stress and jacobian
          sigma = 0.
          jacobi = I6
!         return to end of the subroutine
!         warning message
          call error(21)
          return
      endif
!
!
!
      if(any(forest < 0.)) then
!         did not converge
          cpconv = 0
!         enter dummy stress and jacobian
          sigma = 0.
          jacobi = I6
!         return to end of the subroutine
!         warning message
          call error(21)
          return
      endif
!
      if(substructure < 0.) then
!         did not converge
          cpconv = 0
!         enter dummy stress and jacobian
          sigma = 0.
          jacobi = I6
!         return to end of the subroutine
!         warning message
          call error(21)
          return
      endif
!
!
!     Find the effective hardening increment
!
!     Calculate total and forest density
      call totalandforest(phaid,
     + nscrew, nslip, gnd_t,
     + ssd, ssdtot, forest,
     + forestproj, slip2screw, rhotot,
     + sumrhotot, rhofor)
!
!
!
!     Calculate crss
      call slipresistance(phaid, nslip, gf, G12,
     + burgerv, sintmat1, sintmat2,
     + tauc, rhotot, sumrhotot, rhofor,
     + substructure, tausolute, loop,
     + hardeningmodel, hardeningparam,
     + irradiationmodel, irradiationparam,
     + mattemp, tauceff)
!
!
!
!!     Numerical calculation of derivative dtauc/dgamma
!      Hab = 0.
!      do is = 1, nslip
!!
!          if (abs(gammadot(is))>sqrt(smallnum)) then
!!
!              Hab(is,is) = (tauceff(is)-tauceff_t(is))/dt/
!     + abs(gammadot(is))
!!
!!
!          end if
!!
!      end do
!     
!
!
!
!
!
!     Euler solution
      Nab = 0.
      dgamma = 0.
      do a = 1, nslip
!
          do b = 1, nslip
!
              dummy0 = 0.
              do i = 1, 6
                  dummy0 = dummy0 + ddemsd(i,a)*Schmidvec(b,i)
              end do
!
              Nab(a,b) = theta*dt*dgammadot_dtau(a)*dummy0 -
     + hintmat1(a,b)*dgammadot_dtauc(a)*sign(1.0,gammadot_t(b))*theta*dt
!
!
          end do
!
          Nab(a,a) = Nab(a,a) + 1.
!
!         Given quantities in vector form
          dgamma(a) = dot_product(ddemsd(:,a),dstran)*
     + dgammadot_dtau(a)*theta*dt + gammadot_t(a)*dt
!
!
      end do
!
!
!     Solve for the slip increments
      call nolapinverse(Nab,Mab,nslip)
!
!
!
!
!     Check for the inversion
      if(any(Mab /= Mab)) then
!         did not converge
          cpconv = 0
!         enter dummy stress and jacobian
          sigma = 0.
          jacobi = I6
!         return to end of the subroutine
!         warning message
          call error(14)
          return
      endif
!
!
!
!
!     Slip increments
      dgamma = matmul(Mab,dgamma)
!
!
!!     residual
!      res = dgamma -(1-theta)*dt*gammadot_t -theta*dt*gammadot
!!
!!
!      if (maxval(abs(res))>tolerance) then
!          cpconv=0
!      end if
!
!
!     Slip rates
      gammadot = dgamma/dt
!
!
!
!
!
!
!     Redo the slip calculations
!
!
!     Plastic velocity gradient
      Lp = 0.
      do is = 1, nslip
          Lp = Lp + gammadot(is)*Schmid_0(is,:,:)
      end do
!
!
!
!
!     Total slip over time per slip system
      gammasum = 0.
      do is = 1, nslip
!
          gammasum(is) = gammasum_t(is) +
     + abs(gammadot(is))*dt
!
      enddo
!
!
!     Total slip
      totgammasum = totgammasum_t +
     + sum(abs(gammadot))*dt     
!
!
!
!     variables for plastic part of the deformation gradient
      dummy33 = I3 - Lp*dt
      call inv3x3(dummy33,dummy33_,dummy0)
!
!     plastic part of the deformation gradient
      Fp = matmul(dummy33_,Fp_t)
!
!     determinant
      call deter3x3(Fp,detFp)
!
!
!
!
!     check wheter the determinant is negative
!     or close zero
      if (detFp <= smallnum) then
!         did not converge
          cpconv = 0
!         enter dummy stress and jacobian
          sigma = 0.
          jacobi = I6
!         return to end of the subroutine
!         warning message
          call error(19)
          return
      else
!         Scale Fp with its determinant to make it isochoric
          Fp = Fp / detFp**(1./3.)
!
      end if     
!
!
!     calculate von mises invariant plastic strain rate
      pdot=sqrt(2./3.*sum(Dp*Dp))
!
!
!     Von-Mises equivalent total plastic strain
      evmp = evmp_t + pdot*dt
!
!
!
!
!
!
!
!
!
!     Spin increment (3x3 matrix)
      domega33 = 0.
      domega33(1,2) = domega(1)
      domega33(1,3) = -domega(2)
      domega33(2,1) = -domega(1)
      domega33(2,3) = domega(3)
      domega33(3,1) = domega(2)
      domega33(3,2) = -domega(3)
!
!
!     Calculate plastic spin from slip
      domega33_p=0.
      do is = 1, nslip
!
          W33 = 0.
          W33(1,2) = W(1,is)
          W33(1,3) = -W(2,is)
          W33(2,1) = -W(1,is)
          W33(2,3) = W(3,is)
          W33(3,1) = W(2,is)
          W33(3,2) = -W(3,is)        
!
          domega33_p = domega33_p + W33*dgamma(is)
!
!
      end do
!
!
!     Elastic spin
      domega33_e = domega33 - domega33_p
!
!
!
!
!
!     Orientation update
      dgmatinv = matmul(domega33_e,gmatinv_t)
!
!
      gmatinv = gmatinv_t + dgmatinv
!
!
!
!
!     Elastic strains in the crystal reference
!
!
!     Elastic strain increment
      dstranp=0.
      do is = 1, nslip
          dstranp(:) = dstranp(:) + Schmidvec(is,:)*dgamma(is)
      end do
!
!
!     Subtract the plastic strain increment from total
      dstrane = dstran - dstranp
!
!     undo shear corrections
      dstrane(4:6) = 0.5*dstrane(4:6)
!
!
!     Convert to 3x3 matrix
      call vecmat6(dstrane,dstrane33)
!
!
!     Elastic strains in the crystal reference
      dummy33_ = matmul(transpose(gmatinv),dstrane33)
      dummy33 = matmul(dummy33_,gmatinv)
!
!     Vectorize
      call matvec6(dummy33,dummy6)
!
!     Shear corrections
      dummy6(4:6) = 2.0*dummy6(4:6)
!
!     Add the strain increment to the former value
      Eec=Eec_t + dummy6
!
!
!
!
!
!
!
!
!     Stress update
      dsigma = matmul(Cs,dstran) -
     + sigma_t*dev - matmul(ddemsd,dgamma)
!
      sigma = sigma_t + dsigma
!
!
!     End of Newton-Raphson loop
      end do
!
!
!     Material tangent
!
!     Step-1: Calculate ddgamma_ddeps
      do is=1,nslip
!
          ddgdde(is,:) = dt*theta*dgammadot_dtau(is)*ddemsd(:,is)
!
      end do
!
!
!
!
!     Step-2: Calculate overall expression
      jacobi = Cs - matmul(ddemsd,matmul(Mab,ddgdde))
!
!
!     Correction for large deformations
      do i=1,3
          do j=1,3
              jacobi(i,j) = jacobi(i,j) - sigma(i)
              jacobi(i+3,j) = jacobi(i+3,j) - sigma(i+3)
          end do
      end do
!
!
      jacobi = jacobi / (1.+dev)
!
!     Make it symmetric to help convergence
      jacobi = 0.5*(jacobi + transpose(jacobi))
!
!
!
      return
      end subroutine CP_Huang_fullyimplicit
!