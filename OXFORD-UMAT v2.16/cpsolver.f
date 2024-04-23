!     Oct. 1st, 2022
!     Eralp Demir
!     This module contains the solver schemes for crystal plasticity
!
      module cpsolver
      implicit none
      contains
!
!     This subroutine deals with
!     global variables and their assignments
!     before entering the main solver:
!     1. assigns the global variables to locals
!     before calling the CP-solver
!     2. calls fge-predictor scheme before as
!     spare solution
!     3. calls implicit or explicit CP-solvers
!     4. assigns the results to the glboal state variables
      subroutine solve(noel, npt, dfgrd1, dfgrd0,
     + temp, dtemp, dt, matid, pnewdt, nstatv, statev,
     + sigma, jacobi)
!
      use globalvariables, only : dt_t,
     + statev_gmatinv, statev_gmatinv_t,
     + statev_gammasum_t, statev_gammasum, statev_ssdtot_t,
     + statev_gammadot_t, statev_gammadot, statev_ssdtot,
     + statev_Fp_t, statev_Fp, statev_sigma_t, statev_sigma,
     + statev_jacobi_t, statev_jacobi, statev_Fth_t, statev_Fth,
     + statev_tauc_t, statev_tauc, statev_maxx_t, statev_maxx,
     + statev_Eec_t, statev_Eec, statev_gnd_t, statev_gnd,
     + statev_ssd_t, statev_ssd, statev_loop_t, statev_loop,
     + statev_forest_t, statev_forest, statev_evmp_t, statev_evmp,
     + statev_substructure_t, statev_substructure, statev_sigma_t2,
     + statev_totgammasum_t, statev_totgammasum, 
     + statev_tausolute_t, statev_tausolute, forestproj_all,
     + numslip_all, numscrew_all, phaseid_all, Schmid_0_all,
     + dirc_0_all, norc_0_all, caratio_all, cubicslip_all,
     + Cc_all, gf_all, G12_all, alphamat_all, burgerv_all,
     + sintmat1_all, sintmat2_all, hintmat1_all, hintmat2_all,
     + slipmodel_all, slipparam_all, creepmodel_all, creepparam_all,
     + hardeningmodel_all, hardeningparam_all, irradiationmodel_all,
     + irradiationparam_all, backstressparam_all, slip2screw_all,
     + statev_backstress_t, statev_backstress, statev_plasdiss_t,
     + statev_plasdiss, statev_tauceff,
     + I3, I6, smallnum
!
      use userinputs, only: constanttemperature, temperature,
     + predictor, maxnslip, maxnparam, maxxcr, cutback,
     + phi, maxnloop, stateupdate
!
!
      use usermaterials, only: materialparam
!
      use crss, only: slipresistance, totalandforest
!
      use useroutputs, only: assignoutputs, nstatv_outputs
!
      use errors, only: error
!
      use utilities, only: inv3x3, nolapinverse, matvec6, vecmat6,
     + rotord4sig, gmatvec6
!
      implicit none
!
!     element no
      integer, intent(in) :: noel
!
!     ip no
      integer, intent(in) :: npt
!
!
!     current deformation gradient
      real(8), intent(in) :: dfgrd1(3,3)
!
!     former deformation gradient
      real(8), intent(in) :: dfgrd0(3,3)
!
!     ABAQUS temperature
      real(8), intent(in) :: temp
!
!     ABAQUS temperature increment
      real(8), intent(in) :: dtemp
!
!     time step
      real(8), intent(in) :: dt
!
!     material id
      integer, intent(in) :: matid
!
!     time factor
      real(8), intent(inout) :: pnewdt
!
!     number of state variables - for postprocessing
      integer, intent(in) :: nstatv
!
!     state variables - postprocessing
      real(8), intent(inout) :: statev(nstatv)
!
!     Cauchy stress
      real(8), intent(out) :: sigma(6)
!
!     Material tangent
      real(8), intent(out) :: jacobi(6,6)
!
!     Local variables used within this subroutine
!
!
!     phase-id
      integer :: phaid
!
!     number of slip systems
      integer :: nslip
!
!
!     Convergence flag (initially set to zero!)
!
!     Flag for crystal plasticity explicit/implicit solver
      integer :: cpconv
!
!     Flag for Euler solver convergence
      integer :: cpconv0
!
!
!     Local state variables with known dimensions
!     crystal to sample transformation at the former time step
      real(8) :: gmatinv_t(3,3)
!     crystal to sample transformation at the former time step
      real(8) :: gmatinv(3,3)
!     stress at the former time step
      real(8) :: sigma_t(6)
!     rss/crsss ratio at the former time step
      real(8) :: maxx_t
!     rss/crsss ratio at the current time step 
      real(8) :: maxx
!     elastic strains in the crystal reference at the former time step
      real(8) :: Eec_t(6)
!     elastic strains in the crystal reference at the current time step
      real(8) :: Eec(6)
!     plastic deformation gradient at the former time step
      real(8) :: Fp_t(3,3)
!     plastic deformation gradient at the current time step
      real(8) :: Fp(3,3)
!     thermal deformation gradient at the former time step
      real(8) :: Fth_t(3,3)
!     thermal deformation gradient at the current time step
      real(8) :: Fth(3,3)
!     inverse of the thermal deformation gradient at the former time step
      real(8) :: invFth_t(3,3)
!     inverse of the thermal deformation gradient at the current time step
      real(8) :: invFth(3,3)
!     determinant of the thermal deformation gradient
      real(8) :: detFth, detFth_t
!     mechanical deformation gradient at the former time step
      real(8) :: F_t(3,3)
!     mechanical deformation gradient at the current time step
      real(8) :: F(3,3)
!
!     Variables for velocity gradient calculation
!     Fdot
      real(8) :: Fdot(3,3)
!     velocity gradient at the current time step
      real(8) :: L(3,3)
!     inverse of the deformation gradient
      real(8) :: Finv(3,3)
!     determinant of the deformation gradient
      real(8) :: detF
!
!     Local state variable arrays
!     sum of slip per slip system
      real(8) :: gammasum_t(numslip_all(matid))
      real(8) :: gammasum(numslip_all(matid))
!     slip rates
      real(8) :: gammadot_t(numslip_all(matid))
      real(8) :: gammadot(numslip_all(matid))
!     slip resistance
      real(8) :: tauc_t(numslip_all(matid))
      real(8) :: tauc(numslip_all(matid))
!     effective overall slip resistance
!     (tauc0 + GND + SSD + solute + substructure + forest + etc.)
      real(8) :: tauceff_t(numslip_all(matid))
      real(8) :: tauceff(numslip_all(matid))
!     Note the size of GND is different
!     gnd density (nslip + nscrew)
      real(8) :: gnd_t(numslip_all(matid)+numscrew_all(matid))
      real(8) :: gnd(numslip_all(matid)+numscrew_all(matid))
!
!     ssd density (nslip)
      real(8) :: ssd_t(numslip_all(matid))
      real(8) :: ssd(numslip_all(matid))
!     loop density (maxnloop)
      real(8) :: loop_t(maxnloop)
      real(8) :: loop(maxnloop)
!     total forest dislocation density - derived from other terms
      real(8) :: rhofor_t(numslip_all(matid))
      real(8) :: rhofor(numslip_all(matid))
!     total density
      real(8) :: rhotot_t(numslip_all(matid))
      real(8) :: rhotot(numslip_all(matid))
!     forest dislocation density as a state variable
      real(8) :: forest_t(numslip_all(matid))
      real(8) :: forest(numslip_all(matid))
!
!
!     Scalar state variables
!     equivalent Von-Mises plastic strain
      real(8) :: evmp_t
      real(8) :: evmp
!
!     cumulative slip
      real(8) :: totgammasum_t
      real(8) :: totgammasum
!
!     plastic dissipation power
      real(8) :: plasdiss_t
      real(8) :: plasdiss
!
!     solute strength
      real(8) :: tausolute_t
      real(8) :: tausolute
!     substructure density
      real(8) :: substructure_t
      real(8) :: substructure
!     total density
      real(8) :: ssdtot_t
      real(8) :: ssdtot
!     scalar cumulartive density
      real(8) :: sumrhotot_t
      real(8) :: sumrhotot     
!
!     material-related local variables
!     number screw systems
      integer :: nscrew
!     slip model flag
      integer :: smodel
!     creep model flag
      integer :: cmodel
!     hardening model flag   
      integer :: hmodel
!     irradiation mdoel flag
      integer :: imodel
!     cubic slip flag
      integer :: cubicslip
!     material temperature
      real(8) :: mattemp
!     c/a ratio for hcp materials
      real(8) :: caratio
!     compliance at the crystal reference
      real(8) :: Cc(6,6)
!     geometric factor
      real(8) :: gf
!     shear modulus
      real(8) :: G12
!     Poisson's ratio
      real(8) :: v12
!     thermal expansion coefficient
      real(8) :: alphamat(3,3)
!     slip parameters
      real(8) :: sparam(maxnparam)
!     creep parameters
      real(8) :: cparam(maxnparam)
!     hardening parameters
      real(8) :: hparam(maxnparam)
!     irradiation parameters
      real(8) :: iparam(maxnparam)
!     Backstress parameters
      real(8) :: bparam(maxnparam)
!     Burgers vector
      real(8) :: burgerv(numslip_all(matid))
!     Interaction matrices
!     Strength interaction between dislocations
      real(8) :: sintmat1(numslip_all(matid),numslip_all(matid))
!     Strength interaction dislocation loops related with irradiation
      real(8) :: sintmat2(numslip_all(matid),numslip_all(matid))
!     Latent hardening
      real(8) :: hintmat1(numslip_all(matid),numslip_all(matid))
!     Hardening interaction matrix between dislocations
      real(8) :: hintmat2(numslip_all(matid),numslip_all(matid))
!
!
!     slip direction and slip plane normal at the crystal reference (undeformed)
      real(8) :: dirc_0(numslip_all(matid),3)
      real(8) :: norc_0(numslip_all(matid),3)
!     slip direction and slip plane normal at the sample reference
      real(8) :: dirs_t(numslip_all(matid),3)
      real(8) :: nors_t(numslip_all(matid),3)
!
!     Forest projection operators for GND and SSD
      real(8) :: forestproj(numslip_all(matid),
     + numslip_all(matid)+numscrew_all(matid))
!
!     Slip to screw system map
      real(8) :: slip2screw(numscrew_all(matid),numslip_all(matid))
!
!     Schmid Dyadic
      real(8) :: Schmid_0(numslip_all(matid),3,3)
      real(8) :: Schmid(numslip_all(matid),3,3)
      real(8) :: Schmidvec(numslip_all(matid),6)
      real(8) :: SchmidxSchmid(numslip_all(matid),6,6)
      real(8) :: sdir(3), ndir(3), SNij(3,3), NSij(3,3)
      real(8) :: nsi(6), sni(6)
!
!
!
!     strain calculations
!     total strain increment
      real(8) ::  dstran(6), dstran33(3,3)
!     thermal strain increment
      real(8) ::  dstranth33(3,3), dstranth(6)
!     total spin increment
      real(8) ::  domega(3)
!     total spin (rate) 3x3 matrix
      real(8) ::  W(3,3), dW33(3,3)
!
!     Elasticity transformation
!     temporary array for elastic stiffness calculation
      real(8) :: rot4(6,6)
!     elasticity matrix at the deformed reference
      real(8) :: Cs(6,6)     
!
!     Trial stress calculation
!     trial stress
      real(8) :: sigmatr(6)
!
!     Backstress at current time
      real(8) :: X(numslip_all(matid))
!     Backstress at former time
      real(8) :: X_t(numslip_all(matid))
!     Backstress backup solution
      real(8) :: X0(numslip_all(matid))
!
!     stress with rotations
      real(8) :: sigmarot_t(6)
!
!     stress (3x3)
      real(8) ::  sigma33_t(3,3)
!
!
!     Resolved shear stress calculation
!     GUESS values
      real(8) :: sigma0(6), jacobi0(6,6)
!     value of resolved shear stress
      real(8) :: tau0(numslip_all(matid))
!
!     value of resolved shear stress
      real(8) :: tau_t(numslip_all(matid))
!
!     value of trial resolved shear stress
      real(8) :: tautr(numslip_all(matid))
!     absolute value of resolved shear stress
      real(8) :: abstautr(numslip_all(matid))
!
!     dummy variables
      real(8) :: dummy3(3), dummy33(3,3)
      real(8) :: dummy33_(3,3)
      real(8) :: dummy6(6), dummy66(6,6)
      integer :: dum1, dum2, dum3
!     unused variables
      real(8) :: notused1(maxnslip)
      real(8) :: notused2(maxnslip)
      real(8) :: notused3(maxnslip)
!     In case materials subroutine entered everytime
!     Strength interaction between dislocations
      real(8) :: notused4(maxnslip,maxnslip)
!     Strength interaction dislocation loops related with irradiation
      real(8) :: notused5(maxnslip,maxnslip)
!     Latent hardening
      real(8) :: notused6(maxnslip,maxnslip)
!     Hardening interaction matrix between dislocations
      real(8) :: notused7(maxnslip,maxnslip)
!     Initial forest and substructure densities
      real(8) :: notused8, notused9
!
!     counter
      integer :: is, i, j
!
!
!
!     Reset sparse arrays
      notused1=0.;notused2=0.;notused3=0.
      notused4=0.;notused5=0.
      notused6=0.;notused7=0.
!
!
!
!     convergence check initially
!     if sinh( ) in the slip law has probably blown up
!     then try again with smaller dt
      if(any(dfgrd1 /= dfgrd1)) then
!         Set the outputs to zero initially
          sigma = 0.
          jacobi = I6
!         cut back time
          pnewdt = cutback
!         warning message in .dat file
          call error(11)
!         go to the end of subroutine
          return
      end if
!     
!
!
!     phase-id
      phaid = phaseid_all(matid)
!
!     Number of slip systems
      nslip = numslip_all(matid)
!
!     Number of screw systems
      nscrew = numscrew_all(matid)
!
!
!
!
!     undeformed slip direction
      dirc_0 = dirc_0_all(matid,1:nslip,1:3)
      
!     undeformed slip plane normal
      norc_0 = norc_0_all(matid,1:nslip,1:3)      
!
!
!     Forest projection for GND
      forestproj = forestproj_all(matid,1:nslip,1:nslip+nscrew)
!
!
!     Slip to screw system mapping
      slip2screw = slip2screw_all(matid,1:nscrew,1:nslip)
!
!     Initial Schmid tensor
      Schmid_0 = Schmid_0_all(matid,1:nslip,:,:)
!
!     Assign the global state variables
!     to the local variables
!     at the former time step
      gammasum_t = statev_gammasum_t(noel,npt,1:nslip)
      gammadot_t = statev_gammadot_t(noel,npt,1:nslip)
      tauc_t = statev_tauc_t(noel,npt,1:nslip)
      gnd_t = statev_gnd_t(noel,npt,1:nslip+nscrew)
      ssd_t = statev_ssd_t(noel,npt,1:nslip)
      loop_t = statev_loop_t(noel,npt,1:maxnloop)
      ssdtot_t = statev_ssdtot_t(noel,npt)
      forest_t = statev_forest_t(noel,npt,1:nslip)
      substructure_t = statev_substructure_t(noel,npt)
      evmp_t = statev_evmp_t(noel,npt)
      plasdiss_t = statev_plasdiss_t(noel,npt)
      totgammasum_t = statev_totgammasum_t(noel,npt)
      tausolute_t = statev_tausolute_t(noel,npt)
      sigma_t = statev_sigma_t(noel,npt,:)
      Fp_t = statev_Fp_t(noel,npt,:,:)
      Fth_t = statev_Fth_t(noel,npt,:,:)
      Eec_t = statev_Eec_t(noel,npt,:)
!     Crystal orientations at former time step
      gmatinv_t = statev_gmatinv_t(noel,npt,:,:)
!     Backstress at former time step
      X_t = statev_backstress_t(noel,npt,1:nslip)
!
!
!     Material parameters are constant
      caratio = caratio_all(matid)
      cubicslip = cubicslip_all(matid)
      Cc = Cc_all(matid,:,:)
      gf = gf_all(matid)
      G12 = G12_all(matid)
      alphamat = alphamat_all(matid,:,:)
      burgerv = burgerv_all(matid,1:nslip)
!
!
      sintmat1 = sintmat1_all(matid,1:nslip,1:nslip)
      sintmat2 = sintmat2_all(matid,1:nslip,1:nslip)
      hintmat1 = hintmat1_all(matid,1:nslip,1:nslip)
      hintmat2 = hintmat2_all(matid,1:nslip,1:nslip)
!
!
!
      smodel = slipmodel_all(matid)
      sparam = slipparam_all(matid,1:maxnparam)
      cmodel = creepmodel_all(matid)
      cparam = creepparam_all(matid,1:maxnparam)
      hmodel = hardeningmodel_all(matid)
      hparam = hardeningparam_all(matid,1:maxnparam)
      imodel = irradiationmodel_all(matid)
      iparam = irradiationparam_all(matid,1:maxnparam)
      bparam = backstressparam_all(matid,1:maxnparam)
!
!
!
!
!     Temperature is constant and defined by the user
      if (constanttemperature == 1) then
!
!
!         Assign temperature
          mattemp = temperature
!
!
!
!     Temperature is defined by ABAQUS
!     Material properties are entered every time
!     Because properties can be temperature dependent
      else if (constanttemperature == 0) then
!
!         Use ABAQUS temperature (must be in K)
          mattemp = temp
!
!
!         get material constants
          call materialparam(matid,mattemp,
     + dum1,dum2,dum3,caratio,cubicslip,Cc,
     + gf,G12,v12,alphamat,notused1, ! state variables are NOT updated here
     + notused2,notused3,notused8,notused9,
     + smodel,sparam,cmodel,cparam,
     + hmodel,hparam,imodel,iparam,
     + notused4,notused5,notused6,
     + notused7,bparam) ! Interaction matrices are not updated here
!
!
!  
!    
!
      end if
!
!
!
!
!
!
!
!     Slip directions in the sample reference
      call rotateslipsystems(phaid,nslip,caratio,
     + gmatinv_t,dirc_0,norc_0,dirs_t,nors_t)
!
!
!     Calculate Schmid tensors and Schmid dyadic
      Schmid=0.; SchmidxSchmid=0.
      do is=1,nslip
!
!         Slip direction
          sdir = dirs_t(is,:)
!         Slip plane normal
          ndir = nors_t(is,:)
!
          do i=1,3
              do j=1,3
                  SNij(i,j) = sdir(i)*ndir(j)
                  NSij(j,i) = ndir(j)*sdir(i)
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
!         Vectorized Schmid tensor
          Schmidvec(is,1:6) = sni
!
          do i=1,6
              do j=1,6
                  SchmidxSchmid(is,i,j)=sni(i)*nsi(j)
              enddo
          enddo
!
      enddo
!
!
!
!     Calculate total and forest density
      call totalandforest(phaid, nscrew, nslip,
     + gnd_t, ssd_t,
     + ssdtot_t, forest_t,
     + forestproj, slip2screw, rhotot_t,
     + sumrhotot_t, rhofor_t)
!
!
!
!     Calculate crss
      call slipresistance(phaid, nslip, gf, G12,
     + burgerv, sintmat1, sintmat2, tauc_t,
     + rhotot_t, sumrhotot_t, rhofor_t, substructure_t,
     + tausolute_t, loop_t, hmodel, hparam, imodel, iparam,
     + mattemp, tauceff_t)
!
!
!
!
!     Elastic stiffness in the sample reference
!
!     Rotation matrix - special for symmetric 4th rank transformation
      call rotord4sig(gmatinv_t,rot4)
!
!
!
!     Elasticity tensor in sample reference
      dummy66=matmul(rot4,Cc)
      Cs = matmul(dummy66,transpose(rot4))
!
!!     To avoid numerical problems
!      Cs = (Cs + transpose(Cs))/2.
!
!
!
!
!
!     CALCULATION OF THERMAL STRAINS
!
!     No thermal strains
      if (constanttemperature == 1) then
!
!
          dstranth = 0.
!
          F = dfgrd1
!
          F_t = dfgrd0
!
          Fth = Fth_t
!
!     Thermal strains if temperature change is defined by ABAQUS
      else
!
!         Thermal eigenstrain in the crystal reference system
          dstranth33 = dtemp*alphamat
!
!         Transform the thermal strains to sample reference
          dstranth33 = matmul(matmul(gmatinv_t,dstranth33),
     + transpose(gmatinv_t))
!
!
!
!
!
!         Convert to a vector
          call matvec6(dstranth33,dstranth)
!
!         Shear corrections
          dstranth(4:6) = 2.0*dstranth(4:6)
!
!
!
!         Thermal deformation gradient
          Fth = Fth_t + dstranth33 
!
!         Invert the thermal distortions
          call inv3x3(Fth,invFth,detFth)
!
          call inv3x3(Fth_t,invFth_t,detFth_t)
!
!         Take out the thermal distortions from the total deformation
          F = matmul(dfgrd1,invFth)
!
          F_t = matmul(dfgrd0,invFth_t)
!
!
!
      end if
!     
!
!
!
!     MECHANICAL PART OF THE DEFORMATION GRADIENT
!     Calculate velocity gradient
!     Rate of deformation gradient
      Fdot = (F - F_t) / dt
!
!     Inverse of the deformation gradient
      call inv3x3(F,Finv,detF)
!
!     Velocity gradient
      L = matmul(Fdot,Finv)
!
!
!
!
!
!
!     CALCULATION OF TOTAL & MECHANICAL STRAINS      
!
!     Total stain increment from velocity gradient
      dstran33=(L+transpose(L))*0.5*dt
!
!
!
      call matvec6(dstran33,dstran)
!
!     Shear corrections
      dstran(4:6) = 2.0*dstran(4:6)
!
!
!     Total spin
      W=(L-transpose(L))*0.5
!     
!
!
!     Total spin increment - components
!     This is corrected as follows: Eralp - Alvaro 19.02.2023
!     The solution in Huang et al gives the negative -1/2*W
!     We obtained the spin directly from velocity gradient
!     1. It is positive
!     2. It has to be divided by 2
      domega(1) = W(1,2) - W(2,1)
      domega(2) = W(3,1) - W(1,3)
      domega(3) = W(2,3) - W(3,2)
      domega = domega * dt / 2.
!
!
!
!
!     Store the stress before rotation correction
      sigmarot_t = sigma_t
!
!
!     CALCULATION OF TRIAL STRESS
!
!     Trial stress
      sigmatr =  sigma_t + matmul(Cs,dstran)
!
!
!
!     3x3 stress tensor
      call vecmat6(sigma_t,sigma33_t)
!
!
!
!         Co-rotational stress
          sigma33_t = sigma33_t +
     + (matmul(W,sigma33_t) -
     + matmul(sigma33_t,W))*dt
!
!
!
!     Vectorize the initial guess
      call matvec6(sigma33_t,sigma_t)
!
!
!
!
!     CALCULATE RESOLVED SHEAR STRESS ON SLIP SYSTEMS
!     rss and its sign
      do is = 1, nslip
          tau_t(is) = dot_product(Schmidvec(is,:),sigma_t)
      end do
!
!
!
!     CALCULATE TRIAL-RESOLVED SHEAR STRESS ON SLIP SYSTEMS
!     rss and its sign
      do is = 1, nslip
          tautr(is) = dot_product(Schmidvec(is,:),sigmatr)
          abstautr(is) = abs(tautr(is))
      end do
!
!
!     maximum ratio of rss to crss
      maxx = maxval(abstautr/tauceff_t)
!
!
!
!
!     DECISION FOR USING CRYSTAL PLASTICITY
!     BASED ON THRESHOLD VALUE
!
!     Elastic solution
      if (maxx <= maxxcr) then
!
!
!
!
!         stress
          sigma = sigmatr
!
!         material tangent
          jacobi = Cs
!
!
!         Assign the global state variables
!         For NO SLIP condition
          totgammasum=totgammasum_t
          gammasum=gammasum_t
          gammadot=0.
          tauceff=tauceff_t
          tauc=tauc_t
          gnd=gnd_t
          ssd=ssd_t
          loop = loop_t
          ssdtot=ssdtot_t
          forest=forest_t
          substructure=substructure_t
          evmp=evmp_t
          plasdiss=plasdiss_t
          tausolute=tausolute_t
          Fp=Fp_t
          X=X_t
          cpconv=1
          cpconv0=0
!
!         Update orietnations
!         All the orientation changes are elastic - rotations
          dW33=0.
          dW33(1,2) = domega(1)
          dW33(1,3) = -domega(2)
          dW33(2,1) = -domega(1)
          dW33(2,3) = domega(3)
          dW33(3,1) = domega(2)
          dW33(3,2) = -domega(3)
!
          gmatinv = gmatinv_t + matmul(dW33,gmatinv_t)
!
!         Elastic strains in the crystal lattice
!         Add the former elastic strains
!
!         Undo shear corrections
          dummy6 = dstran
          dummy6(4:6) = 0.5*dummy6(4:6)
!
!         Convert the strain into matrix
          call vecmat6(dummy6,dummy33)
!
!         Elastic strains in the crystal reference
          dummy33_ = matmul(transpose(gmatinv),dummy33)
          dummy33 = matmul(dummy33_,gmatinv)
!
!
!         Vectorize
          call matvec6(dummy33,dummy6)  
!
!         Shear corrections
          dummy6(4:6) = 2.0*dummy6(4:6)
!
          Eec=Eec_t+dummy6
!
!
!
!
!     Solve using crystal plasticity
      else
! 
!
!
!         Guess if Forward Gradient Predictor scheme is not active
          if (predictor == 0) then
!
              sigma0 = (1.-phi)*sigma_t + phi*sigmatr
!
!
!
!         This part is added by Chris Hardie (11/05/2023)   
!         Former stress scheme
          elseif (predictor == 1) then
!
!
!
              if (dt_t > 0.0) then
!
              do i = 1, 6
                  sigma0(i) =
     + statev_sigma_t(noel,npt,i) +
     + (statev_sigma_t(noel,npt,i) -
     + statev_sigma_t2(noel,npt,i))*dt/dt_t
              end do
!
              else
                  sigma0 = sigma_t
              end if
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
          end if
!
!
!         CALCULATE RESOLVED SHEAR STRESS ON SLIP SYSTEMS
!         rss and its sign
          do is = 1, nslip
              tau0(is) = dot_product(Schmidvec(is,:),sigma0)
          end do
!
!

!
          call CP_Dunne(matid, phaid,
     + nslip, nscrew, mattemp, Cs, gf, G12,
     + burgerv, cubicslip, caratio,
     + Fp_t, gmatinv_t, Eec_t,
     + gammadot_t, gammasum_t,
     + totgammasum_t, evmp_t, plasdiss_t,
     + sigma0, tau0, sigmatr, 
     + forestproj, slip2screw,
     + Schmid_0, Schmid, 
     + Schmidvec, SchmidxSchmid,
     + smodel, sparam, cmodel, cparam,
     + imodel, iparam, hmodel, hparam,
     + bparam, sintmat1, sintmat2,
     + hintmat1, hintmat2,
     + tauceff_t, tauc_t, rhotot_t,
     + sumrhotot_t, ssdtot_t,
     + rhofor_t, forest_t, substructure_t,
     + gnd_t, ssd_t, loop_t, X_t,
     + dt, L, dstran,
     + Fp, gmatinv, Eec,
     + gammadot, gammasum,
     + totgammasum, evmp, plasdiss,
     + tauceff, tauc, tausolute,
     + ssdtot, ssd, loop, X,
     + forest, substructure,
     + sigma, jacobi, cpconv)
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
!         STILL NOT CONVERGED THE CUTBACK TIME
!         Convergence check
!         Use cpsolver did not converge
!         Diverged! - use time cut backs
          if (cpconv == 0) then         
!
!             Set the outputs to zero initially
              sigma = 0.!statev_sigma_t(noel,npt,:)
              jacobi = 0.!statev_jacobi_t(noel,npt,:,:)
!             Set time cut back and send a message
              pnewdt = cutback
!
          endif
!
!
!
!
      endif
!
!
!
!
!
!
!     Assign the global state variables          
      statev_gammasum(noel,npt,1:nslip)=gammasum
      statev_gammadot(noel,npt,1:nslip)=gammadot
      statev_tauc(noel,npt,1:nslip)=tauc
      statev_ssd(noel,npt,1:nslip)=ssd
      statev_loop(noel,npt,1:maxnloop)=loop
      statev_ssdtot(noel,npt)=ssdtot
      statev_forest(noel,npt,1:nslip)=forest
      statev_substructure(noel,npt)=substructure
      statev_evmp(noel,npt)=evmp
      statev_maxx(noel,npt)=maxx
      statev_totgammasum(noel,npt)=totgammasum
      statev_tausolute(noel,npt)=tausolute
      statev_sigma(noel,npt,1:6)=sigma
      statev_jacobi(noel,npt,1:6,1:6)=jacobi
      statev_Fp(noel,npt,1:3,1:3)=Fp
      statev_Fth(noel,npt,1:3,1:3)=Fth
      statev_Eec(noel,npt,1:6)=Eec
!     Crystal orientations at former time step
      statev_gmatinv(noel,npt,1:3,1:3)=gmatinv
!     Backstress
      statev_backstress(noel,npt,1:nslip)=X
!     Plastic dissipation power density
      statev_plasdiss(noel,npt)=plasdiss
!     Overall CRSS
      statev_tauceff(noel,npt,1:nslip)=tauceff
!
!     Write the outputs for post-processing
!     If outputs are defined by the user
      if (nstatv_outputs>0) then
!
          call assignoutputs(noel,npt,nstatv,statev)
!
      end if
!
!
!
!
!
!
      return
      end subroutine solve
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
!
!
!
!
!
!
!     Semi-implicit / Explicit state update rule
!     Solution using state variables at the former time step
      subroutine CP_Dunne(matid, phaid,
     + nslip, nscrew, mattemp, Cs, gf, G12,
     + burgerv, cubicslip, caratio,
     + Fp_t, gmatinv_t, Eec_t,
     + gammadot_t, gammasum_t,
     + totgammasum_t, evmp_t, plasdiss_t,
     + sigma0, tau0,
     + sigmatr, forestproj, slip2screw,
     + Schmid_0, Schmid,
     + Schmidvec, SchmidxSchmid,
     + slipmodel, slipparam,
     + creepmodel, creepparam,
     + irradiationmodel, irradiationparam,
     + hardeningmodel, hardeningparam,
     + backstressparam, 
     + sintmat1, sintmat2,
     + hintmat1, hintmat2,
     + tauceff_t, tauc_t, rhotot_t,
     + sumrhotot_t, ssdtot_t, rhofor_t,
     + forest_t, substructure_t,
     + gnd_t, ssd_t, loop_t, X_t,
     + dt, L, dstran,
     + Fp, gmatinv, Eec,
     + gammadot, gammasum,
     + totgammasum, evmp, plasdiss,
     + tauceff, tauc, tausolute,
     + ssdtot, ssd, loop, X,
     + forest, substructure,
     + sigma, jacobi, cpconv)
!
      use globalvariables, only : I3, I6, smallnum
!
      use userinputs, only : maxniter, maxnparam, maxnloop,
     + tauctolerance , SVDinversion,
     + backstressmodel, stateupdate, inversebackup
!
      use innerloop, only : Dunne_innerloop, Hardie_innerloop
!
      use utilities, only : vecmat6, matvec6,
     + nolapinverse, deter3x3, inv3x3, trace3x3,
     + vecmat9, matvec9, nolapinverse, SVDinverse
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
      use crss, only: slipresistance, totalandforest
!
      use errors, only : error
!
      implicit none
!
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
!     Plastic dissipation power density at the former time step
      real(8), intent(in) :: plasdiss_t
!     Cauchy stress guess
      real(8), intent(in) :: sigma0(6)
!     rss guess
      real(8), intent(in) :: tau0(nslip)
!     trial stress
      real(8), intent(in) :: sigmatr(6)
!     Forest projections
      real(8), intent(in) :: forestproj(nslip,nslip+nscrew)
!     Forest projections
      real(8), intent(in) :: slip2screw(nscrew,nslip)
!     Initial Schmid tensor
      real(8), intent(in) :: Schmid_0(nslip,3,3)  
!     Schmid tensor
      real(8), intent(in) :: Schmid(nslip,3,3)  
!     Vectorized Schmid tensor
      real(8), intent(in) :: Schmidvec(nslip,6) 
!     Schmid dyadic
      real(8), intent(in) :: SchmidxSchmid(nslip,6,6)
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
!     total dislocation density over all slip systems at the former time step
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
!     loop defect density per slip system at the former time step
      real(8), intent(in) :: loop_t(maxnloop)
!     backstress at former time step
      real(8), intent(in) :: X_t(nslip)
!     time increment
      real(8), intent(in) :: dt
!     total velocity gradient at the current time step
      real(8), intent(in) :: L(3,3)
!     mechanical strain increment
      real(8), intent(in) :: dstran(6)
!
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
!     Plastic dissipation power density at the current time step
      real(8), intent(out) :: plasdiss
!     Current values of state variables
!     overall crss
      real(8), intent(out) :: tauceff(nslip)
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
!     loop defect density per slip system at the current time step
      real(8), intent(out) :: loop(maxnloop)
!     crss at the current time step
      real(8), intent(out) :: X(nslip)
!     Cauchy stress
      real(8), intent(out) :: sigma(6)
!     material tangent
      real(8), intent(out) :: jacobi(6,6)  
!     convergence flag
      integer, intent(out) :: cpconv
!
!
!
!     Local variables used within this subroutine    
!
!     plastic velocity gradient for slip
      real(8) Lp_s(3,3), Dp_s(3,3)
!     plastic velocity gradient for creep
      real(8) Lp_c(3,3), Dp_c(3,3)
!     plastic velocity gradient
      real(8) Lp(3,3), Dp(3,3)
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
!     derivative of slip rates wrto crss for slip
      real(8) dgammadot_dtauc_s(nslip)
!     derivative of slip rates wrto crss for creep
      real(8) dgammadot_dtauc_c(nslip)
!
!
!     rss at the former time step
      real(8) :: tau(nslip)
!    absolute RSS and sign (used in Hardie scheme)
      real(8) :: abstau(nslip), signtau(nslip)
!
!     trial absolute RSS and sign (used in Hardie scheme)
      real(8) :: tautr(nslip), abstautr(nslip), signtautr(nslip)
!
!     Jacobian of the Newton-Raphson loop
!     and its inverse
      real(8)  :: dpsi_dsigma(6,6), invdpsi_dsigma(6,6)
!     residual of the Newton-Raphson loop
!     vector and scalar
      real(8) :: psinorm, psi(6)
!
!     Von-Mises equivalent plastic strain rate and increment
      real(8) :: pdot
!
!     stress increment
      real(8) :: dsigma(6)
!
!     stress 3x3 matrix
      real(8) :: sigma33(3,3)
!
!     plastic part of the deformation gradient
      real(8) :: detFp, invFp(3,3)
!
!     elastic part of the deformation gradient
      real(8) :: Fe(3,3)
!
!     elastic part of the velocity gradient
      real(8) :: Le(3,3)
!
!     elastic spin
      real(8) :: We(3,3)
!
!     increment in rotation matrix
      real(8) :: dR(3,3)
!
!     Von-Mises stress
      real(8) :: sigmaii, vms, sigmadev(3,3)
!
!     Co-rotational stress update
      real(8) :: dotsigma33(3,3)
!
!     Cauchy stress at former time step in 3x3
      real(8) :: sigma33_t(3,3)
!
!     Total mechanical strain increment
      real(8) :: dstran33(3,3)
!
!     plastic strain increment
      real(8) :: dstranp33(3,3)
!
!     elastic strain increment
      real(8) :: dstrane33(3,3)
!
!     crss increment
      real(8) :: dtauc(nslip)
!
!
!     ssd density increment
      real(8) :: dssd(nslip)
!
!     ssd density increment
      real(8) :: dloop(maxnloop)
!
!     backstress increment
      real(8) :: dX(nslip)    
!
!     total ssd density increment
      real(8) :: dssdtot
!
!     forest dislocation density increment
      real(8) :: dforest(nslip)
!
!     substructure dislocation density increment
      real(8) :: dsubstructure
!
!     Residues
      real(8) :: dtauceff(nslip), tauceff_old(nslip)
!
!
!     overall forest density
      real(8) :: rhofor(nslip)
!
!     overall total density
      real(8) :: rhotot(nslip)
!
!     overall total scalar density
      real(8) :: sumrhotot
!
!     Overall residue
      real(8) :: dtauceffnorm
!
!     error flag for svd inversion
      integer :: err
!
!     other variables
      real(8) :: dummy3(3), dummy33(3,3),
     + dummy33_(3,3), dummy6(6), dummy0
      integer :: is, il, iter, oiter
      integer :: iterinverse
!
!
!     Set convergence flag to "converged"
      cpconv = 1
!
!
!     Initial guess for NR scheme
!     Stress at the former time step
      sigma = sigma0
      tau = tau0
!
!
!     State assignments
      gammasum=gammasum_t
      tauc = tauc_t
      tauceff = tauceff_t
      ssdtot = ssdtot_t
      ssd = ssd_t
      loop = loop_t
      rhofor = rhofor_t
      X = X_t
      forest = forest_t
      substructure = substructure_t
!
!     Reset variables for the inner iteration    
      oiter = 0
      dtauceffnorm = 1.
!
!
!     Outer loop for state update
      do while ((dtauceffnorm >= tauctolerance)
     +.and.(oiter < maxniter).and.(cpconv==1))
!
!
!
!
          call Dunne_innerloop(
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
!
!         convergence check
          if (iter == maxniter) then
!             did not converge
              cpconv = 0
          end if
!
!
!         Check for NaN in the slip rate vector
          if(any(gammadot/=gammadot)) then
!             did not converge
              cpconv = 0
          endif
!
!         Check for Inf in the slip rate vector
          if(any(gammadot*0./=gammadot*0.)) then
!             did not converge
              cpconv = 0
          endif
!
!         Check for NaN in the stress vector
          if(any(sigma/=sigma)) then
!             did not converge
              cpconv = 0
          endif
!
!
!         Inverse scheme start here!!!
!         Try inverse scheme if not converged
          if ((cpconv==0).and.(inversebackup==1)) then
!
!             Reset convergence flag
              cpconv=1
!
!             Calculate trial-resolved shear stress on slip systems
!             Absolute value of trial-RSS and its sign
              do is = 1, nslip
                  tautr(is) = dot_product(Schmidvec(is,:),sigmatr)
                  signtautr(is) = sign(1.0,tautr(is))
                  abstautr(is) = abs(tautr(is))
              end do    
!
!             Initial guess for NR scheme
!             Stress at the former time step
              sigma = sigma0
              tau = tau0
!             Absolute value of RSS and its sign             
              do is = 1, nslip
                  signtau(is) = sign(1.0,tau0(is))
                  abstau(is) = abs(tau0(is))
              end do  
!
!             Reverse slip scheme
              call Hardie_innerloop(
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
!
!
!
!             Recalculate RSS on slip systems
!             RSS and its sign
              do is = 1, nslip
                  tau(is) = dot_product(Schmidvec(is,:),sigma)
              end do
!
!
!             Call the Dunne solver again
              call Dunne_innerloop(
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
!             assign jacobi and stress
              if (cpconv == 0) then
                  return
              end if
!
          end if
!         End of inverse solution
!
!
!
!         Calculate Von Mises invariant plastic strain rate
          pdot=sqrt(2./3.*sum(Dp*Dp))
!
!
!         Total plastic strain increment
          evmp = evmp_t + pdot*dt
!
!         Total slip over time per slip system
          gammasum = 0.
          do is = 1, nslip
!
              gammasum(is) = gammasum_t(is) +
     + gammadot(is)*dt
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
!         convergence check
          if (iter == maxniter) then
!             did not converge
              cpconv = 0
              return
          end if
!
!
!
!         Check for NaN in the slip rate vector
          if(any(gammadot/=gammadot)) then
!             did not converge
              cpconv = 0
              return
          endif
!
!         Check for Inf in the slip rate vector
          if(any(gammadot*0./=gammadot*0.)) then
!             did not converge
              cpconv = 0
              return
          endif
!
!         Check for NaN in the stress vector
          if(any(sigma/=sigma)) then
!             did not converge
              cpconv = 0
              return
          endif
!
!
!
!
!
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
     + rhofor_t,rhotot_t,substructure_t,
     + tausolute,dtauc,dssdtot,dforest,
     + dsubstructure,dssd,dloop)
!
!
!
!
!         Update the hardening states
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
!         Recalculate total and forest density
          call totalandforest(phaid,
     + nscrew, nslip, gnd_t,
     + ssd, ssdtot, forest,
     + forestproj, slip2screw, rhotot,
     + sumrhotot, rhofor)
!
!
!
!         Store the former value of effective tauc
          tauceff_old = tauceff
!
!         Recalculate slip resistance for the next iteration
!         Calculate crss
          call slipresistance(phaid, nslip,
     + gf, G12, burgerv, sintmat1, sintmat2,
     + tauc, rhotot, sumrhotot, rhofor,
     + substructure, tausolute, loop,
     + hardeningmodel, hardeningparam,
     + irradiationmodel, irradiationparam,
     + mattemp, tauceff)
!
!
!         Calculate the change of state
!         with respect to the former increment
          dtauceff = abs(tauceff-tauceff_old)
!
!         Explicit state update
          if (stateupdate==0) then
              dtauceffnorm = 0.
!         Semi-implicit state update
          elseif (stateupdate==1) then
              dtauceffnorm = maxval(dtauceff)
          end if
      
!
!
!
!         Check if the statevariables going negative due to softening
!         This may happen at high temperature and strain rates constants going bad
          if(any(tauc < 0.)) then
!             did not converge
              cpconv = 0
              return
          endif
!
          if(any(ssd < 0.)) then
!             did not converge
              cpconv = 0
              return
          endif
!
!         Loop density set to zero if negative
          do il = 1, maxnloop
              if(loop(il) < 0.) then
!
                  loop(il) = 0.
!
              endif
          enddo
!
!
          if(any(forest < 0.)) then
!             did not converge
              cpconv = 0
              return
          endif
!
          if(substructure < 0.) then
!             did not converge
              cpconv = 0
              return
          endif
!
!
!         Update backstress if local model is selected
          if (backstressmodel==1) then
!
              call backstressmodel1(backstressparam,
     + nslip,X,gammadot,dt,dX)
!
!             Update the backstress
              X = X_t + dX
!
          end if
!
!
!         increment iteration no.
          oiter = oiter + 1
!
!     End of state update (outer loop)
      end do
!
!
!
!     convergence check
      if (oiter == maxniter) then
!         did not converge
          cpconv = 0
          return
      end if
!
!
!
!
!     calculate jacobian
      jacobi = matmul(invdpsi_dsigma,Cs)    
!
!
!
!     Check for NaN in the jacobi matrix
      if(any(jacobi/=jacobi))  then
!         did not converge
          cpconv = 0
          return
      endif 
!
!
!
!
!     convert it to 3x3 marix
      call vecmat6(sigma,sigma33)
!
!     Trace of stress
      call trace3x3(sigma33,sigmaii)
!
!     deviatoric stress
      sigmadev = sigma33 - sigmaii*I3/3.
!
!     Von-Mises stress
      vms = sqrt(3./2.*(sum(sigmadev*sigmadev)))
!
!
!     variables for plastic part of the deformation gradient
      !dummy33 = I3 - Lp*dt
      !call inv3x3(dummy33,dummy33_,dummy0)
      dummy33_ = I3 + Lp*dt
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
          return
      else
!         Scale Fp with its determinant to make it isochoric
          Fp = Fp / detFp**(1./3.)
!
      end if
!
!
!
!
!
!
!     Elastic part of the velocity gradient
      Le = L - Lp
!
!     Elastic spin
      We = (Le - transpose(Le)) / 2.
!
!
!
!     stress rate due to spin
      dotsigma33 = matmul(We,sigma33) - matmul(sigma33,We)
!
!
!     Update co-rotational sress state
      sigma33 = sigma33 + dotsigma33*dt
!
!
!     Vectorize stress
      call matvec6(sigma33,sigma)
!
!
!
!
!
!     Orientation update  
!
!     Intermediate variable
      dR = I3 - We*dt
!
!     Invert or transpose since rotations are orthogonal
      dR = transpose(dR)
!
!
!
!     Update the crystal orientations
      gmatinv = matmul(dR, gmatinv_t)
!
!
!
!
!     Undo shear corrections
      dummy6 = dstran
      dummy6(4:6) = 0.5*dummy6(4:6)
!
!     Convert the strain into matrix
      call vecmat6(dummy6,dstran33)
!
!     Elastic strain increment
      dstrane33 = dstran33-dstranp33
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
      Eec=Eec_t+dummy6
!
!
!     Plastic dissipation power density
      plasdiss=0.
      do is = 1, nslip
!
          plasdiss = plasdiss +
     + tau(is)*gammadot(is)*dt
!
      enddo
!
!
!     Sum over the fomer value
      plasdiss = plasdiss_t + plasdiss
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
      return
      end subroutine CP_Dunne
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
!
!
!     Implicit state update rule
!     Solution using the updated state variables
      subroutine rotateslipsystems(iphase,nslip,caratio,
     + gmatinv,dirc,norc,dirs,nors)
      implicit none
      integer, intent(in) :: iphase
      integer, intent(in) :: nslip
      real(8), intent(in) :: caratio
      real(8), intent(in) :: gmatinv(3,3)
      real(8), intent(in) :: dirc(nslip,3)
      real(8), intent(in) :: norc(nslip,3)
      real(8), intent(out) :: dirs(nslip,3)
      real(8), intent(out) :: nors(nslip,3)
!
      integer :: i, is
      real(8) :: tdir(3), tnor(3)
      real(8) :: tdir1(3), tnor1(3)
      real(8) :: dirmag, normag
!
      dirs=0.;nors=0.
!
!
      do is=1,nslip ! rotate slip directions 
!
          tdir = dirc(is,:)
          tnor = norc(is,:)
!
!
!
          tdir1 = matmul(gmatinv,tdir)
          tnor1 = matmul(gmatinv,tnor)
!
          dirmag = norm2(tdir1)
          normag = norm2(tnor1)
!
!
          dirs(is,:) = tdir1/dirmag
          nors(is,:) = tnor1/normag
!
!
      end do
!
!
!
!
      return
      end subroutine rotateslipsystems
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
!
      end module cpsolver