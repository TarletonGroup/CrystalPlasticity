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
     + statev_gmatinv, statev_gmatinv_t, statev_gmatinv_0,
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
     + numslip_all, numscrew_all, phaseid_all,
     + dirc_0_all, norc_0_all, caratio_all, cubicslip_all,
     + Cc_all, gf_all, G12_all, v12_all, alphamat_all, burgerv_all,
     + sintmat1_all, sintmat2_all, hintmat1_all, hintmat2_all,
     + slipmodel_all, slipparam_all, creepmodel_all, creepparam_all,
     + hardeningmodel_all, hardeningparam_all, irradiationmodel_all,
     + irradiationparam_all, backstressparam_all, slip2screw_all,
     + statev_backstress_t, statev_backstress,
     + I3, I6, smallnum, largenum
!
      use userinputs, only: constanttemperature, temperature,
     + predictor, maxnslip, maxnparam, maxxcr, maxnloop, 
     + cutback, phi, stateupdate
!
!
      use usermaterials, only: materialparam
!
      use crss, only: slipresistance
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
!     Flag for Euler solver convergence
      integer :: cpconv0
!
!     Local state variables with known dimensions
!     crystal to sample transformation initially
      real(8) :: gmatinv_0(3,3)
!     crystal to sample transformation at the former time step
      real(8) :: gmatinv_t(3,3)
!     crystal to sample transformation at the former time step
      real(8) :: gmatinv(3,3), gmatinv0(3,3)
!     stress at the former time step
      real(8) :: sigma_t(6)
!     rss/crsss ratio at the former time step
      real(8) :: maxx_t
!     rss/crsss ratio at the current time step 
      real(8) :: maxx
!     elastic strains in the crystal reference at the former time step
      real(8) :: Eec_t(6)
!     elastic strains in the crystal reference at the current time step
      real(8) :: Eec(6), Eec0(6)
!     plastic deformation gradient at the former time step
      real(8) :: Fp_t(3,3)
!     plastic deformation gradient at the current time step
      real(8) :: Fp(3,3), Fp0(3,3)
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
      real(8) :: gammasum(numslip_all(matid)),
     + gammasum0(numslip_all(matid))
!     slip rates
      real(8) :: gammadot_t(numslip_all(matid))
      real(8) :: gammadot(numslip_all(matid)),
     + gammadot0(numslip_all(matid))
!     slip resistance
      real(8) :: tauc_t(numslip_all(matid))
      real(8) :: tauc(numslip_all(matid)),
     + tauc0(numslip_all(matid))
!     effective overall slip resistance
!     (tauc0 + GND + SSD + solute + substructure + forest + etc.)
      real(8) :: tauceff_t(numslip_all(matid))
!     Note the size of GND is different
!     gnd density (nslip + nscrew)
      real(8) :: gnd_t(numslip_all(matid)+numscrew_all(matid))
      real(8) :: gnd(numslip_all(matid)+numscrew_all(matid))
!
!     ssd density (nslip)
      real(8) :: ssd_t(numslip_all(matid))
      real(8) :: ssd(numslip_all(matid)),
     + ssd0(numslip_all(matid))
!     loop density (maxnloop)
      real(8) :: loop_t(maxnloop)
      real(8) :: loop(maxnloop), loop0(maxnloop)
!     total forest dislocation density - derived from other terms
      real(8) :: rhofor_t(numslip_all(matid))
      real(8) :: rhofor(numslip_all(matid)),
     + rhofor0(numslip_all(matid))
!     total density
      real(8) :: rhotot_t(numslip_all(matid))
      real(8) :: rhotot(numslip_all(matid))
!     forest dislocation density as a state variable
      real(8) :: forest_t(numslip_all(matid))
      real(8) :: forest(numslip_all(matid)),
     + forest0(numslip_all(matid))
!
!
!     Scalar state variables
!     equivalent Von-Mises plastic strain
      real(8) :: evmp_t
      real(8) :: evmp, evmp0
!
!     cumulative slip
      real(8) :: totgammasum_t
      real(8) :: totgammasum, totgammasum0
!     solute strength
      real(8) :: tausolute_t
      real(8) :: tausolute, tausolute0
!     substructure density
      real(8) :: substructure_t
      real(8) :: substructure, substructure0
!     total density
      real(8) :: ssdtot_t
      real(8) :: ssdtot, ssdtot0
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
!     slip direction and slip plane normal at the sample reference (undeformed)
      real(8) :: dirs_0(numslip_all(matid),3)
      real(8) :: nors_0(numslip_all(matid),3)
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
      real(8) ::  sigmatr(6)
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
!     BACKUP values - FG method
      real(8) ::  sigma0(6), jacobi0(6,6)
!     value of resolved shear stress
      real(8) :: tau0(numslip_all(matid))
!
!
!     value of resolved shear stress
      real(8) :: tau_t(numslip_all(matid))
!
!     value of trial resolved shear stress
      real(8) :: tautr(numslip_all(matid))
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
!
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
      totgammasum_t = statev_totgammasum_t(noel,npt)
      tausolute_t = statev_tausolute_t(noel,npt)
      sigma_t = statev_sigma_t(noel,npt,:)
      Fp_t = statev_Fp_t(noel,npt,:,:)
      Fth_t = statev_Fth_t(noel,npt,:,:)
      Eec_t = statev_Eec_t(noel,npt,:)
!     Crystal orientations at former time step
      gmatinv_t = statev_gmatinv_t(noel,npt,:,:)
!     Crystal orientations initially
      gmatinv_0 = statev_gmatinv_0(noel,npt,:,:)
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
      v12 = v12_all(matid)
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
     + notused2,notused3, !forest and substructure also need to be added
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
!     Slip directions in the undeformed sample reference
      call rotateslipsystems(phaid,nslip,caratio,
     + gmatinv_0,dirc_0,norc_0,dirs_0,nors_0)
!
!
!     Calculate undeformed Schmid tensor
      Schmid_0=0.
      do is=1,nslip
!
!         Slip direction
          sdir = dirs_0(is,:)
!         Slip plane normal
          ndir = nors_0(is,:)
!
          do i=1,3
              do j=1,3
                  SNij(i,j) = sdir(i)*ndir(j)
                  Schmid_0(is,i,j) = SNij(i,j)
              enddo
          enddo
!
      end do
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
                  NSij(i,j) = ndir(i)*sdir(j)
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
     + gnd_t, ssd_t, ssdtot_t, forest_t,
     + forestproj, slip2screw,
     + rhotot_t, sumrhotot_t, rhofor_t)
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
      end do
!
!
!     maximum ratio of rss to crss
      maxx = maxval(abs(tautr)/tauceff_t)
!
!
!
!
!     DECISION FOR USING CRYSTAL PLASTICITY
!     BASED ON THRESHOLD VALUE
!
!     Elastic solution
      if (maxx < maxxcr) then
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
          tauc=tauc_t
          
          gnd=gnd_t
          ssd=ssd_t
          loop = loop_t
          ssdtot=ssdtot_t
          forest=forest_t
          substructure=substructure_t
          evmp=evmp_t
          tausolute=tausolute_t
          Fp=Fp_t
          X=X_t
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
!         Else Forward Gradient Predictor scheme computes sigma0
          elseif (predictor == 1) then
!
!
!
! 
!
              call CP_Huang_linear(
     + matid, phaid, nslip, nscrew,
     + mattemp, Cs, gf, G12,
     + burgerv, cubicslip, caratio,
     + tau_t, Fp_t, gmatinv_t, Eec_t,
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
!
!         This part is added by Chris Hardie (11/05/2023)   
!         Former stress scheme
          elseif (predictor == 2) then
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
          end if
!
      
!         CALCULATE RESOLVED SHEAR STRESS ON SLIP SYSTEMS
!         rss and its sign
          do is = 1, nslip
              tau0(is) = dot_product(Schmidvec(is,:),sigma0)
          end do
!
!         Reset convergence flag
          cpconv = 0
!
!         Explicit time integration of states
          if (stateupdate==0) then
!
!
!
!             Solve using explicit crystal plasticity solver
              call CP_Dunne_explicit(matid, phaid,
     + nslip, mattemp, Cs, gf, G12,
     + burgerv, cubicslip, caratio,
     + Fp_t, gmatinv_t, Eec_t,
     + gammadot_t, gammasum_t,
     + totgammasum_t, evmp_t,
     + sigma0, tau0, sigmatr,
     + Schmid_0, Schmid, 
     + Schmidvec, SchmidxSchmid,  
     + smodel, sparam, cmodel, cparam,
     + imodel, iparam, hmodel, hparam,
     + bparam, sintmat1, sintmat2,
     + hintmat1, hintmat2,
     + tauceff_t, tauc_t,
     + rhotot_t, sumrhotot_t,
     + ssdtot_t, rhofor_t,
     + forest_t, substructure_t,
     + ssd_t, loop_t, X_t,
     + dt, L, dstran,
     + Fp, gmatinv, Eec,
     + gammadot, gammasum,
     + totgammasum, evmp,
     + tauc, tausolute,
     + ssdtot, ssd, loop, X,
     + forest, substructure,
     + sigma, jacobi, cpconv)
!
!
!
!         Semi-implicit time integration of states
          elseif (stateupdate==1) then
!
!
              call CP_Dunne_semiimplicit(matid, phaid,
     + nslip, nscrew, mattemp, Cs, gf, G12,
     + burgerv, cubicslip, caratio,
     + Fp_t, gmatinv_t, Eec_t,
     + gammadot_t, gammasum_t,
     + totgammasum_t, evmp_t,
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
     + totgammasum, evmp,
     + tauc, tausolute,
     + ssdtot, ssd, loop, X,
     + forest, substructure,
     + sigma, jacobi, cpconv)
!
!
!
!         Fully-implicit time integration of states
          elseif (stateupdate==2) then
!
!
!             Fully implicit state update
              call CP_Huang_fullyimplicit(matid, phaid,
     + nslip, nscrew, mattemp, Cs, gf, G12,
     + burgerv, cubicslip, caratio,
     + tau_t, Fp_t, gmatinv_t, Eec_t,
     + gammadot_t, gammasum_t,
     + totgammasum_t, evmp_t, sigmarot_t,
     + forestproj, slip2screw,
     + dirc_0, norc_0,
     + Schmidvec, Schmid_0,
     + Schmid, SchmidxSchmid,
     + smodel, sparam,
     + cmodel, cparam,
     + imodel, iparam,
     + hmodel, hparam,
     + bparam, sintmat1, sintmat2,
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
!
!
          endif
!
!
!
!
!
!         If stress-based Dunne crystal plasticity DOES NOT CONVERGE
          if (cpconv == 0) then
!             If Forward Gradient Predictor solution available
!             Use the Eulersolver result (if turned ON)
              if (predictor == 1) then
!
!
!                 If Forward Gradient Predictor solution converged
                  if (cpconv0 == 1) then
!
!                     USE Forward Gradient Predictor SOLUTION AS THE SOLUTION
                      Fp = Fp0
                      gmatinv = gmatinv0
                      Eec = Eec0
                      gammadot = gammadot0
                      gammasum = gammasum0
                      totgammasum = totgammasum0
                      evmp = evmp0
                      tauc= tauc0

                      tausolute = tausolute
                      ssdtot = ssdtot0
                      ssd = ssd0
                      loop = loop0
                      X = X0
                      forest = forest0
                      substructure = substructure0
                      sigma = sigma0
                      jacobi =jacobi0
                      cpconv = cpconv0
!
!
!
                  endif
!
              endif
!
!
!
          endif
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
              call error(14)
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
!
!
!
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
!     Explicit state update rule
!     Solution using state variables at the former time step
      subroutine CP_Dunne_explicit(matid, phaid,
     + nslip, mattemp, Cs, gf, G12,
     + burgerv, cubicslip, caratio,
     + Fp_t, gmatinv_t, Eec_t,
     + gammadot_t, gammasum_t,
     + totgammasum_t, evmp_t,
     + sigma0, tau0, sigmatr,
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
     + sumrhotot_t, ssdtot_t,
     + rhofor_t, forest_t, substructure_t,
     + ssd_t, loop_t, X_t,
     + dt, L, dstran,
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
      use userinputs, only : maxniter, tolerance,
     + maxnparam, maxnloop, SVDinversion,
     + backstressmodel, crit
!
      use utilities, only : vecmat6, matvec6,
     + nolapinverse, deter3x3, inv3x3, trace3x3,
     + vecmat9, matvec9, nolapinverse,
     + SVDinverse
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
!     Cauchy stress guess
      real(8), intent(in) :: sigma0(6)
!     rss guess
      real(8), intent(in) :: tau0(nslip)
!     trial stress
      real(8), intent(in) :: sigmatr(6)
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
      real(8), intent(in) :: ssd_t(nslip)
!     defect loop density per slip system at the former time step
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
!     Local variables used within this subroutine
!
!     plastic velocity gradient / stretch rate for slip
      real(8) Lp_s(3,3), Dp_s(3,3)
!     plastic velocity gradient / stretch rate for creep
      real(8) Lp_c(3,3), Dp_c(3,3)
!     total plastic velocity gradient / stretch rate
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
!     derivative of slip rates wrto rss for slip
      real(8) dgammadot_dtauc_s(nslip)
!     derivative of slip rates wrto rss for creep
      real(8) dgammadot_dtauc_c(nslip)
!
!     rss at the former time step
      real(8) :: tau(nslip)
!
!     Jacobian of the Newton-Raphson loop
!     and its inverse
      real(8)  :: dpsi_dsigma(6,6), invdpsi_dsigma(6,6)
!     residual of the Newton-Raphson loop
!     vector and scalar
      real(8) :: psinorm, psi(6)
!
!     plastic strain increment
      real(8) :: dstranp33(3,3), dstranp(6)
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
!     Total mechanical strain increment
      real(8) :: dstran33(3,3)
!
!     elastic strain increment
      real(8) :: dstrane33(3,3)
!
!     crss increment
      real(8) :: dtauc(nslip)
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
!     threshold for stress increment
      real(8) :: dsigma_cr
!
!     error flag for svd inversion
      integer :: err
!
!
!     derivative of tauc with respect to slip rate
      real(8) :: ddtauc_ddgamma(nslip,nslip)
!
!     other variables
      real(8) :: dummy3(3), dummy33(3,3),
     + dummy33_(3,3), dummy6(6), dummy0
      integer :: is, il, iter, i
!   
!     
!
!
!
!     Set convergence flag
      cpconv = 1
!
!
!     Reset variables for iteration   
      iter = 0
      psinorm = 1.
!
!
!     Initial guess for NR scheme
!     Stress at the former time step
      sigma = sigma0
      tau = tau0
!
!
!     Threshod calculation
      dsigma_cr = crit*sum(tauceff_t(1:nslip))/nslip
!
!     Newton-Raphson (NR) iteration to find stress increment
      do while ((psinorm > tolerance).and.(iter < maxniter))
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
     + tau,X_t,tauceff_t,rhofor_t,burgerv,dt,
     + nslip,phaid,mattemp,slipparam,
     + irradiationmodel,irradiationparam,
     + cubicslip,caratio,Lp_s,Dp_s,Pmat_s,
     + gammadot_s,dgammadot_dtau_s,
     + dgammadot_dtauc_s)
!
!
!         double exponent law (exponential law)
          elseif (slipmodel == 2) then
!
!
             call doubleexpslip(Schmid_0,Schmid,SchmidxSchmid,
     + tau,X_t,tauceff_t,burgerv,dt,nslip,phaid,
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
     + tau,X_t,tauceff_t,burgerv,dt,
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
              call expcreep(Schmid_0,Schmid,SchmidxSchmid,
     + tau,X_t,tauceff_t,dt,nslip,phaid,
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
          Lp = Lp_s + Lp_c
          Dp = Dp_s + Dp_c
          Pmat = Pmat_s + Pmat_c
          gammadot = gammadot_s + gammadot_c
!
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
!         Check infinity
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
!         Then invert (double precision version)
          call nolapinverse(dpsi_dsigma,invdpsi_dsigma,6)
!
!
!
!
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
!
!                 Invert
                  call SVDinverse(dpsi_dsigma,6,invdpsi_dsigma,err)
!
              else
!
                  err = 1
!
              endif
!
!
!
!
              if (err==1) then
!                 did not converge
                  cpconv = 0
                  return
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
          endif
!
!
!
!
!
!
!
!
!
!         residual (predictor - corrector)
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
!         stress increment threshold
!         stress component checks
          do i=1,6
              if (dabs(dsigma(i)).gt.dsigma_cr) then
                  dsigma(i) = dsign(1.0,dsigma(i))*dsigma_cr
              endif
          enddo
!
!         stress update
          sigma = sigma + dsigma
!
!         convert it to 3x3 marix
          call vecmat6(sigma,sigma33)
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
!     End of NR iteration
      end do
!
!
!
!
!
!
!
!
!     convergence check
      if (iter == maxniter) then
!         did not converge
          cpconv = 0
          return
      end if
!
!
!     calculate jacobian
      jacobi = matmul(invdpsi_dsigma,Cs)
!
!
!
!     Check for NaN in the stress vector
      if(any(sigma/=sigma)) then
!         did not converge
          cpconv = 0
          return
      endif
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
!     calculate von mises invariant plastic strain rate
      pdot=sqrt(2./3.*sum(Dp*Dp))
!
!
!     Von-Mises equivalent total plastic strain
      evmp = evmp_t + pdot*dt
!
!     Total slip over time per slip system
      gammasum = 0.
      do is =1, nslip
!
          gammasum(is) = gammasum_t(is) +
     + gammadot(is)*dt
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
!
!
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
!      dummy33 = I3 - Lp*dt
!      call inv3x3(dummy33,dummy33_,dummy0)
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
!
!
!
!
!
!
!
!     Update the states using hardening laws
       call hardeningrules(phaid,nslip,
     + mattemp,dt,G12,burgerv,
     + totgammasum,gammadot,pdot,
     + irradiationmodel,irradiationparam,
     + hardeningmodel,hardeningparam,
     + hintmat1,hintmat2,
     + tauc_t,ssd_t,loop_t,
     + rhofor_t,rhotot_t,substructure_t,
     + tausolute,dtauc,dssdtot,dforest,
     + dsubstructure,dssd,dloop,
     + ddtauc_ddgamma)
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
!
!
!
!     Check if the statevariables going negative due to softening
!     This may happen at high temperature and strain rates constants going bad
      if(any(tauc < 0.)) then
!         did not converge
          cpconv = 0
          return
      endif
!
      if(any(ssd < 0.)) then
!         did not converge
          cpconv = 0
          return
      endif
!
!
!     Loop density set to zero if negative
      do il = 1, maxnloop
          if(loop(il) < 0.) then
!
              loop(il) = 0.
!
          endif
      enddo
!
!
!
!
      if(any(forest < 0.)) then
!         did not converge
          cpconv = 0
          return
      endif
!
      if(substructure < 0.) then
!         did not converge
          cpconv = 0
          return
      endif
!
!
      X = X_t
!     Update backstress if local AF model is selected
      if (backstressmodel==1) then
!
          call backstressmodel1(backstressparam,
     + nslip,X_t,gammadot,dt,dX)
!
!         Update the backstress
          X = X_t + dX
!
      end if
!
      return
      end subroutine CP_Dunne_explicit
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
!     Solution using state variables at the former time step
      subroutine CP_Dunne_semiimplicit(matid, phaid,
     + nslip, nscrew, mattemp, Cs, gf, G12,
     + burgerv, cubicslip, caratio,
     + Fp_t, gmatinv_t, Eec_t,
     + gammadot_t, gammasum_t,
     + totgammasum_t, evmp_t,
     + sigma0, tau0,
     + sigmatr, forestproj, 
     + slip2screw,
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
     + totgammasum, evmp,
     + tauc, tausolute,
     + ssdtot, ssd, loop, X,
     + forest, substructure,
     + sigma, jacobi, cpconv)
!
      use globalvariables, only : I3, I6, smallnum
!
      use userinputs, only : 
     + maxniter, maxnparam, maxnloop,
     + tolerance, tauctolerance, SVDinversion,
     + backstressmodel, crit
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
      use crss, only: slipresistance
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
!
!     Jacobian of the Newton-Raphson loop
!     and its inverse
      real(8)  :: dpsi_dsigma(6,6), invdpsi_dsigma(6,6)
!     residual of the Newton-Raphson loop
!     vector and scalar
      real(8) :: psinorm, psi(6)
!
!     plastic strain increment
      real(8) :: dstranp33(3,3), dstranp(6)
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
!     Current values of state variables
!     overall crss
      real(8) :: tauceff(nslip)
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
!     threshold for stress increment
      real(8) :: dsigma_cr
!
!     error flag for svd inversion
      integer :: err
!
!     derivative of tauc with respect to slip rate
      real(8) :: ddtauc_ddgamma(nslip,nslip)
!
!     other variables
      real(8) :: dummy3(3), dummy33(3,3),
     + dummy33_(3,3), dummy6(6), dummy0
      integer :: is, il, iter, oiter, i
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
!     Threshod calculation
      dsigma_cr = crit * sum(tauceff_t(1:nslip))/nslip
!
!     Outer loop for state update
      do while ((dtauceffnorm > tauctolerance).and.(oiter < maxniter))
!
!
!         Reset variables for the inner iteration
          psinorm = 1.
          iter = 0
!
!         Newton-Raphson (NR) iteration to find stress increment
          do while ((psinorm > tolerance).and.(iter < maxniter))
!
!
!             Slip models to find slip rates
!
!             none
              if (slipmodel == 0) then
!
                  Lp_s = 0.
                  Dp_s = 0.
                  Pmat_s = 0.
                  gammadot_s = 0.
!
!
!             sinh law
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
!             exponential law
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
!             power law
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
!             Slip due to creep
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
!             Sum the effects of creep and slip rates
              Lp = Lp_s + Lp_c
              Dp = Dp_s + Dp_c
              Pmat = Pmat_s + Pmat_c
              gammadot = gammadot_s + gammadot_c
!
!
!
!             Check for NaN in the slip rates
              if(any(gammadot/=gammadot)) then
!                 did not converge
                  cpconv = 0
                  return
              endif
!
!
!             Check infinity
              if(any(gammadot*0./=gammadot*0.)) then
!                 did not converge
                  cpconv = 0
                  return
              endif
!
!
!             Check for the Pmat
              if(any(Pmat /= Pmat)) then
!                 did not converge
                  cpconv = 0
                  return
              endif
!
!
!
!             Plastic strain increment
              dstranp33 = Dp*dt
              call matvec6(dstranp33,dstranp)
              dstranp(4:6)=2.*dstranp(4:6)
!
!
!
!
!             Tangent-stiffness calculation
!             Jacobian of the Newton loop (see Dunne, Rugg, Walker, 2007)
              dpsi_dsigma = I6 + matmul(Cs, Pmat)
!
!
!
!             Then invert (double precision version)
              call nolapinverse(dpsi_dsigma,invdpsi_dsigma,6)
!
!
!
!             If inversion is not successfull!
!             Check for the inverse
              if(any(invdpsi_dsigma /= invdpsi_dsigma)) then
!
!                 Try using singular value decomposition
!                 If singular value decomposition is ON
                  if (SVDinversion==1) then
!
!                     Invert
                      call SVDinverse(dpsi_dsigma,6,invdpsi_dsigma,err)
!
!
!
                  else
!
!
!                     did not converge
                      err = 1
!
!
!
                  end if
!
!                 Check again and if still not successfull
                  if(err==1) then
!                     did not converge
                      cpconv = 0
                      return
                  end if
!
!
              endif
!
!
!
!             residual (predictor - corrector scheme)
              psi = sigmatr - sigma - matmul(Cs,dstranp)
!
!             norm of the residual
              psinorm = sqrt(sum(psi*psi))
!
!
!             stress increment
              dsigma = matmul(invdpsi_dsigma,psi)
!
!             stress increment threshold
!             stress component checks
              do i=1,6
                  if (dabs(dsigma(i)).gt.dsigma_cr) then
                      dsigma(i) = dsign(1.0,dsigma(i))*dsigma_cr
                  endif
              enddo
!
!             stress update
              sigma = sigma + dsigma
!
!             convert it to 3x3 marix
              call vecmat6(sigma,sigma33)
!
!
!             calculate resolved shear stress on slip systems
!             rss and its sign
              do is = 1, nslip
                  tau(is) = dot_product(Schmidvec(is,:),sigma)
              end do     
!
!
!             increment iteration no.
              iter = iter + 1
!
!         End of NR iteration (inner loop)
          end do
!
!
!
!
!
!
!
!
!
!         calculate Von Mises invariant plastic strain rate
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
     +    sum(abs(gammadot))*dt
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
!         convergence check
          if (iter == maxniter) then
!             did not converge
              cpconv = 0
              return
          end if
!
!
!
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
     + tauc,ssd,loop,
     + rhofor,rhotot,substructure,
     + tausolute,dtauc,dssdtot,dforest,
     + dsubstructure,dssd,dloop,
     + ddtauc_ddgamma)
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
     + forestproj, slip2screw,
     + rhotot, sumrhotot, rhofor)
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
          dtauceffnorm = maxval(dtauceff)
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
!
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
      return
      end subroutine CP_Dunne_semiimplicit
!
!
!
!
!
!
!
!     Forward Gradient Predictor scheme
      subroutine CP_Huang_linear(matid, phaid,
     + nslip, nscrew, mattemp, Cs, gf, G12,
     + burgerv, cubicslip, caratio,
     + tau_t, Fp_t, gmatinv_t, Eec_t,
     + gammadot_t, gammasum_t,
     + totgammasum_t, evmp_t, sigma_t,
     + forestproj, slip2screw, dirs_t, nors_t,
     + Schmidvec, Schmid_0,
     + Schmid, SchmidxSchmid,
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
      use globalvariables, only : I3, I6, 
     + smallnum, largenum
!
      use userinputs, only : theta, maxnparam, maxnloop,
     + backstressmodel, crit
!
      use utilities, only : vecmat6, matvec6,
     + nolapinverse, deter3x3, inv3x3, trace3x3,
     + vecmat9, matvec9, SVDgeninverse
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
!     RSS at former time step
      real(8), intent(in) :: tau_t(nslip)
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
!     undeformed slip direction in crsytal reference frame
      real(8), intent(in) :: dirs_t(nslip,3)
!     undeformed slip plane normal in crystal reference frame
      real(8), intent(in) :: nors_t(nslip,3)
!     Schmid tensor - vectorized
      real(8), intent(in) :: Schmidvec(nslip,6)
!     Initial Schmid tensor - unused
      real(8), intent(in) :: Schmid_0(nslip,3,3)  
!     Schmid tensor - unused
      real(8), intent(in) :: Schmid(nslip,3,3)
!     Schmid dyadic - unused
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
      real(8) :: tau(nslip), dtau(nslip)
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
      real(8) :: Pmat_s(6,6), Pmat_c(6,6)
      real(8) :: gammadot_s(nslip), gammadot_c(nslip)    
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
      real(8) :: dgamma(nslip)
!     Total spin increment
      real(8) :: domega33(3,3)
!     Plastic spin increment
      real(8) :: domegap33(3,3)
!     Elastic spin increment
      real(8) :: domegae33(3,3)
!     Rotation matrix increment
      real(8) :: dgmatinv(3,3)
!
!     Increment in Cauchy stress
      real(8) :: dsigma(6)
!
!
!     Determinant of plastic deformation gradient
      real(8) :: detFp
!     Strain increment and related variables
      real(8) :: dstranp(6), dstrane(6)
      real(8) :: dstrane33(3,3), dstranp33(3,3)
!     crss increment
      real(8) :: dtauc(nslip)
!     ssd density increment
      real(8) :: dssd(nslip)
!     loop density increment
      real(8) :: dloop(maxnloop)
!     backstress increment
      real(8) :: dX(nslip)
!     total ssd density increment
      real(8) :: dssdtot
!     forest dislocation density increment
      real(8) :: dforest(nslip)
!     substructure dislocation density increment
      real(8) :: dsubstructure
!
!
!
!
!     overall crss
      real(8) :: tauceff(nslip)
!     forest density
      real(8) :: rhofor(nslip)
!     total density
      real(8) :: rhotot(nslip)
!     total scalar density
      real(8) :: sumrhotot
!
!     Derivative of tauc
      real(8) :: ddtauc_ddgamma(nslip,nslip)
!
!     Linear solution
      real(8) :: Amat(nslip+1,6), invAmat(6,nslip+1)
      real(8) :: bvec(nslip+1), dsigma_cr
      integer :: err, flag
!
!     Dummy variables
      real(8) :: dummy33(3,3), dummy33_(3,3)
      real(8) :: dummy0, dummy6(6), dummy66(6,6),
     + dummynslip(nslip,nslip) 
!
!     Variables uses within the subroutine
      integer :: is, js, il, i, j, k, a, b
!
!
!
!     Initiate convergence flag
      cpconv = 1
!
!     Threshod calculation
      dsigma_cr = crit*sum(tauceff_t(1:nslip))/nslip
!
!
!     Volumetric change in the strain
      dev = dstran(1) + dstran(2) + dstran(3)
!
!
!     Plastic spin dyadic
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
!
!
!
!     Calculate beta and ddemsd
      beta=0.; ddemsd=0.
      do is = 1, nslip
!
!         Symbolic math result
          beta(1,is) = -2.*W(2,is)*sigma_t(5) + 2.*W(1,is)*sigma_t(4)
          beta(2,is) = -2.*W(1,is)*sigma_t(4) + 2.*W(3,is)*sigma_t(6)
          beta(3,is) = -2.*W(3,is)*sigma_t(6) + 2.*W(2,is)*sigma_t(5)
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
!     Slip models to find slip rates
!
!     none
      if (slipmodel == 0) then
!
          Lp_s = 0.
          Dp_s = 0.
          Pmat_s = 0.
          gammadot_s = 0.
          dgammadot_dtau_s = 0.
          dgammadot_dtauc_s = 0.
!
!     sinh law
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
!     exponential law
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
!     power law
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
!     Slip due to creep     
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
!     Sum the effects of creep and slip rates
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
!     Check for the slip rates
      if(any(gammadot /= gammadot)) then
!         did not converge
          cpconv = 0
          return
      endif      
!
!     Check infinity
      if(any(gammadot*0./=gammadot*0.)) then
!         did not converge
          cpconv = 0
          return
      endif
!
!     Plastic strain-related quantities for hardening calculations
!
!
!
!     calculate von mises invariant plastic strain rate
      pdot=sqrt(2./3.*sum(Dp*Dp))
!
!
!
!     Total slip over time per slip system
      gammasum = 0.
      do is =1, nslip
!
          gammasum(is) = gammasum_t(is) +
     + gammadot(is)*dt
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
!
      X = X_t
!     Update backstress
      if (backstressmodel==1) then
!
          call backstressmodel1(backstressparam,
     + nslip,X_t,gammadot,dt,dX)
!
!         Update the backstress
          X = X_t + dX
!
      end if
!

!
!
!     Update the states using hardening laws
       call hardeningrules(phaid,nslip,
     + mattemp,dt,G12,burgerv,
     + totgammasum,gammadot,pdot,
     + irradiationmodel,irradiationparam,
     + hardeningmodel,hardeningparam,
     + hintmat1,hintmat2,
     + tauc_t,ssd_t,loop_t,
     + rhofor_t,rhotot_t,substructure_t,
     + tausolute,dtauc,dssdtot,dforest,
     + dsubstructure,dssd,dloop,
     + ddtauc_ddgamma)
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
     + ddtauc_ddgamma(a,b)*dgammadot_dtauc(a)*theta*dt
!     + *sign(1.0,gammadot(b))
!
!
          end do
!
          Nab(a,a) = Nab(a,a) + 1.
!
!         Given quantities in vector form
          dgamma(a) = dot_product(ddemsd(1:6,a),dstran(1:6))*
     + dgammadot_dtau(a)*theta*dt + gammadot(a)*dt
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
!
!     Stress incremenet
      dsigma = matmul(Cs,dstran) -
     + sigma_t*dev - matmul(ddemsd,dgamma)
!
!
!
!     Elastic strain increment (trial value)
      dstranp=0.
      do is = 1, nslip
          dstranp(:) = dstranp(:) + Schmidvec(is,:)*dgamma(is)
      end do
!
!
!     Subtract the plastic strain increment from total
      dstrane = dstran - dstranp
!
!
!
!     Increment in RSS
      dtau = matmul(ddemsd,dstrane)
!
!
!
!!     stress component checks
!      flag=0
!      do is = 1, nslip
!          if (dabs(dtau(is)).gt.tauc_t(is)) then
!              flag=1
!              dtau(is) = sign(1.,dtau(is))*tauc_t(is)
!          endif
!      enddo      
!!
!!
!!
!      if (flag==1) then
!!         Solve for stress increments
!          Amat(1:nslip,1:6) = transpose(ddemsd)
!          Amat(nslip+1,1:6) = [1., 1., 1., 0., 0., 0.]
!          call SVDgeninverse(Amat,nslip+1,6,invAmat,err)
!!
!!
!!         Check for the inversion
!          if(err == 1) then
!!             did not converge
!              cpconv = 0
!              return
!          endif
!!      
!!         vectorize the given
!          bvec(1:nslip) = dtau
!          bvec(nslip+1) = 0.
!          dsigma = matmul(invAmat,bvec)
!!
!!
!      end if

!
!
!
!     update stress
      sigma = sigma_t + dsigma      
!
!
!
!
!
!     update RSS
      tau = tau_t + dtau
!
!
!
!     Update slip resistance
      dtauc = matmul(ddtauc_ddgamma,dgamma)
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
!     Plastic strain-related quantities for hardening calculations
!
!!     This does not make any change, commented out!
!      gammadot = dgamma/dt
!
      Lp=0.; Dp=0.;
      do is=1, nslip
!         Update plastic velocity gradient         
          Lp = Lp + gammadot(is)*Schmid_0(is,:,:)
!
!         Update plastic strain rate
          Dp = Dp + gammadot(is)*Schmid(is,:,:)
      end do
      Dp = (transpose(Dp)+Dp)/2.
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
!     Total slip over time per slip system
      gammasum = gammasum_t +! dgamma
     + gammadot*dt !dgamma
!
!
!
!     Total slip
      totgammasum = totgammasum_t +! sum(abs(dgamma))
     + sum(abs(gammadot*dt)) !sum(abs(dgamma))
!
!
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
      domegap33=0.
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
          domegap33 = domegap33 + W33*dgamma(is)
!
!
      end do
!
!
!     Elastic spin
      domegae33 = domega33 - domegap33
!
!
!
!
!
!
!
!
!     Elastic strains in the crystal reference
!
!
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
!     Orientation update
      dgmatinv = matmul(dstrane33+domegae33,gmatinv_t)
!
!
      gmatinv = gmatinv_t + dgmatinv
!
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
      end subroutine CP_Huang_linear
!
!
!
!  
!
!      
!
!     Fully implicit state update
      subroutine CP_Huang_fullyimplicit(matid, phaid,
     + nslip, nscrew, mattemp, Cs, gf, G12,
     + burgerv, cubicslip, caratio,
     + tau_t, Fp_t, gmatinv_t, Eec_t,
     + gammadot_t, gammasum_t,
     + totgammasum_t, evmp_t, sigma_t,
     + forestproj, slip2screw,
     + dirc_0, norc_0,
     + Schmidvec, Schmid_0,
     + Schmid, SchmidxSchmid,
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
     + vecmat9, matvec9
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
!     RSS at former time step
      real(8), intent(in) :: tau_t(nslip)
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
      real(8), intent(in) :: dirc_0(nslip,3)
!     deformed slip plane normal in sample reference frame
      real(8), intent(in) :: norc_0(nslip,3)
!     Schmid tensor - vectorized
      real(8), intent(in) :: Schmidvec(nslip,6)
!     Initial Schmid tensor - unused
      real(8), intent(in) :: Schmid_0(nslip,3,3)  
!     Schmid tensor - unused
      real(8), intent(in) :: Schmid(nslip,3,3)
!     Schmid dyadic - unused
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
      real(8) :: Pmat_s(6,6), Pmat_c(6,6)
      real(8) :: gammadot_s(nslip), gammadot_c(nslip)    
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
!     Variables used in the solution for slip increments
      real(8) :: Nab(nslip,nslip)
      real(8) :: Mab(nslip,nslip)
!     The derivative of shear rates with respecto to rss
      real(8) :: ddgdde(nslip,6)
!     Slip increments
      real(8) :: dgamma(nslip)
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
      real(8) :: dsigma(6)
!
!
!     Determinant of plastic deformation gradient
      real(8) :: detFp
!     Strain increment and related variables
      real(8) :: dstranp(6), dstrane(6)
      real(8) :: dstrane33(3,3), dstranp33(3,3)
!     crss increment
      real(8) :: dtauc(nslip)
!     ssd density increment
      real(8) :: dssd(nslip)
!     loop density increment
      real(8) :: dloop(maxnloop)
!     backstress increment
      real(8) :: dX(nslip)
!     total ssd density increment
      real(8) :: dssdtot
!     forest dislocation density increment
      real(8) :: dforest(nslip)
!     substructure dislocation density increment
      real(8) :: dsubstructure
!
!
!
!     rss
      real(8) :: tau(nslip), dtau(nslip)
!
!     overall crss
      real(8) :: tauceff(nslip)
!     forest density
      real(8) :: rhofor(nslip)
!     total density
      real(8) :: rhotot(nslip)
!     total scalar density
      real(8) :: sumrhotot
!
!     Added variables for the implicit solver
!     deformed slip direction in sample reference frame at current time step
      real(8) :: ns(3)
!     deformed slip plane normal in sample reference frame at current time step
      real(8) :: nn(3)
!     Slip deformation tensor
      real(8) :: P(6,nslip)
!     Schmid tensor (3x3) and vector (dummy)
!      real(8) :: SNij(3,3), NSij(3,3), sni(6), nsi(6)
!     Hardening derivative
      real(8) :: ddtauc_ddgamma(nslip,nslip)
      real(8) :: ddtauc_ddgamma_dgamma(nslip)
!     Old increments
      real(8) :: dsigma_0(6), dgamma_0(nslip)
      real(8) :: dtau_0(nslip), dtauc_0(nslip)
      real(8) :: dgmatinv_0(3,3)
!     Linear solution
      real(8) :: dsigma_1(6), dgamma_1(nslip)
      real(8) :: dtau_1(nslip), dtauc_1(nslip)
      real(8) :: dgmatinv_1(3,3), gammadot_1(nslip)
      real(8) :: Fp_1(3,3), Eec_1(6), evmp_1
      real(8) :: tausolute_1, dssdtot_1
      real(8) :: dssd_1(nslip), dloop_1(maxnloop)
      real(8) :: dX_1(nslip), dforest_1(nslip)
      real(8) :: dsubstructure_1, jacobi_1(6,6)
!     
!     residual
      real(8) :: psi(nslip), psinorm
!
!     Dummy variables
      real(8) :: dummy33(3,3), dummy33_(3,3)
      real(8) :: dummy0, dummy6(6), dummy66(6,6)
!
!     Variables uses within the subroutine
      integer :: is, js, il, i, j, k, a, b
      integer :: nitrtn
      integer :: debug, debugwait
!     
!     Turn VS debugger on/off
      debug=0
      do while (debug==1)
          debugwait = 1
      end do
!
!     Convergence flag
      cpconv = 1
!
!     Orientation matrix
      gmatinv = gmatinv_t
!
!
!     Stress
      sigma = sigma_t
!     States (initial values)
      tauc = tauc_t
      X = X_t
      rhofor = rhofor_t
      tauceff = tauceff_t
      rhotot = rhotot_t
      sumrhotot = sumrhotot_t
      ssdtot = ssdtot_t
      forest = forest_t
      substructure = substructure_t
      ssd = ssd_t
      loop = loop_t
!     Assign initial RSS
      tau = tau_t
      
!     Total slip per slip system
      gammasum = gammasum_t
!
!     Set the old increments to zero
      dsigma_0=0.; dgamma_0=0.
      dtau_0=0.; dtauc_0=0.
      dgmatinv_0=0.
!
!
!
!
!     Total spin increment (3x3 matrix)
      domega33 = 0.
      domega33(1,2) =  domega(1)
      domega33(1,3) = -domega(2)
      domega33(2,1) = -domega(1)
      domega33(2,3) =  domega(3)
      domega33(3,1) =  domega(2)
      domega33(3,2) = -domega(3)
!
!
!     Volumetric change in the strain
      dev = dstran(1) + dstran(2) + dstran(3)
!
!     
!
!
!
!
!
!     Number of iterations
      nitrtn = 0
!     Initial residue
      psinorm=1.
!
!     Iteration starts here!
      do while ((psinorm > tolerance).and.(nitrtn < maxniter))
          
!         Plastic spin dyadic
          W = 0.
!         Plastic deformatin tensor
          P = 0.
          do is = 1, nslip
!
!             Slip direction
              ns = matmul(gmatinv,dirc_0(is,:))
!             Slip plane normal
              nn = matmul(gmatinv,norc_0(is,:))
!
!             Deformation tensor
              P(1,is) = ns(1)*nn(1)
              P(2,is) = ns(2)*nn(2)
              P(3,is) = ns(3)*nn(3)
              P(4,is) = ns(2)*nn(1)+ns(1)*nn(2)
              P(5,is) = ns(1)*nn(3)+ns(3)*nn(1)
              P(6,is) = ns(3)*nn(2)+ns(2)*nn(3)
!
!             Spin tensor
              W(1,is) = 0.5*(ns(1)*nn(2)-ns(2)*nn(1))
              W(2,is) = 0.5*(ns(3)*nn(1)-ns(1)*nn(3))
              W(3,is) = 0.5*(ns(2)*nn(3)-ns(3)*nn(2))
!
!
          enddo
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
              beta(1,is)=-2.*W(2,is)*sigma(5)+2.*W(1,is)*sigma(4)
              beta(2,is)=-2.*W(1,is)*sigma(4)+2.*W(3,is)*sigma(6)
              beta(3,is)=-2.*W(3,is)*sigma(6)+2.*W(2,is)*sigma(5)
              beta(4,is)=-W(1,is)*sigma(1) + W(1,is)*sigma(2) - 
     + W(2,is)*sigma(6) + W(3,is)*sigma(5)
              beta(5,is)=-W(2,is)*sigma(3) + W(2,is)*sigma(1) + 
     + W(1,is)*sigma(6) - W(3,is)*sigma(4)
              beta(6,is)=-W(3,is)*sigma(2) - W(1,is)*sigma(5) + 
     + W(2,is)*sigma(4) + W(3,is)*sigma(3)
!
              ddemsd(:,is) = matmul(Cs,P(:,is)) + beta(:,is)
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
     + Lp_s,Dp_s,Pmat_s,gammadot_s,
     + dgammadot_dtau_s,dgammadot_dtauc_s) 
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
              dgammadot_dtau_c = 0.
              dgammadot_dtauc_c = 0.
!
          elseif (creepmodel == 1) then
!
!
!
              call expcreep(Schmid_0,Schmid,SchmidxSchmid,
     + tau,X,tauceff,dt,nslip,phaid,
     + mattemp,creepparam,gammasum,
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
              return
          endif      
!          
!
!         Check infinity
          if(any(gammadot*0./=gammadot*0.)) then
!             did not converge
              cpconv = 0
              return
          endif
!
!
!         calculate von mises invariant plastic strain rate
          pdot=sqrt(2./3.*sum(Dp*Dp))
!
!
!
!
!
!         Total slip
          totgammasum = totgammasum_t +
     + sum(abs(gammadot))*dt
!
!
!
!
          X = X_t; dX = 0.
!         Update backstress
          if (backstressmodel==1) then
!
              call backstressmodel1(backstressparam,
     + nslip,X_t,gammadot,dt,dX)
!
!             Update the backstress
              X = X_t + dX
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
     + tauc,ssd,loop,
     + rhofor,rhotot_t,substructure,
     + tausolute,dtauc,dssdtot,dforest,
     + dsubstructure,dssd,dloop,
     + ddtauc_ddgamma)
!
!
!
!
!
!
!
!         Solution
          Nab = 0.
          dgamma = 0.
          do a = 1, nslip
!
              do b = 1, nslip
!
                  dummy0 = 0.
                  do i = 1, 6
                      dummy0 = dummy0 + ddemsd(i,a)*P(i,b)
                  end do
!
                  Nab(a,b) = theta*dt*dgammadot_dtau(a)*dummy0 -
     + ddtauc_ddgamma(a,b)*dgammadot_dtauc(a)*theta*dt
!     + *sign(1.0,gammadot(b))
!
                  if (nitrtn>0) then
!
                      Nab(a,b) = Nab(a,b) -
     + dgammadot_dtauc(a)*ddtauc_ddgamma_dgamma(b)*theta*dt 
!
                  end if
!
!
              end do
!
              Nab(a,a) = Nab(a,a) + 1.
!
              if (nitrtn==0) then
                  
!                 Given quantities in vector form
                  dgamma(a) = dot_product(ddemsd(:,a),dstran)*
     + dgammadot_dtau(a)*theta*dt + gammadot(a)*dt
!
              else

                  dgamma(a) = theta*dt*(gammadot(a)-gammadot_1(a)) +
     + dt*gammadot_1(a) - dgamma_0(a)

              end if
!
          end do
!
!
!         Solve for the slip increments
          call nolapinverse(Nab,Mab,nslip)
!              
!
!         Check for the inversion
          if(any(Mab /= Mab)) then
!             did not converge
              cpconv = 0
              return
          endif
!
!
!         Slip increments
          dgamma = matmul(Mab,dgamma)          
!
!         Update increment
          dgamma = dgamma + dgamma_0
!
!         Update total slip
          gammasum = gammasum - dgamma_0 + dgamma
!
!
!
!
!         Total slip
          totgammasum = totgammasum_t +
     + sum(abs(gammadot))*dt     
!
!
!         
!
!
!!         Plastic velocity gradient
!          Lp = 0.
!          do is = 1, nslip
!              Lp = Lp + gammadot(is)*Schmid_0(is,:,:)
!          end do
!!
!!
!!
!!
!!         variables for plastic part of the deformation gradient
!          !dummy33 = I3 - Lp*dt
!          !call inv3x3(dummy33,dummy33_,dummy0)
!          dummy33_ = I3 + Lp*dt
!!
!!         plastic part of the deformation gradient
!          Fp = matmul(dummy33_,Fp_t)
!!
!!         determinant
!          call deter3x3(Fp,detFp)
!!
!!
!!
!!
!!         check wheter the determinant is negative
!!         or close zero
!          if (detFp <= smallnum) then
!!             did not converge
!              cpconv = 0
!              return
!          else
!!             Scale Fp with its determinant to make it isochoric
!              Fp = Fp / detFp**(1./3.)
!!
!          end if     
!
!
!         calculate von mises invariant plastic strain rate
          pdot=sqrt(2./3.*sum(Dp*Dp))
!
!
!!         Von-Mises equivalent total plastic strain
!          evmp = evmp_t + pdot*dt
!
!
!
!
!         Update the states using hardening laws
          call hardeningrules(phaid,nslip,
     + mattemp,dt,G12,burgerv,
     + totgammasum,dgamma/dt,pdot,
     + irradiationmodel,irradiationparam,
     + hardeningmodel,hardeningparam,
     + hintmat1,hintmat2,
     + tauc,ssd,loop,
     + rhofor,rhotot_t,substructure,
     + tausolute,dtauc,dssdtot,dforest,
     + dsubstructure,dssd,dloop,
     + ddtauc_ddgamma)
!
!
          dtauc = matmul(ddtauc_ddgamma,dgamma)
!              
!
!
!         Calculate plastic spin from slip
          domega33_p=0.
          do is = 1, nslip
!
              W33 = 0.
              W33(1,2) =  W(1,is)
              W33(1,3) = -W(2,is)
              W33(2,1) = -W(1,is)
              W33(2,3) =  W(3,is)
              W33(3,1) =  W(2,is)
              W33(3,2) = -W(3,is)        
!
              domega33_p = domega33_p + W33*dgamma(is)
!
!
          end do
!
!
!         Elastic spin
          domega33_e = domega33 - domega33_p
!
!
!
!
!
!         Orientation update
          dgmatinv = matmul(domega33_e,gmatinv_t)
!
!
          gmatinv = gmatinv_t + dgmatinv
!
!
!
!
!         Elastic strains in the crystal reference
!
!
!         Elastic strain increment
          dstranp=0.
          do is = 1, nslip
              dstranp = dstranp + P(:,is)*dgamma(is)
          end do
!
!
!         Subtract the plastic strain increment from total
          dstrane = dstran - dstranp
!
!
!         increment in RSS
          dtau = matmul(transpose(ddemsd), dstrane)
!
!         undo shear corrections
          dstrane(4:6) = 0.5*dstrane(4:6)
!
!         Convert to 3x3 matrix
          call vecmat6(dstrane,dstrane33)
!
!         
!         Elastic strains in the crystal reference
          dummy33_ = matmul(transpose(gmatinv),dstrane33)
          dummy33 = matmul(dummy33_,gmatinv)
!
!         Vectorize
          call matvec6(dummy33,dummy6)
!
!         Shear corrections
          dummy6(4:6) = 2.0*dummy6(4:6)
!
!         Add the strain increment to the former value
          Eec=Eec_t + dummy6
!
!
!
!
!         Stress update
          dsigma = matmul(Cs,dstran)-sigma*dev-matmul(ddemsd,dgamma)
!
          sigma = sigma + dsigma - dsigma_0
!
!
!
!         Material tangent
!
!         Step-1: Calculate ddgamma_ddeps
          do is=1,nslip
!
              ddgdde(is,:) = dt*theta*dgammadot_dtau(is)*ddemsd(:,is)
!
          end do
!
!
!
!
!         Step-2: Calculate overall expression
          jacobi = Cs - matmul(ddemsd,matmul(Mab,ddgdde))
!
!
!         Correction for large deformations
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
!         Make it symmetric to help convergence
          jacobi = 0.5*(jacobi + transpose(jacobi))
!
!
!         Store the linear iteration results
          if (nitrtn==0) then
              dsigma_1 = dsigma
              dtau_1 = dtau
              dtauc_1 = dtauc
              dgamma_1 = dgamma
              gammadot_1 = gammadot
              dgmatinv_1 = dgmatinv
              Eec_1 = Eec
              Fp_1 = Fp
              evmp_1 = evmp
              tausolute_1 = tausolute
              dssdtot_1 = dssdtot
              dssd_1 = dssd
              dloop_1 = dloop
              dX_1 = dX
              dforest_1 = dforest
              dsubstructure_1 = dsubstructure
              jacobi_1 = jacobi
          end if
!
!         State update
          tau = tau + dtau - dtau_0
          tauc = tauc + dtauc - dtauc_0
!
!         Store the states
          dsigma_0 = dsigma
          dtau_0 = dtau
          dtauc_0 = dtauc
          dgamma_0 = dgamma
          dgmatinv_0 = dgmatinv         
!    
!
!         Update the other states
!
!
          ssd = ssd_t + dssd
!
          ssdtot = ssdtot_t + dssdtot
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
!         Find the effective hardening increment
!
!         Calculate total and forest density
          call totalandforest(phaid,
     + nscrew, nslip, gnd_t,
     + ssd, ssdtot, forest,
     + forestproj, slip2screw,
     + rhotot, sumrhotot, rhofor)
!
!
!
!         Calculate crss
          call slipresistance(phaid, nslip, gf, G12,
     + burgerv, sintmat1, sintmat2,
     + tauc, rhotot, sumrhotot, rhofor,
     + substructure, tausolute, loop,
     + hardeningmodel, hardeningparam,
     + irradiationmodel, irradiationparam,
     + mattemp, tauceff)
!          
!
!     Residual calculation
!
!     Enter the slip and creep laws
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
     + Lp_s,Dp_s,Pmat_s,gammadot_s,
     + dgammadot_dtau_s,dgammadot_dtauc_s) 
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
              dgammadot_dtau_c = 0.
              dgammadot_dtauc_c = 0.
!
          elseif (creepmodel == 1) then
!
!
!
              call expcreep(Schmid_0,Schmid,SchmidxSchmid,
     + tau,X,tauceff,dt,nslip,phaid,
     + mattemp,creepparam,gammasum,
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
!
!
!
!         Residual
          psi = theta*dt*gammadot + dt*(1.-theta)*gammadot_1 - dgamma 
!
!         Norm of the residual
          psinorm = sqrt(sum(psi*psi))
!
!         Increment iteration
          nitrtn = nitrtn + 1
!
!         Compute derivative x dgamma
          ddtauc_ddgamma_dgamma = matmul(ddtauc_ddgamma,dgamma)
!
!
      end do
!     End of iterations
!
!
!
!
!
!     If iterations do not converge
      if ((nitrtn==maxniter).or.(cpconv==0)) then
!
          cpconv = 1
          sigma = sigma_t + dsigma_1
          gammadot = dgamma_1/dt
          gammasum = gammasum_t + dgamma_1
          totgammasum = totgammasum_t +
     + sum(abs(dgamma_1))
          Fp = Fp_1
          Eec = Eec_1
          evmp = evmp_1
          tau = tau_t + dtau_1
          tauc = tauc_t + dtauc_1
          gmatinv = gmatinv_t + dgmatinv_1
          ssdtot = ssdtot_t + dssdtot_1
          ssd = ssd_t + dssd_1
          forest = forest_t + dforest_1
          substructure = substructure_t + dsubstructure_1
          tausolute = tausolute_1
          loop = loop_t + dloop_1
          X = X_t + dX_1
          jacobi = jacobi_1
!
      end if
!
!
! REDO PLASTIC CALCULATIONS
!******************************************************
!         Plastic velocity gradient
          Lp = 0.
          do is = 1, nslip
              Lp = Lp + gammadot(is)*Schmid_0(is,:,:)
          end do
!
!
!
!
!
!         Total slip
          totgammasum = totgammasum_t +
     + sum(abs(gammadot))*dt     
!
!
!
!
!         variables for plastic part of the deformation gradient
          !dummy33 = I3 - Lp*dt
          !call inv3x3(dummy33,dummy33_,dummy0)
          dummy33_ = I3 + Lp*dt
!
!         plastic part of the deformation gradient
          Fp = matmul(dummy33_,Fp_t)
!
!         determinant
          call deter3x3(Fp,detFp)
!
!
!
!
!         check wheter the determinant is negative
!         or close zero
          if (detFp <= smallnum) then
!             did not converge
              cpconv = 0
              return
          else
!             Scale Fp with its determinant to make it isochoric
              Fp = Fp / detFp**(1./3.)
!
          end if     
!
!
!         calculate von mises invariant plastic strain rate
          pdot=sqrt(2./3.*sum(Dp*Dp))
!
!
!         Von-Mises equivalent total plastic strain
          evmp = evmp_t + pdot*dt
!
!
!******************************************************
!
!
!
!      
      return
      end subroutine CP_Huang_fullyimplicit
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
!     Calculates the total and forest density
!     from SSD and GND densities using summation
!     and forest projections
      subroutine totalandforest(iphase, nscrew, nslip,
     + gnd, ssd, ssdtot, forest, 
     + forestproj, slip2screw,
     + rhotot, sumrhotot, rhofor)
      use utilities, only : vecprod
      implicit none
      integer, intent(in) :: iphase
      integer, intent(in) :: nscrew
      integer, intent(in) :: nslip
      real(8), intent(in) :: gnd(nslip+nscrew)
      real(8), intent(in) :: ssd(nslip)
      real(8), intent(in) :: ssdtot
      real(8), intent(in) :: forest(nslip)
      real(8), intent(in) :: forestproj(nslip,nslip+nscrew)
      real(8), intent(in) :: slip2screw(nscrew,nslip)
      real(8), intent(out) :: rhotot(nslip)
      real(8), intent(out) :: sumrhotot
      real(8), intent(out) :: rhofor(nslip)
!
!     local variables
!
      real(8) vec(3)
      integer i, j
!
!
!
!
!
!     Compute forest density
!     Effect of screws are ignored!!!
!     Add gnd and sdd contributions (if exist)
!     Forest projections
      rhofor = forest + matmul(forestproj, abs(gnd)) +
     + matmul(forestproj(1:nslip,1:nslip), ssd) + ! edges
     + matmul(forestproj(1:nslip,nslip+1:nslip+nscrew),
     + matmul(slip2screw,ssd)) ! screws
!
!

      rhotot = ssd +
     + abs(gnd(1:nslip)) + ! edges
     + matmul(transpose(slip2screw), 
     + abs(gnd(nslip+1:nslip+nscrew))) ! screws
!
!
!   
!     Scalar sum
      sumrhotot = sqrt(sum(gnd*gnd)) + ssdtot
!
!
      return
      end subroutine totalandforest
!
!
!
!
      end module cpsolver