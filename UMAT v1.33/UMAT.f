! 
! Nov., 1st, 2022 - 1st working version
! May 25th, 2023 - Release v1.26
!
! Crystal Plasticity UMAT for DBF-project
! University of Oxford
! Eralp Demir
!
! Rewrite of UMAT by Ed Tarleton and Nicolo Grilli
! which was based on UEL by Fionn Dunne 2007
!
! Major changes:
!     - Initialization routines for computational efficiency
!     - Forward-gradient Euler predictor scheme for better convergence
!     - Option for implicit state update
!     - Multiple materials with different phases can co-exist in a mesh
!     - Modular code for flexible constitutive model development
!     - GND calculations:
!         o Multiple options for GND calculation (i.e. curlFp or slip gradients)
!         o Option for element center homogenization
!         o Different 2D/3D element types for GND calculation
!     - 2D plane stress/strain problems
!
      include "userinputs.f"
      include "globalvariables.f"
      include "irradiation.f"
      include "errors.f"
      include "utilities.f"
      include "meshprop.f"
      include "usermaterials.f"
      include "useroutputs.f"
      include "initializations.f"
      include "crss.f"
      include "slip.f"
      include "creep.f"
      include "hardening.f"
      include "cpsolver.f"
      include "straingradients.f"
      include "backstresses.f"
!
!
!
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
!     Subroutine for initialization 
      use initializations, only : initialize_once
      use userinputs, only : numel, gndmodel, backstressmodel
      use globalvariables, only: 
     + gnd_init, statev_gmatinv, statev_gmatinv_t,
     + statev_gammasum_t, statev_gammasum,
     + statev_gammadot_t, statev_gammadot,
     + statev_Fp_t, statev_Fp, statev_sigma_t, statev_sigma,
     + statev_jacobi_t, statev_jacobi, statev_Fth_t, statev_Fth,
     + statev_tauc_t, statev_tauc, statev_maxx_t, statev_maxx,
     + statev_Eec_t, statev_Eec, statev_gnd_t, statev_gnd,
     + statev_ssd_t, statev_ssd, statev_forest_t, statev_forest,
     + statev_substructure_t, statev_substructure, statev_evmp_t,
     + statev_evmp, statev_totgammasum_t, statev_totgammasum,
     + statev_ssdtot_t, statev_ssdtot, statev_loop_t, statev_loop,
     + statev_Lambda_t, statev_Lambda, statev_sigma_t2,
     + statev_tausolute_t, statev_tausolute, time_old, dt_t
      use meshprop, only: initialize_gradientoperators
      use useroutputs, only: write_statev_legend
      use straingradients, only: gndmodel1, gndmodel2, 
     + gndmodel3, gndmodel4
      use backstresses, only: backstressmodel1
!
      implicit none
!
!
!
      INTEGER,                        INTENT(IN) ::
     + LOP,
     + LRESTART,
     + KSTEP,
     + KINC
      REAL(8), DIMENSION(1),          INTENT(IN) ::
     + DTIME
      REAL(8), DIMENSION(2),          INTENT(IN) ::
     + TIME
!
!
!
!
!
!     At the start of the analysis (only ONCE!)
!     If there are initializations/calculations that are needed once,
!     and that are independepent of element properties, you may use here.
      if (LOP==0) then
!
          write(*,*) '--------------------------------'
!
!         Call the initializations
          call initialize_once
!
!         Display upon completion
          write(*,*) 'initialization once at the
     + beginning is complete!'
!
!
          write(*,*) '--------------------------------'
!
      endif
!
!
!
!
!
!     At the end of each increment
      if (LOP==2) then
!
!
!
          write(*,*) '--------------------------------'
!
          write(*,*) 'At the end of increment: ', KINC
!
!
!
!
!         End of first increment
!
          if (KINC == 1) then
!
!
!             write element number
              write(*,*) 'Total element number in the mesh: ', numel
!
!             write the legend to a text file
              call write_statev_legend
!             message for initializatoin
              write(*,*) '5. "STATEV_legend.txt" file is ready!'
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
!         GND calculations
!
!         Update and calculate GNDs (nonlocal calculations)
!         This is done at the end of calculations.
!         GNDs that belong to the PREVIOUS time step are used.
!         Initially GNDs are assumed to have "0" values.
!         If there is a GND model defined for this case
          if (gndmodel/=0) then
!
!             Initialization of gradient operators
!             This is done only ONCE!
              if (gnd_init==0) then
!
!
!                 Calculate the gradient mapping for each element once here!
                  call initialize_gradientoperators
!
!                 Set the flag as initiliazed
                  gnd_init = 1
!
                  write(*,*) '6. Gradient operators initialized!'
!
!
              endif
!
!             Calculate GNDs (if initialization complete)
!
!             Curl followed by L2 minimization - SVD -
!             constrained to active slip systems
              if (gndmodel==1) then
!
                  call gndmodel1
!
!             Rate form followed by direct projections
              elseif (gndmodel==2) then
!
                  call gndmodel2(DTIME(1))
!
!             Slip gradients
              elseif (gndmodel==3) then
!
                  call gndmodel3(DTIME(1))
!
!             Gurtin's incompatibility measure - L2 minimization - SVD
!             constrained to active slip systems (same as model-1)
              elseif (gndmodel==4) then
!
                  call gndmodel4
!
!
              endif
!
!
              write(*,*) 'GNDs are computed!'
!
!
              if (backstressmodel==1) then
!
                  call backstressmodel1
!
                  write(*,*) 'Backstresses are computed!'
!
              end if
!
          endif
!
!
!         Update the time
          time_old = time(2)
!         Store former converged time increment
          dt_t = DTIME(1)
!
!
!         Update state variables
!         Assign the current results to the former values
          statev_gmatinv_t(:,:,:,:) = statev_gmatinv(:,:,:,:)
          statev_evmp_t(:,:) = statev_evmp(:,:)
          statev_gammasum_t(:,:,:) = statev_gammasum(:,:,:)
          statev_totgammasum_t(:,:) = statev_totgammasum(:,:)
          statev_gammadot_t(:,:,:) = statev_gammadot(:,:,:)
          statev_Fp_t(:,:,:,:) = statev_Fp(:,:,:,:)
          statev_Fth_t(:,:,:,:) = statev_Fth(:,:,:,:)
          statev_sigma_t2(:,:,:) = statev_sigma_t(:,:,:)
          statev_sigma_t(:,:,:) = statev_sigma(:,:,:)
          statev_jacobi_t(:,:,:,:) = statev_jacobi(:,:,:,:)
          statev_tauc_t(:,:,:) = statev_tauc(:,:,:)
          statev_maxx_t(:,:) = statev_maxx(:,:)
          statev_Eec_t(:,:,:) = statev_Eec(:,:,:)
          statev_gnd_t(:,:,:) = statev_gnd(:,:,:)
          statev_ssd_t(:,:,:) = statev_ssd(:,:,:)
          statev_ssdtot_t(:,:) = statev_ssdtot(:,:)
          statev_forest_t(:,:,:) = statev_forest(:,:,:)
          statev_loop_t(:,:,:) = statev_loop(:,:,:)
          statev_substructure_t(:,:) = statev_substructure(:,:)
          statev_tausolute_t(:,:) = statev_tausolute(:,:)
          statev_Lambda_t(:,:,:) = statev_Lambda(:,:,:)
!
          write(*,*) 'States are updated!'
!
          write(*,*) '--------------------------------'
!
!
      endif
!
!
!
      RETURN
      END
!
!
!
!
!
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
      use userinputs, only: cutback, pastefront
      use globalvariables, only: ip_init
      use initializations, only: initialize_atfirstinc
      use cpsolver, only: solve
      implicit none
!   
!
      INTEGER,                        INTENT(IN) ::
     + NDI,       ! Number of direct stress components at this point
     + NSHR,      ! Number of engineering shear stress components
     + NTENS,     ! Size of the stress / strain components (NDI + NSHR)
     + NSTATV,    ! Number of solution-dependent state variables
     + NPROPS,    ! User-defined number of material constants
     + NOEL,      ! Element number
     + NPT,       ! Integration point number
     + LAYER,     ! Layer number (shell elements etc.)
     + KSPT,      ! Section point within the current layer
     + KSTEP,     ! Step number (load step)
     + KINC       ! Increment number (time increment)
      CHARACTER(LEN=80),              INTENT(IN) ::
     + CMNAME     ! Usee-specified material name, left justified
      REAL(8),                        INTENT(IN) ::
     + DTIME,     ! Time increment
     + TEMP,      ! Temperature at the start of the increment
     + DTEMP,     ! Increment of temperature
     + CELENT     ! Temperature
      REAL(8), DIMENSION(1),          INTENT(IN) ::
     + PREDEF,
     + DPRED
      REAL(8), DIMENSION(2),          INTENT(IN) ::
     + TIME       ! Step time/total time at begin, of the current inc.
      REAL(8), DIMENSION(3),          INTENT(IN) ::
     + COORDS
      REAL(8), DIMENSION(NTENS),       INTENT(IN) ::
     + STRAN,     ! Total strains at beginning of the increment
     + DSTRAN     ! Strain increments
      REAL(8), DIMENSION(NPROPS),      INTENT(IN) ::
     + PROPS
      REAL(8), DIMENSION(3,3),         INTENT(IN) ::
     + DROT,      ! Rotation increment matrix
     + DFGRD0,    ! Deformation gradient at the former time increment
     + DFGRD1     ! Deformation gradient at the current time increment
      REAL(8),                         INTENT(INOUT) ::
     + PNEWDT,    ! Ratio of suggested new time increment to current time increment
     + SSE,       ! Specific elastic strain engergy
     + SPD,       ! Specific plastic dissipation
     + SCD,       ! Specific creep dissipation
     + RPL,       ! Volumetric heat generation
     + DRPLDT     ! Variation of RPL with respect to the temperature
      REAL(8), DIMENSION(NTENS),       INTENT(INOUT) ::
     + STRESS     ! Cauchy stress vector at the current time increment
      REAL(8), DIMENSION(NSTATV),      INTENT(INOUT) ::
     + STATEV     ! Solution-dependent state variables
      REAL(8), DIMENSION(NTENS),       INTENT(OUT) ::
     + DDSDDT,    ! Variation of the stress increments with respect to the temperature
     + DRPLDE     ! Variation of RPL with respect to the strain increments
      REAL(8), DIMENSION(NTENS,NTENS), INTENT(OUT) ::
     + DDSDDE     ! Material tangent at the current time increment
!
!     Local array for 3D crystal plasticity solver
      real(8) :: sigma(6), jacobi(6,6)
!     Material-id needed for array allocations
      integer :: matid, debug, debugwait
!
!     Turn VS debugger on/off
      debug=0
      do while (debug==1)
          debugwait = 1
      end do
!
!     To avoid warning messages set the following to zero
      DDSDDT = 0.
      DRPLDE = 0.
!
!
!
!
!
!     Initialize some variables at the first increment
!     If not initialized only
      if (ip_init(NOEL,NPT)==0) then
!
          call initialize_atfirstinc(NOEL,NPT,COORDS,
     + NPROPS,PROPS,TEMP,NSTATV)
!
          write(*,*) 'element: ', NOEL
          write(*,*) 'ip: ', NPT
          write(*,*) 'initialized!'
!
      endif
!
!     material/phase identifier
      matid=int(PROPS(5))
!
!
!
!     call the main crystal plasticity solver
      call solve(NOEL,NPT,DFGRD1,DFGRD0,
     + TEMP,DTEMP,DTIME,matid,
     + PNEWDT,NSTATV,STATEV,
     + sigma,jacobi)
!
!
!
!     Increase the time step if desired or,
!     let ABAQUS decide!
      if (pastefront > 1.) then
!     If there is no cutback introduced
          if (PNEWDT /= cutback) then
              PNEWDT = pastefront
          end if
      endif
!
!
!
!
!
!
!     3D case
      if (NTENS==6) then
!
          STRESS = sigma
          DDSDDE = jacobi
!
!     Plane Strain case
      elseif (NTENS==4) then
!
          STRESS=0.
          STRESS=0.
          STRESS(1) = sigma(1)
          STRESS(2) = sigma(2)
          STRESS(3) = sigma(3)
          STRESS(4) = sigma(4)
!
          DDSDDE=0.
          DDSDDE(1,1) = jacobi(1,1)
          DDSDDE(1,2) = jacobi(1,2)
          DDSDDE(1,3) = jacobi(1,3)
          DDSDDE(2,1) = jacobi(2,1)
          DDSDDE(2,2) = jacobi(2,2)
          DDSDDE(2,3) = jacobi(2,3)
          DDSDDE(3,1) = jacobi(3,1)
          DDSDDE(3,2) = jacobi(3,2)
          DDSDDE(3,3) = jacobi(3,3)
          DDSDDE(4,4) = jacobi(4,4)
!
!     Plane Stress case
      elseif (NTENS==3) then
!
          STRESS=0.
!         Correction by Alvaro
          STRESS(1) = sigma(1)-sigma(3)
          STRESS(2) = sigma(2)-sigma(3)
!
          STRESS(3) = sigma(4)
!
!
          DDSDDE=0.
          DDSDDE(1,1) = jacobi(1,1)
          DDSDDE(1,2) = jacobi(1,2)
          DDSDDE(2,1) = jacobi(2,1)
          DDSDDE(2,2) = jacobi(2,2)
          DDSDDE(3,3) = jacobi(4,4)
!
!
      endif      
!
!
!
      RETURN
      END