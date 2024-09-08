! 
! November, 01st, 2022 - 1st working version
! May, 20th, 2024 - Release v2.20
!
! Crystal Plasticity UMAT
! University of Oxford
! Eralp Demir
!
! Sponsored by Design-by-Fundamentals (DbF) project under STEP program of UKAEA
!
! Rewrite of UMAT by Ed Tarleton and Nicolo Grilli
! Originally based on the UEL by Fionn Dunne 2007
! Deformation twinning is not included
!
! Major changes:
!     - Initialization routines for computational efficiency
!     - Reverse slip formulation by C. Hardie
!     - Option for implicit state update
!     - Multiple materials with different phases can co-exist in a mesh
!     - Modular code for flexible constitutive model development
!     - GND calculations:
!         o Restricted solution to the active slip systems
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
      include "crss.f"
      include "initializations.f"
      include "slip.f"
      include "slipreverse.f"
      include "creep.f"
      include "innerloop.f"
      include "hardening.f"
      include "backstress.f"
      include "cpsolver.f"
      include "straingradients.f"
!
!
!
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
!     Subroutine for initialization 
      use initializations, only : initialize_variables,
     + initialize_once 
      use userinputs, only : gndmodel, backstressmodel
      use globalvariables, only: numel, numpt, 
     + init_once, statev_gmatinv, statev_gmatinv_t,
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
     + statev_tausolute_t, statev_tausolute, time_old, dt_t,
     + statev_backstress, statev_backstress_t,
     + statev_plasdiss, statev_plasdiss_t, 
     + ip_count, calculategradient, ip_init, grad_init
!
      use straingradients, only: gndmodel1, gndmodel2, 
     + gndmodel3, gndmodel4
      use backstress, only: backstressmodel2
      use meshprop, only: initialize_gradientoperators
      use useroutputs, only: write_statev_legend   
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
      integer :: i
!
!
!
!     at the start of the analysis (only ONCE!)
!     if there are initializations/calculations that are needed once,
!     and that are independepent of element properties, you may use here.
      if (LOP==0) then
!
!         0. Initialize variables (instead of using data statements)
          call initialize_variables
!
!
      end if
!
!     Initialization of gradient operators
!     Check if the element has more than one integration point
!     And the gradient initialization was done before or not
      if ((calculategradient==1).and.(grad_init==0)) then
!
!         Check if IP coordinates were collected
          if ((sum(ip_init)==numel*numpt).and.
     + (numel>0).and.(numpt>0)) then
!
!
!             Calculate the gradient mapping for each element once.
              call initialize_gradientoperators
!
              grad_init=1
              write(*,*) '7. Gradient operators initialized!'
!
          endif
!
      end if
!
!
!
!
!
      if (LOP==1) then
!
!
!         in case of force BC
          if ((KINC==1).and.(KSTEP==1)) then
!
!             
!             
!             check if the one-time initialization is done or not
              if ((init_once==0).and.(numel>0).and.(numpt>0)) then
!            
!
!
!                 check the number of elements
                  i = sum(ip_count(1,:))
                  if (i>numpt) numpt = i
!
!                 check the number of elements
                  i = sum(ip_count(:,1))
                  if (i>numel) numel = i
!
!                 
!                 write element number
                  write(*,*) 'total elements in the mesh: ', numel
!
!                 write integration point number
                  write(*,*) 'integration points per element: ', numpt
!
!
                  write(*,*) '--------------------------------'
!
!                 call the initializations
                  call initialize_once
!
!                 display upon completion
                  write(*,*) 'one-time initialization is complete!'
!
!
                  write(*,*) '--------------------------------'
!
!
!                 write the legend to a text file
                  call write_statev_legend
!
!                 message for initializatoin
                  write(*,*) '6. "STATEV_legend.txt" file is ready!'
!
!                 set the one-time initilaziation flag (at the very end)
                  init_once=1
!
              end if
!
!
!
!
!
          end if
!
!
      end if
!
!
!
!     at the end of each increment
      if (LOP==2) then
!
!
!         in case of displacement BC
          if ((KINC==1).and.(KSTEP==1)) then
!
!             
!             
!             check if the one-time initialization is done or not
              if ((init_once==0).and.(numel>0).and.(numpt>0)) then
!            
!
!
!                 check the number of elements
                  i = sum(ip_count(1,:))
                  if (i>numpt) numpt = i
!
!                 check the number of elements
                  i = sum(ip_count(:,1))
                  if (i>numel) numel = i
!
!                 
!                 write element number
                  write(*,*) 'total elements in the mesh: ', numel
!
!                 write integration point number
                  write(*,*) 'integration points per element: ', numpt
!
!
                  write(*,*) '--------------------------------'
!
!                 call the initializations
                  call initialize_once
!
!                 display upon completion
                  write(*,*) 'one-time initialization is complete!'
!
!
                  write(*,*) '--------------------------------'
!
!
!                 write the legend to a text file
                  call write_statev_legend
!
!                 message for initializatoin
                  write(*,*) '6. "STATEV_legend.txt" file is ready!'
!
!
!
!                 set the one-time initilaziation flag (at the very end)
                  init_once=1
!
!             
              end if
!
!
!
!
!
!
          end if

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
!
!         GND calculations
!         do this only if the initiliazation complete and 
!         element gradients are calculatable (calculategradient==1)
          if ((init_once==1).and.(grad_init==1).and.
     + (calculategradient==1).and.(KINC>1)) then
!
!             update and calculate GNDs (nonlocal calculations)
!             this is done at the end of calculations.
!             GNDs that belong to the PREVIOUS time step are used.
!             initially GNDs are assumed to have "0" values.
!
!             calculate GNDs (if initialization complete)
              if (gndmodel>0) then
!
!                 curl followed by L2 minimization - SVD -
!                 constrained to active slip systems
                  if (gndmodel==1) then
!
                      call gndmodel1
!
!                 Cermelli-Gurtin's incompatibility measure - L2 minimization - SVD
!                 constrained to active slip systems (same as model-1)
                  elseif (gndmodel==2) then
!
                      call gndmodel2
!
!                 rate form followed by direct projections
                  elseif (gndmodel==3) then
!
                      call gndmodel3(DTIME(1))
!
!                 slip gradients
                  elseif (gndmodel==4) then
!
                      call gndmodel4(DTIME(1))
!
!
                  endif
!
!
                  write(*,*) 'GNDs are computed!'
!
!
                  if (backstressmodel==2) then
!
                      call backstressmodel2
!
                      write(*,*) 'Backstresses are computed!'
!
                  end if
!
              end if
!
          end if
!
!         update the time
          time_old = time(2)
!         store former converged time increment
          dt_t = DTIME(1)
!
!
!         update state variables
!         assign the current results to the former values
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
          statev_backstress_t(:,:,:) = statev_backstress(:,:,:)
          statev_plasdiss_t(:,:) = statev_plasdiss(:,:)
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
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
!
      use userinputs, only: cutback, pastefront
      use globalvariables, only: numel, numpt, numtens,
     + largenum, ip_init, init_once, ip_count
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
     + KINC       ! Increment number (time increment)
      INTEGER, DIMENSION(4),          INTENT(IN) ::
     + JSTEP     ! Step number (load step)
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
!     local array for 3D crystal plasticity solver
      real(8) :: sigma(6), jacobi(6,6)
!     material-id needed for array allocations
      integer :: matid, debug, debugwait
!     counter
      integer :: i
!
!     turn VS debugger on/off
      debug=0
      do while (debug==1)
          debugwait = 1
      end do
!
!     to avoid warning messages set the following to zero
      DDSDDT = 0.
      DRPLDE = 0.
!
!     outputs
      sigma = 0.
      jacobi = 0.
!
!
!
!     read the number of elements and number of integration points per element
      if (ip_count(NOEL,NPT)<1) then
!

!         Increment the counter
          ip_count(NOEL,NPT)=ip_count(NOEL,NPT)+1
!
!         store the number of elements
          if (NOEL>numel) numel=NOEL
!         store the number of integration points
          if (NPT>numpt) numpt=NPT
!
!         dimension of the problem (will be re-assigned in feprop)
          numtens = ntens
!
          DDSDDE=0.
          do i=1,NTENS
              DDSDDE(i,i)=largenum   
          end do
          STRESS=0.
          PNEWDT=cutback
!
!
!
      end if
!
!
!     if arrays allocated then
!     do the crystal plasticity calculations
      if (init_once==1) then
!
!
!         material/phase identifier
          matid=int(PROPS(5))
!
!         in some multi-body simulations UMAT entry for RGB has happened!
!         in that case, RGB having no material-ID assigned gave array access issues
!         this case deals with that situation
          if (matid > 0) then
!
!             material property initialization
              if (ip_init(NOEL,NPT)==0) then
!
                  call initialize_atfirstinc(NOEL,NPT,COORDS,
     + NPROPS,PROPS,TEMP,NSTATV)
!
!
!                 write(*,*) 'element: ', NOEL
!                 write(*,*) 'ip: ', NPT
!                 write(*,*) 'initialized!'
              end if
!
!
!
!
!             call the main crystal plasticity solver
              call solve(NOEL,NPT,DFGRD1,DFGRD0,
     + TEMP,DTEMP,DTIME,matid,
     + PNEWDT,NSTATV,STATEV,
     + sigma,jacobi)
!
!
!
!             increase the time step if desired or,
!             let ABAQUS decide!
              if (pastefront > 1.) then
!             if there is no cutback introduced
                  if (PNEWDT /= cutback) then
                      PNEWDT = pastefront
                  end if
              end if
!
!
!
!             3D case
              if (NTENS==6) then
!
                  STRESS = sigma
                  DDSDDE = jacobi
!
!             plane strain case
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
!             plane stress case
              elseif (NTENS==3) then
!
                  STRESS=0.
!                 correction by Alvaro
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
              end if
!
!         if rigid body (or matid not defined in PROPS)
          else
!
              DDSDDE=0.
              do i=1,NTENS
                  DDSDDE(i,i)=largenum    
              end do
              STRESS=0.              
!
!         end of crystal plasticity solution
          end if
!
!     end of one-time initialization performed case
      end if
!
!
!
      RETURN
      END