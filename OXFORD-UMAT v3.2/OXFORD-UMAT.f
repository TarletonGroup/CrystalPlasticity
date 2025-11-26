! OXFORD-UMAT - Crystal Plasticity Solver
! November 25th, 2025 - Release v3.2
! November 1st, 2022 - 1st working version
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
!     - Lp correction
!     - Jaumann rate correction
!
      include "userinputs.f"
      include "globalvariables.f"
      include "irradiation.f"
      include "errors.f"
      include "utilities.f"
      include "meshprop.f"
      include "usermaterials.f"
      include "miscellaneous.f"
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
      include "UEXTERNALDB.f"
!
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
     + NPROPS,PROPS,TEMP,NSTATV,NTENS,STRESS)
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
     + TEMP,DTEMP,JSTEP(1),TIME(1),DTIME,matid,
     + PNEWDT,NSTATV,STATEV,sigma,jacobi)
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
!                 Corrected stress for plane stress
                  STRESS=0.
                  STRESS(1) = sigma(1)
     + -jacobi(1,3)*sigma(3)/jacobi(3,3)
     + -jacobi(1,5)*sigma(5)/jacobi(5,5)
     + -jacobi(1,6)*sigma(6)/jacobi(6,6)
!
                  STRESS(2) = sigma(2)
     + -jacobi(2,3)*sigma(3)/jacobi(3,3)
     + -jacobi(2,5)*sigma(5)/jacobi(5,5)
     + -jacobi(2,6)*sigma(6)/jacobi(6,6)
!
                  STRESS(3) = sigma(4)
     + -jacobi(4,3)*sigma(3)/jacobi(3,3)
     + -jacobi(4,5)*sigma(5)/jacobi(5,5)
     + -jacobi(4,6)*sigma(6)/jacobi(6,6)
!
!                 Corrected jacobian for plane stress
                  DDSDDE=0.
                  DDSDDE(1,1) = jacobi(1,1) -
     + jacobi(1,3)*jacobi(3,1)/jacobi(3,3)
                  DDSDDE(1,2) = jacobi(1,2) -
     + jacobi(1,3)*jacobi(3,2)/jacobi(3,3)
                  DDSDDE(1,3) = jacobi(1,4) -
     + jacobi(1,3) * jacobi(3,4) / jacobi(3,3)
                  DDSDDE(2,1) = jacobi(2,1) -
     + jacobi(2,3)*jacobi(3,1)/jacobi(3,3)
                  DDSDDE(2,2) = jacobi(2,2) -
     + jacobi(2,3)*jacobi(3,2)/jacobi(3,3)
                  DDSDDE(2,3) = jacobi(2,4) -
     + jacobi(2,3) * jacobi(3,4) / jacobi(3,3)
                  DDSDDE(3,1) = jacobi(4,1) -
     + jacobi(4,3) * jacobi(3,1) / jacobi(3,3)
                  DDSDDE(3,2) = jacobi(4,2) -
     + jacobi(4,3) * jacobi(3,2) / jacobi(3,3)
                  DDSDDE(3,3) = jacobi(4,4) -
     + jacobi(4,3) * jacobi(3,4) / jacobi(3,3)
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
      END SUBROUTINE UMAT