!     Sept. 26th, 2022
!     Eralp Demir
!
!     User inputs
!     Global variables that the user entry is required
!     Users asked to change the model and numerical parameters here!
!
!
      module userinputs
      implicit none
!
!
!
!    FILE INPUTS
!     -------------------------------------------------------------------------------
!
!     INP file name
!     In case the file read is not successfull
      character(len=*), parameter, public :: inpfilename=
     +'/data/engs-dbf-project/engs2474/myscripts/example/Job-1.inp'
!     
!     NUMERICAL INPUTS
!     -------------------------------------------------------------------------------
!
!     Explicit/Implicit state update
!     0: implicit
!     1: explicit
!     By default set to be ON
      integer, parameter, public :: explicit = 1
!
!
!     Forward-Gradient predictor scheme
!     0: off
!     1: on
!     By default set to be OFF
      integer, parameter, public :: FGpredictor = 0
!
!     Factor used in forward gradient scheme
!     This is a factor used within the predictor scheme
!     0 : Euler solution
!     1 : Implicit solution
!     Default value is set to 0.5
      real(8), parameter, public :: theta = 0.5
!
!     Initial guess weight factor
!     This is used only when the predictor is turned off or,
!     when the predictor does not converge!
!     0 : fully plastic guess (stress at former time step, sigma_t)
!     1 : fully elastic guess (trial stress, sigma_tr)
!     guess = (1 - phi) * sigma_t + phi * sigma_tr
!     For sinh law slip 0.5 is recommended
!     For power law slip 0.0 is recommended
!     Default value is set to 0.1
      real(8), parameter, public :: phi = 0.1
!
!     Threshold for rss to crss ratio
!     This threshold becomes redundant when the sinh law is used!
!     Elastic solution will be used for lower values
!     Default value is set to 0.1
!     Do not choose this value greater than 0.5 for power law!
      real(8), parameter, public :: maxxcr = 0.1
!
!     Convergence tolerance for Cauchy stress of Newton-Raphson loop (absolute)
!     Default value is set to 1d-8 MPa
      real(8), parameter, public :: tolerance = 1.d-8
!
!     Convergence tolerance for state update of Newton-Raphson loop (absolute)
!     Default value is set to 1d-4 MPa
      real(8), parameter, public :: tauctolerance = 1.d-4
!
!     Maximum number of iterations
!     Default value is set to 200
      integer, parameter, public :: maxniter = 200
!
!     Increase the time step if there is no convergence
!     Default value is set to 0.5
      real(8), parameter, public :: cutback = 0.5
!
!     Increase the time step if there is no cutback
!     If selected as '0.' then ABAQUS determines the time step,
!     which is more efficient!
!     Default value is set to 1.25
!     Recommended value is 1.25 if manual setting is desired!
      real(8), parameter, public :: pastefront = 1.25
!
!     Inversion using Quad precision (16-bytes) in CP solver
!     0: OFF
!     1: ON
!     By default set to be ON
      integer, parameter, public :: quadprec = 1
!
!     Inversion by Singular Value Decomposition in CP solver
!     This is used when the inversion in CP solver cannot be
!     accomplished due to matrix to be inverted being close to 
!     singular (rank deficit)!
!     0: OFF
!     1: ON
!     By default set to be ON
      integer, parameter, public :: SVDinversion = 1
!
!
!     The size of parameter array
!     Please do not change this value
!     The PROPS correspondence must be updated in case it has changed
!     Default value is set to 14 (should give enough freedom)
      integer, parameter, public :: maxnparam = 14
!
!
!
!
!     MATERIAL INPUTS
!     -------------------------------------------------------------------------------
!
!     The material inputs are defined in "usermaterials.f"
!     Maximum number of materials available in the library
!     Materials refer to phases with different properties
!     Default value is set to the maximum available material id
      integer, parameter, public :: maxnmaterial = 10   
!
!     state variables per slip system
!     If turned on Forward-Gradient Predictor will be switched off
!     0: slip rate and strain hardening using individual slip systems
!     1: slip rate and strain hardening using the average properties
      integer,  parameter, public :: useaveragestatevars = 0
!
!
!     GND calculation identifier
!     0: No GNDs
!     1: curl(Fp) followed L2 method (Rank-deficit matrix inversion)
!     2: curl(Fp) followed L2 method (KKT solution)
!     3: curl(Fp) followed L2 method (A-matrix column elimination suggested by Chris Hardie) 
!     4: curl(gdotnFp) followed by projections
!     5: ||curl(gdotnFp)|| scalar GND
!     6: slip gradients followed by projections
!     Default value is set to 0
      integer, parameter, public :: gndmodel = 0
!
!     GND threshold
!     This is used to ignore GNDs below a threshold value
!     Used in gndmodels: 1-2-4-5-6 (except gndmodel=3)
!     Acts on the absolute GND increments
      real(8), parameter, public :: gndthreshold = 2.d-10
!
!     GND homogenization flag
!     0: GNDs at every integration point individually using extrapolation
!     1: GNDs at the element centers only (the same value for the integration points)
!     Not used when no GNDs are required
!     Default value is set to 0
      integer, parameter, public :: gndhomogenization = 0
!     
!     GND linear approximation flag - applicable to only quadratic elements
!     0: use the exact interpolation function of the element type (linear/quadratic)
!     1: use linear interpolation functions (only valid for quadratic elements)
      integer, parameter, public :: gndlinear = 0
!
!     Temperature flag
!     0: Temperature is defined by the field variable in ABAQUS in [K]
!     Material properties are entered once only at the initialization
!     1: Temperature is defined by the user (constant)
!     Material subroutine is entered every time step
      integer, parameter, public :: constanttemperature = 1
!
!
!     Temperature in Kelvins [K]
      real(8), parameter, public :: temperature = 298. 
!
!
!     Maximum number of slip systems in the mesh
!     among all possible materials in the mesh
!     Default value is set to 12
!     This is used to allocate arrays
!     Using a smaller number reduces the memory allocated
!     1. BCC material: 12/24
!     2. FCC material: 12/18 (if cubic slip is active)
!     3. HCP material: 3/6/12/24/30
!     4. Multiple phases: choose the highest value amongst the phases
      integer, parameter, public :: maxnslip = 12
!
!
!     The material inputs are defined in "usermaterials.f"
!     Maximum number of defect types for irradiation model-2
!     This is used to allocate arrays
      integer, parameter, public :: maxnloop = 3
!
!
!
      end module userinputs