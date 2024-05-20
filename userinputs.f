!     Sept. 26th, 2022
!     Eralp Demir
!
!     User inputs
!     Global variables that the user entry is required
!     Users asked to change the model and numerical parameters here
!
!
      module userinputs
      implicit none
!
!
!  
!
!     NUMERICAL INPUTS
!     -------------------------------------------------------------------------------      
!     Explicit/Implicit state update
!     0: explicit
!     1: semi-implicit
!     By default set to be ON
      integer, parameter, public :: stateupdate = 0
!
!
!     Predictor schemes
!     0: Weighted average of fully-plastic and fully-elastic guess
!     1: Stress extrapolation (Chris Hardie)
!     By default set to "1"
      integer, parameter, public :: predictor = 0
!
!
!
!     Inverted Newton Loop as backup predictor/solution 
!     Hardie et al. 2023: https://www.sciencedirect.com/science/article/pii/S0749641923002577
!     0: Don't use
!     1: Use
      integer, parameter, public :: inversebackup = 0
!
!     Tolerance for the reverse scheme
      real(8), parameter, public :: inversetolerance = 10.0 ! MPa.s
!
!
!
!     Factor used in forward gradient scheme
!     This is a factor used within the predictor scheme
!     0 : Euler solution
!     1 : Implicit solution
!     Default value is set to 1.
      real(8), parameter, public :: theta = 1.
!
!     Initial guess weight factor
!     This is used only when the predictor is turned off or,
!     when the predictor does not converge!
!     0 : fully plastic guess (stress at former time step, sigma_t)
!     1 : fully elastic guess (trial stress, sigma_tr)
!     guess = (1 - phi) * sigma_t + phi * sigma_tr
!     For sinh law slip 0.5 is recommended
!     For power law slip 0.0 is recommended
!     Default value is set to 0.
      real(8), parameter, public :: phi = 0.
!
!     Threshold for rss to crss ratio
!     This threshold becomes redundant when the sinh law is used!
!     Elastic solution will be used for lower values
!     Default value is set to 0.1
!     Do not choose this value greater than 0.0 for power law!
      real(8), parameter, public :: maxxcr = 0.0
!
!     Convergence tolerance for Cauchy stress of Newton-Raphson loop (absolute)
!     Default value is set to 1d-6 MPa
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
      real(8), parameter, public :: pastefront = 0.
!
!
!     Inversion by Singular Value Decomposition in CP solver
!     This is used when the inversion in CP solver cannot be
!     accomplished due to matrix to be inverted being close to
!     singular (rank deficit)!
!     0: OFF
!     1: ON
!     By default set to be OFF (Problems may occur with array allocation)
      integer, parameter, public :: SVDinversion = 0
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
!     In-grain orientation scatter
!     Need to provide .vox file name that has the following format
!     Default value is set to "0"
      integer, parameter, public :: readmaterialfile = 0
!
!     Material file that tells orientation of each element
!     The file has the following format (for each element):
!     phi1, PHI, phi2, x-div., y-div., z-div., grainID, phaseID
      character(len=*), parameter, public :: 
     + voxfilename = 'Job-1.vox'
!     
!     Folder location
!     This is necessary if there is an extra input file
!     Ex. vox file for in-grain orientation scatter
      character(len=*), parameter, public :: foldername=
     +'C:\Users\engs2474\Documents\Eralp
     + \Oxford\UMAT\Test Cases\0. version test'   
!
!     The material inputs are defined in "usermaterials.f"
!     Maximum number of materials available in the library
!     Materials refer to phases with different properties
!     Default value is set to the maximum available material id
      integer, parameter, public :: maxnmaterial = 12
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
!     1: (curlxFp)^T followed L2 method - SVD - Active SS restriction (new)
!     2: Gurtin's measure for incompatibility - L2 - SVD - Active SS restriction (Cermeli, P. and Gurtin, M.E., 2001)
!     3: curl(gdotnFp) followed by projections (H.Dai, PhD Thesis, 1997.)
!     4: slip gradients followed by projections (Gerken, J.M. and Dawson, P.R., 2008.)
!     Default value is set to 0
      integer, parameter, public :: gndmodel = 0
!
!     GND threshold
!     This is used to ignore GNDs below a threshold value
!     Used in all gndmodels
!     Acts on the absolute GND increments
      real(8), parameter, public :: gndthreshold = 2.d-10
!
!     Slip rate threshold for GND calculation
!     This is used as the absolute limit below which slip is ignored
!     Used in gndmodels: 1 and 2 (SVD inversion)
!     Acts on the absolute value of total slip
!     Only valid for models 1 and 2
      real(8), parameter, public :: slipthreshold = 1.d-10
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
!     Backstress flag
!     0: No Backstress
!     1: Local backstress model (Armstron-Frederick)
!     2: Non-local backstress model based on GNDs
!     (2 will only be effect in case any of the GND model is active)
!     Default value is set to 0
      integer, parameter, public :: backstressmodel = 0
!
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
!     1. BCC material: 12/24/48
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
!     Maximum number of elements in the mesh
      integer, parameter, public :: maxnumel = 10000000
!
!
!     Maximum number of integration points per element
      integer, parameter, public :: maxnumpt = 27
!
!
!
      end module userinputs