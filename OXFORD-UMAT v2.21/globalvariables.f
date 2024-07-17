!     Sept. 20th, 2022
!     Eralp Demir
!     University of Oxford
!
!     Global variables used within the code
      module globalvariables
      implicit none
!
!
!
!     MATERIAL-LINKED GLOBAL VARIABLES
!     -------------------------------------------------------------------------------       
!     Global variable for Euler angles
      real(8), allocatable, public :: Euler(:,:)
      
!     Global variable storing grain id
!     Different grains can be present in the mesh
      integer, allocatable, public ::	 featureid(:,:)
!  
!     Global variable storing material id
      integer, allocatable, public ::	 materialid(:,:)
!
!     Different phases can be present in the mesh
      integer, allocatable, public ::	 phaseid(:,:)
!
!
!
!     Global variable storing phase-id
!     This refers to the crystal structure
!     Different phases can be present in the mesh
      integer, allocatable, public ::	 phaseid_all(:)
!
!     Global variable for number of slip systems
!     Number of slip systems could be different for each ip
      integer, allocatable, public ::	 numslip_all(:)
!
!     Global variable for number of slip systems
!     Required for GND calculations
!     Number of slip systems could be different for each ip
      integer, allocatable, public ::	 numscrew_all(:)
!
!     Slip model identifier
!     0: no slip (if only creep is considered)
!     1: sinh law
!     2: double exponent law
!     3: power law
      integer, allocatable, public :: slipmodel_all(:)
!
!     Slip parameters
!     Array size adjustable
      real(8), allocatable, public :: slipparam_all(:,:)
!     For sinh law (slipmodel=1)
!     1: alpha0 / constant value for alpha /  [1/s] / if a non-zero value is defined,
!     the alpha calculation will be ignored
!     2: beta0 / constant value for beta / [1/MPa] / if a non-zero value is defined,
!     the beta calculation will be ignored
!     3: psi / fraction of mobile dislocations / [-] / alpha calculation
!     4: rhom0 / reference mobile dislocation density / [1/micrometer^2] / alpha calculation
!     5: DeltaF / activation energy for slip / [J/mol] / alpha calculation
!     6: nu0 / attempt frequency / [1/s]  / alpha calculation
!     7: gamma0 / reference slip strain / [-] / beta calculation
!
!     For double exponent law (slipmodel=2)
!     1: gammadot0 / refence slip rate / [1/s]
!     2: p / inner exponent / [-]
!     3: q / outer exponent / [-]
!     4: DeltaFoct / activation energy for octahedral slip / [J/mol]
!     5: DeltaFcub / activation energy for cubic slip / [J/mol]
!
!     For power law (slipmodel=3)
!     1: gammadot0 / refence slip rate / [1/s]
!     2: n / power exponent / [-]
!
!
!     Creep model identifier
!     0: no creep
!     1: sinh slip strain
!     2: exponential slip strain
!     3: power law slip strain
!     4: sinh slip strain + climb strain
!     Default value is set to 0
      integer, allocatable, public :: creepmodel_all(:)
!     Creep parameters
!     Array size adjustable
      real(8), allocatable, public :: creepparam_all(:,:)
!
!     Hardening identifier
!     0: Hardening is off (tauc constant)
!     1: Kocks-Mecking hardening
!     2: Voce type hardening
!     Default value is set to 0
      integer, allocatable, public :: hardeningmodel_all(:)
!     Hardening parameters
!     Array size adjustable
      real(8), allocatable, public :: hardeningparam_all(:,:)
!
!     Irradiation identifier
!     0: Irradiaion is off
!     1: Irradiation is on
!     Default value is set to 0
      integer, allocatable, public :: irradiationmodel_all(:)
!     Irradiation model parameters
!     Array size adjustable
      real(8), allocatable, public :: irradiationparam_all(:,:)
!     1: tau_s0: initial solute strength
!     2: gamma_s: saturation value of the cumulative slip
!     3: psi / fraction of mobile density for irradiated material for sinh law
!
!
!     Backstress model parameters
!     Array size adjustable (maxnmaterial, maxnparam)
      real(8), allocatable, public :: backstressparam_all(:,:)
!
!     The following variables are stored for every element and ip
!     INITALLY - only once at the beginning
!     This is done for the following reasons
!     1. In order to reduce computation time
!     2. Various materials can be present in the mesh
!
!     Global variable for caratio
      real(8), allocatable, public ::	 caratio_all(:)
!
!     Global variable for cubicslip
      integer, allocatable, public ::	 cubicslip_all(:)
!
!     Global variable for elasticity in crystal reference
      real(8), allocatable, public ::	 Cc_all(:,:,:)
!
!     Global variable for geometric factor in Taylor equation
      real(8), allocatable, public ::	 gf_all(:)
!
!     Global variable for shear modulus
      real(8), allocatable, public ::	 G12_all(:)
!
!     Global variable for Poisson's ratio
      real(8), allocatable, public ::	 v12_all(:)
!
!     Global variable for thermal expansion matrix
      real(8), allocatable, public ::	 alphamat_all(:,:,:)
!
!     Global variable for Burgers vector
      real(8), allocatable, public ::	 burgerv_all(:,:)
!
!
!     INTERACTION MATRICES
!     Global variables for dislocation strength interaction matrix
      real(8), allocatable, public ::	 sintmat1_all(:,:,:)
!
!     Global variable for dislocation irradiation loop strength interaction matrix
      real(8), allocatable, public ::	 sintmat2_all(:,:,:)
!
!     Global variable for latent hardening interaction
      real(8), allocatable, public ::	 hintmat1_all(:,:,:)
!
!     Global variable for dislocation hardening interaction
      real(8), allocatable, public ::	 hintmat2_all(:,:,:)
!
!     Screw systems
      integer, allocatable, public ::	 screw_all(:,:)
!
!
!
!     -------------------------------------------------------------------------------
!
!
!
!
!     TIME (SOLUTION) DEPENDENT GLOBAL VARIABLES
!     -------------------------------------------------------------------------------
!     time at former time step
      real(8), public :: time_old
!
!     time increment at former time step
      real(8), public :: dt_t
!
!     Global initialization flag for elemental initialization
!     Orientation assigment
!     Initial slip direction calculations
!     Global coordinates
      integer, allocatable, public ::	initialized_firstinc(:,:)
!
!
!     Global initialization flag for IP initialization
      integer, allocatable, public ::	 ip_init(:,:)
!
!     Global one-time initialization flag
      integer, public              ::	 init_once
!
!     Global one-time initialization flag
      integer, public              ::	 grad_init
!
!     Global flag for IP count (10m elements, 100 ips)
      integer, allocatable, public ::	 ip_count(:,:)
!
!     Crystal orientation at current time step
      real(8), allocatable, public :: statev_gmatinv(:,:,:,:)
!     Crystal orientation at former time step
      real(8), allocatable, public :: statev_gmatinv_t(:,:,:,:)
!     Crystal orientation initially
      real(8), allocatable, public :: statev_gmatinv_0(:,:,:,:)
!
!     Total cumulative Von-Mises equivalent plastic strain at current time step
      real(8), allocatable, public :: statev_evmp(:,:)
!     Total cumulative Von-Mises equivalent plastic strain at former time step
      real(8), allocatable, public :: statev_evmp_t(:,:)
!
!     Plastic dissipation power density at current time step
      real(8), allocatable, public :: statev_plasdiss(:,:)
!     Plastic dissipation power density at former time step
      real(8), allocatable, public :: statev_plasdiss_t(:,:)
!
!     Effective critical resolved shear stress at current time step
      real(8), allocatable, public :: statev_tauceff(:,:,:)
!
!
!     Total cumulative slip at current time step
      real(8), allocatable, public :: statev_totgammasum(:,:)
!     Total cumulative slip at former time step
      real(8), allocatable, public :: statev_totgammasum_t(:,:)
!
!     Cumulative slip at current time step
      real(8), allocatable, public :: statev_gammasum(:,:,:)
!     Cumulative slip at former time step
      real(8), allocatable, public :: statev_gammasum_t(:,:,:)
!
!     Slip rates at current time step
      real(8), allocatable, public :: statev_gammadot(:,:,:)
!     Slip rates at former time step
      real(8), allocatable, public :: statev_gammadot_t(:,:,:)
!
!
!     Plastic part of the deformation gradient at current time step
      real(8), allocatable, public :: statev_Fp(:,:,:,:)
!     Plastic part of the deformation gradient at former time step
      real(8), allocatable, public :: statev_Fp_t(:,:,:,:)
!
!
!     Thermal part of the deformation gradient at current time step
      real(8), allocatable, public :: statev_Fth(:,:,:,:)
!     Thermal part of the deformation gradient at former time step
      real(8), allocatable, public :: statev_Fth_t(:,:,:,:)
!
!
!     Cauchy stress at current time step
      real(8), allocatable, public :: statev_sigma(:,:,:)
!     Cauchy stress at former time step
      real(8), allocatable, public :: statev_sigma_t(:,:,:)
!     Cauchy stress at previous time step
      real(8), allocatable, public :: statev_sigma_t2(:,:,:)
!
!     Cauchy stress at current time step
      real(8), allocatable, public :: statev_jacobi(:,:,:,:)
!     Cauchy stress at former time step
      real(8), allocatable, public :: statev_jacobi_t(:,:,:,:)
!
!
!     Critical resolved shear stress at current time step
      real(8), allocatable, public :: statev_tauc(:,:,:)
!     Critical resolved shear stress at former time step
      real(8), allocatable, public :: statev_tauc_t(:,:,:)
!
!     maximum of the tau/tauc ratio at current time step
      real(8), allocatable, public :: statev_maxx(:,:)
!     maximum of tau/tauc ratio at former time step
      real(8), allocatable, public :: statev_maxx_t(:,:)
!
!     elastic strains at the crystal ref. at current time step
      real(8), allocatable, public :: statev_Eec(:,:,:)
!     elastic strains at the crystal ref. at former time step
      real(8), allocatable, public :: statev_Eec_t(:,:,:)
!
!     Lattice curvature at current time step
      real(8), allocatable, public :: statev_curvature(:,:,:)
!
!     Incompatibility at current time step
      real(8), allocatable, public :: statev_Lambda(:,:,:)
!     Incompatibility at former time step
      real(8), allocatable, public :: statev_Lambda_t(:,:,:)
!
!     GND density per slip system at current time step
      real(8), allocatable, public :: statev_gnd(:,:,:)
!     GND density per slip system at former time step
      real(8), allocatable, public :: statev_gnd_t(:,:,:)
!     GND density per slip system initially - EBSD import
      real(8), allocatable, public :: statev_gnd_0(:,:,:)
!
!     SSD density per slip system at current time step
      real(8), allocatable, public :: statev_ssd(:,:,:)
!     SSD density per slip system at former time step
      real(8), allocatable, public :: statev_ssd_t(:,:,:)
!
!     Total SSD density at current time step
      real(8), allocatable, public :: statev_ssdtot(:,:)
!     Total SSD density at former time step
      real(8), allocatable, public :: statev_ssdtot_t(:,:)
!
!     Forest dislocation density per slip system at current time step
      real(8), allocatable, public :: statev_forest(:,:,:)
!     Forest dislocation density per slip system at former time step
      real(8), allocatable, public :: statev_forest_t(:,:,:)
!
!     Substructure dislocation density for irradiation per slip system at current time step
      real(8), allocatable, public :: statev_substructure(:,:)
!     Substructure dislocation density for irradiation per slip system at former time step
      real(8), allocatable, public :: statev_substructure_t(:,:)
!
!     Solute strength for irradiation per slip system at current time step
      real(8), allocatable, public :: statev_tausolute(:,:)
!     Solute strength for irradiation per slip system at former time step
      real(8), allocatable, public :: statev_tausolute_t(:,:)
!
!
!     Defect loop density for irradiation per slip system at current time step
      real(8), allocatable, public :: statev_loop(:,:,:)
!     Defect loop for irradiation per slip system at former time step
      real(8), allocatable, public :: statev_loop_t(:,:,:)
!
!     Backstress per slip system at current time step
      real(8), allocatable, public :: statev_backstress(:,:,:)
!     Backstress per slip system at former time step
      real(8), allocatable, public :: statev_backstress_t(:,:,:)
!
!     -------------------------------------------------------------------------------
!
!
!
!     MESH INPUTS
!     -------------------------------------------------------------------------------
!
!     Total number of elements in the mesh
!     Adaptive mesh cannot be used!
      integer, public :: numel
!
!     Element type
!     Single element type is possible throughout the mesh
!     Please do not leave any spaces between the characters
!     Use the element types defined in ABAQUS element library
!     Please see meshprop.f for the list of available element types
      character(len=:), allocatable, public :: eltyp
!
!
!     Number of integration points per element
      integer, public ::  numpt
!
!     Number of nodes per element
      integer, public ::  nnpel
!
!
!     Dimensions of the problem:
!     Tensor dimensions in UMAT
      integer, public ::  numtens
!
!     2: 2D / 3: 3D
      integer, public ::  numdim
!
!
!     -------------------------------------------------------------------------------
!
!
!
!     -------------------------------------------------------------------------------
!
!     GLOBAL VARIABLES FOR NON-LOCAL CALCULATIONS
!     -------------------------------------------------------------------------------
!     These variables will be assigned at the initialization
!     A single type of element is assumed throughout the mesh
!
!     integration point domain (area/volume)
      real(8), allocatable, public ::	 ipdomain(:,:)
!
!     integration point weights
      real(8), allocatable, public :: ipweights(:)
!
!     Global integration point coordinates
      real(8), allocatable, public ::	 ipcoords(:,:,:)
!
!     Element specific shape functions and mappings
!     Dependent on element type
!
!     Gradient map for elements
      real(8), allocatable, public :: gradip2ip(:,:,:,:)
!
!     Interpolation matrix
      real(8), allocatable, public :: Nmat(:,:)
!
!     Inverse of interpolation matrix
      real(8), allocatable, public :: invNmat(:,:)
!
!     Shape function derivatives
      real(8), allocatable, public :: dNmat(:,:,:)
!
!     Shape function derivatives at element center
      real(8), allocatable, public :: dNmatc(:,:)
!
!     Element having a single IP
      integer, public :: calculategradient
!
!     -------------------------------------------------------------------------------
!
!
!
!
!     CONSTANTS
!     -------------------------------------------------------------------------------    
!
!	Number pi
	real(8), parameter, public :: pi = 3.14159265359
!
!	Taylor factor for a polycrytal aggregate
	real(8), parameter, public :: TF = 3.1
!
!     Universal gas contant [J/mol/K]
      real(8), parameter, public :: Rgas = 8.31432
!
!     Boltzman contant [m2 kg s-2 K-1 ]
      real(8), parameter, public :: KB = 1.380649d-23
!
!	Small real number
	real(8), parameter, public :: smallnum = 1.0d-20
!
!	Large real number
	real(8), parameter, public :: largenum = 1.0d+50
!     -------------------------------------------------------------------------------
!
!
!
!
!     IDENTITY TENSORS
!     -------------------------------------------------------------------------------
!
!	Identity matrix (3x3)
	real(8), public :: I3(3,3)
!
!	Identity matrix (6x6)
      real(8), public :: I6(6,6)
!
!	Identity matrix (9x9)
      real(8), public :: I9(9,9)
!
!	Permutation symbol (3,3,3)
      real(8), public :: eijk(3,3,3)
!
!
!     SLIP DIRECTIONS
!     ------------------------------------------------------------------------------- 
!
!
!     Global variable for slip directions
!     Slip direction can be different for each material
      real(8), allocatable, public ::	 dirc_0_all(:,:,:)
!
!     Global variable for slip plane normal
!     Slip plane normal can be different for each material
      real(8), allocatable, public ::	 norc_0_all(:,:,:)
!
!     Global variable for transverse direction
!     Transverse direction can be different for each material
      real(8), allocatable, public ::	 trac_0_all(:,:,:)
!
!     Global variable for line direction
!     Line direction can be different for each material
      real(8), allocatable, public ::	 linc_0_all(:,:,:)
!
!     Global variable for initial Schmid tensor at the sample referencec
!     Schmid tensor can be different for each material
      real(8), allocatable, public ::	 Schmid_0_all(:,:,:,:)
!
!     Global variable for forest projections
!     Forest projection for GND and SSD
!     Can be different for each material
      real(8), allocatable, public ::	 forestproj_all(:,:,:)     
!
!     Global variable to map screw dislocations from the corespondent slip systems
!     This is needed for gndmodels 4/5/6 where screw dislocations are computed at each system
!     Therefore, need to account for only the defined set of screw systems
!     Can be different for each material
      real(8), allocatable, public ::	 slip2screw_all(:,:,:)
!
!
!
!     BCC slip directions
      real(8), public ::  dir1(48,3)
!     <111> {110} slip family
      data dir1(1,:)   /-1.,  1.,  1. /
	data dir1(2,:)   / 1., -1.,  1. /
	data dir1(3,:)   / 1.,  1.,  1. /
	data dir1(4,:)   / 1., -1.,  1. /
	data dir1(5,:)   / 1.,  1.,  1. /
	data dir1(6,:)   / 1.,  1., -1. /
	data dir1(7,:)   / 1.,  1.,  1. /
	data dir1(8,:)   /-1.,  1.,  1. /
	data dir1(9,:)   / 1., -1.,  1. /
	data dir1(10,:)   /1.,  1., -1. /
	data dir1(11,:)  /-1.,  1.,  1. /
	data dir1(12,:)  / 1.,  1., -1. /
!     <111> {112} slip family
	data dir1(13,:)  / 1.,  1., -1. /
	data dir1(14,:)  / 1., -1.,  1. /
	data dir1(15,:)  /-1.,  1.,  1. /
	data dir1(16,:)  / 1.,  1.,  1. /
	data dir1(17,:)  / 1., -1.,  1. /
	data dir1(18,:)  / 1.,  1., -1. /
	data dir1(19,:)  / 1.,  1.,  1. /
	data dir1(20,:)  /-1.,  1.,  1. /
	data dir1(21,:)  /-1.,  1.,  1. /
	data dir1(22,:)  / 1.,  1.,  1. /
	data dir1(23,:)  / 1.,  1., -1. /
	data dir1(24,:)  / 1., -1.,  1. /
!     <111> {123} slip family
      data dir1(25,:) / 1.,  1., -1. /
      data dir1(26,:) / 1., -1.,  1. /
      data dir1(27,:) /-1.,  1.,  1. /
      data dir1(28,:) / 1.,  1.,  1. /
      data dir1(29,:) / 1., -1.,  1. /
      data dir1(30,:) / 1.,  1., -1. /
      data dir1(31,:) / 1.,  1.,  1. /
      data dir1(32,:) /-1.,  1.,  1. /
      data dir1(33,:) / 1.,  1., -1. /
      data dir1(34,:) / 1., -1.,  1. /
      data dir1(35,:) /-1.,  1.,  1. /
      data dir1(36,:) / 1.,  1.,  1. /
      data dir1(37,:) / 1., -1.,  1. /
      data dir1(38,:) / 1.,  1., -1. /
      data dir1(39,:) / 1.,  1.,  1. /
      data dir1(40,:) /-1.,  1.,  1. /
      data dir1(41,:) /-1.,  1.,  1. /
      data dir1(42,:) / 1.,  1.,  1. /
      data dir1(43,:) / 1.,  1., -1. /
      data dir1(44,:) / 1., -1.,  1. /
      data dir1(45,:) /-1.,  1.,  1. /
      data dir1(46,:) / 1.,  1.,  1. /
      data dir1(47,:) / 1.,  1., -1. /
      data dir1(48,:) / 1., -1.,  1. /
!
!
!     BCC slip plane normals
      real(8), public ::  nor1(48,3)
!     <111> {110} slip family
      data nor1(1,:)   / 1.,  1.,  0. /
      data nor1(2,:)   / 1.,  1.,  0. /
      data nor1(3,:)   /-1.,  0.,  1. /
      data nor1(4,:)   /-1.,  0.,  1. /
      data nor1(5,:)   /-1.,  1.,  0. /
      data nor1(6,:)   /-1.,  1.,  0. /
      data nor1(7,:)   / 0.,  1., -1. /
      data nor1(8,:)   / 0.,  1., -1. /
      data nor1(9,:)   / 0.,  1.,  1. /
      data nor1(10,:)  / 0.,  1.,  1. /
      data nor1(11,:)  / 1.,  0.,  1. /
      data nor1(12,:)  / 1.,  0.,  1. /
!     <111> {112} slip family
      data nor1(13,:)  / 1.,  1.,  2. /
      data nor1(14,:)  /-1.,  1.,  2. /
      data nor1(15,:)  / 1., -1.,  2. /
      data nor1(16,:)  / 1.,  1., -2. /
      data nor1(17,:)  / 1.,  2.,  1. /
      data nor1(18,:)  /-1.,  2.,  1. /
      data nor1(19,:)  / 1., -2.,  1. /
      data nor1(20,:)  / 1.,  2., -1. /
      data nor1(21,:)  / 2.,  1.,  1. /
      data nor1(22,:)  /-2.,  1.,  1. /
      data nor1(23,:)  / 2., -1.,  1. /
      data nor1(24,:)  / 2.,  1., -1. /
!     <111> {123} slip family
      data nor1(25,:) /  1.,  2.,  3. /
      data nor1(26,:) / -1.,  2.,  3. /
      data nor1(27,:) /  1., -2.,  3. /
      data nor1(28,:) /  1.,  2., -3. /
      data nor1(29,:) /  1.,  3.,  2. /
      data nor1(30,:) / -1.,  3.,  2. /
      data nor1(31,:) /  1., -3.,  2. /
      data nor1(32,:) /  1.,  3., -2. /
      data nor1(33,:) /  2.,  1.,  3. /
      data nor1(34,:) / -2.,  1.,  3. /
      data nor1(35,:) /  2., -1.,  3. /
      data nor1(36,:) /  2.,  1., -3. /
      data nor1(37,:) /  2.,  3.,  1. /
      data nor1(38,:) / -2.,  3.,  1. /
      data nor1(39,:) /  2., -3.,  1. /
      data nor1(40,:) /  2.,  3., -1. /
      data nor1(41,:) /  3.,  1.,  2. /
      data nor1(42,:) / -3.,  1.,  2. /
      data nor1(43,:) /  3., -1.,  2. /
      data nor1(44,:) /  3.,  1., -2. /
      data nor1(45,:) /  3.,  2.,  1. /
      data nor1(46,:) / -3.,  2.,  1. /
      data nor1(47,:) /  3., -2.,  1. /
      data nor1(48,:) /  3.,  2., -1. /      
!
!
!
!     FCC slip directions
      real(8), public ::  dir2(18,3)
!     <110> {111} slip family
	data dir2(1,:)  / 1., -1.,  0. /
	data dir2(2,:)  / 0.,  1., -1. /
	data dir2(3,:)  / 1.,  0., -1. /
	data dir2(4,:)  / 1.,  1.,  0. /
	data dir2(5,:)  / 0.,  1., -1. /
	data dir2(6,:)  / 1.,  0.,  1. /
	data dir2(7,:)  / 1.,  1.,  0. /
	data dir2(8,:)  / 0.,  1.,  1. /
	data dir2(9,:)  / 1.,  0., -1. /
	data dir2(10,:) / 1., -1.,  0. /
	data dir2(11,:) / 0.,  1.,  1. /
	data dir2(12,:) / 1.,  0.,  1. /
!     <110> {100} cubic slip family
      data dir2(13,:) / 0.,  1.,  1. /
      data dir2(14,:) / 0.,  1., -1. /
      data dir2(15,:) / 1.,  0.,  1. /
      data dir2(16,:) / 1.,  0., -1. /
      data dir2(17,:) / 1.,  1.,  0. /
      data dir2(18,:) / 1., -1.,  0. /
!
!
!     FCC slip plane normals
      real(8), public ::  nor2(18,3)
!     <110> {111} slip family
	data nor2(1,:)  / 1.,  1.,  1. /
	data nor2(2,:)  / 1.,  1.,  1. /
	data nor2(3,:)  / 1.,  1.,  1. /
	data nor2(4,:)  /-1.,  1.,  1. /
	data nor2(5,:)  /-1.,  1.,  1. /
	data nor2(6,:)  /-1.,  1.,  1. /
	data nor2(7,:)  / 1., -1.,  1. /
	data nor2(8,:)  / 1., -1.,  1. /
	data nor2(9,:)  / 1., -1.,  1. /
	data nor2(10,:) / 1.,  1., -1. /
	data nor2(11,:) / 1.,  1., -1. /
	data nor2(12,:) / 1.,  1., -1. /
!     <110> {100} cubic slip family
      data nor2(13,:) / 1.,  0.,  0. /
      data nor2(14,:) / 1.,  0.,  0. /
      data nor2(15,:) / 0.,  1.,  0. /
      data nor2(16,:) / 0.,  1.,  0. /
      data nor2(17,:) / 0.,  0.,  1. /
      data nor2(18,:) / 0.,  0.,  1. /
!
!
!     The directions are corrected by Alvaro 23.02.2023
!     HCP slip directions
      real(8), public ::  dir3h(30,4)
!     <-1-1.0>{00.1} / basal systems (independent of c/a-ratio)
      data dir3h(1,:)   / 2., -1., -1.,  0./
      data dir3h(2,:)   /-1.,  2., -1.,  0./
      data dir3h(3,:)   /-1., -1.,  2.,  0./
!     <-1-1.0>{1-1.0} / prismatic systems (independent of c/a-ratio)
      data dir3h(4,:)   / 2., -1., -1.,  0./
      data dir3h(5,:)   /-1.,  2., -1.,  0./
      data dir3h(6,:)   /-1., -1.,  2.,  0./
!     <-1-1.0>{-11.1} / 1st order pyramidal <a> systems (direction independent of c/a-ratio)
      data dir3h(7,:)  /-1.,  2., -1.,  0./
      data dir3h(8,:)  /-2.,  1.,  1.,  0./
      data dir3h(9,:)  /-1., -1.,  2.,  0./
      data dir3h(10,:)  / 1., -2.,  1.,  0./
      data dir3h(11,:)  / 2., -1., -1.,  0./
      data dir3h(12,:)  / 1.,  1., -2.,  0./
!     <11.3>{-10.1} / 1st order pyramidal <c+a> systems (direction independent of c/a-ratio)
      data dir3h(13,:)  /-2.,  1.,  1.,  3./
      data dir3h(14,:)  /-1., -1.,  2.,  3./
      data dir3h(15,:)  /-1., -1.,  2.,  3./
      data dir3h(16,:)  / 1., -2.,  1.,  3./
      data dir3h(17,:)  / 1., -2.,  1.,  3./
      data dir3h(18,:)  / 2., -1., -1.,  3./
      data dir3h(19,:)  / 2., -1., -1.,  3./
      data dir3h(20,:)  / 1.,  1., -2.,  3./
      data dir3h(21,:)  / 1.,  1., -2.,  3./
      data dir3h(22,:)  /-1.,  2., -1.,  3./
      data dir3h(23,:)  /-1.,  2., -1.,  3./
      data dir3h(24,:)  /-2.,  1.,  1.,  3./
!
!     <11.3>{-1-1.2} / 2nd order pyramidal <c+a> systems
      data dir3h(25,:)   /-1., -1.,  2.,  3./
      data dir3h(26,:)   / 1., -2.,  1.,  3./
      data dir3h(27,:)   / 2., -1., -1.,  3./
      data dir3h(28,:)  / 1.,  1., -2.,  3./
      data dir3h(29,:)  /-1.,  2., -1.,  3./
      data dir3h(30,:)  /-2.,  1.,  1.,  3./
!
!      
!     HCP slip plane normals
      real(8), public ::  nor3h(30,4)
!     <-1-1.0>{00.1} / basal systems (independent of c/a-ratio)
      data nor3h(1,:)  /0.,  0.,  0.,  1./
      data nor3h(2,:)  /0.,  0.,  0.,  1./
      data nor3h(3,:)  /0.,  0.,  0.,  1./
!     <-1-1.0>{1-1.0} / prismatic systems (independent of c/a-ratio)
      data nor3h(4,:)  /0.,  1., -1.,  0./
      data nor3h(5,:)  /-1.,  0.,  1.,  0./
      data nor3h(6,:)  /1., -1.,  0.,  0./
!     <-1-1.0>{-11.1} / 1st order pyramidal <a> systems (direction independent of c/a-ratio)
      data nor3h(7,:)  /1.,  0., -1.,  1./
      data nor3h(8,:)  /0.,  1., -1.,  1./
      data nor3h(9,:)  /-1.,  1.,  0., 1./
      data nor3h(10,:)  /-1.,  0.,  1., 1./
      data nor3h(11,:)  /0., -1.,  1.,  1./
      data nor3h(12,:)  /1., -1.,  0.,  1./
!     <11.3>{-10.1} / 1st order pyramidal <c+a> systems (direction independent of c/a-ratio)
      data nor3h(13,:)  /1.,  0., -1.,  1./
      data nor3h(14,:)  /1.,  0., -1.,  1./
      data nor3h(15,:)  /0.,  1., -1.,  1./
      data nor3h(16,:)  /0.,  1., -1.,  1./
      data nor3h(17,:)  /-1.,  1.,  0.,  1./
      data nor3h(18,:)  /-1.,  1.,  0.,  1./
      data nor3h(19,:)  /-1.,  0.,  1.,  1./
      data nor3h(20,:)  /-1.,  0.,  1.,  1./
      data nor3h(21,:)  /0., -1.,  1.,  1./
      data nor3h(22,:)  /0., -1.,  1.,  1./
      data nor3h(23,:)  /1., -1.,  0.,  1./
      data nor3h(24,:)  /1., -1.,  0.,  1./
!     <11.3>{-1-1.2} / 2nd order pyramidal <c+a> systems
      data nor3h(25,:)  /1.,  1., -2.,  2./
      data nor3h(26,:)  /-1.,  2., -1.,  2./
      data nor3h(27,:)  /-2.,  1.,  1.,  2./
      data nor3h(28,:)  /-1., -1.,  2.,  2./
      data nor3h(29,:)  /1., -2.,  1.,  2./
      data nor3h(30,:)  /2., -1., -1.,  2./
!
!
!     Slip direction for alpha-Uranium
!     see McCabe, Tome et al 2010 (Figure 2)
      real(8), public ::  dir4(8,3)
      data dir4(1,:)  /1.0, 0.0, 0.0/
      data dir4(2,:)  /1.0, 0.0, 0.0/
	data dir4(3,:)  /0.43731, -0.89931, 0.0/ ! [1-10](110) -> [a,-b,0](b,a,0)
	data dir4(4,:)  /0.43731,0 .89931, 0.0/ ! [110](1-10) -> [a,b,0](b,-a,0)
	data dir4(5,:)  /0.24074, -0.49507, 0.83483/ ! [1-12](021) -> [a,-b,2c](0,c,b/2)
	data dir4(6,:)  /-0.24074, -0.49507, 0.83483/ ! [-1-12](021) -> [-a,-b,2c](0,c,b/2)
	data dir4(7,:)  /0.24074, 0.49507, 0.83483/ ! [112](0-21) -> [a,b,2c](0,c,-b/2)
	data dir4(8,:)  /0.24074, -0.49507, -0.83483/ ! [1-1-2](0-21) -> [a,-b,-2c](0,c,-b/2)
!
!
!     Slip plane normals of alpha-Uranium
!     see McCabe, Tome et al 2010 (Figure 2)
      real(8), public ::  nor4(8,3)
	data nor4(1,:)  /0.0, 1.0, 0.0/
	data nor4(2,:)  /0.0, 0.0, 1.0/
	data nor4(3,:)  /0.89931, 0.43731, 0.0/ ! [1-10](110) -> [a,-b,0](b,a,0)
	data nor4(4,:)  /0.89931, -0.43731, 0.0/ ! [110](1-10) -> [a,b,0](b,-a,0)
	data nor4(5,:)  /0.0, 0.86013, 0.51008/ ! [1-12](021) -> [a,-b,2c](0,c,b/2)
	data nor4(6,:)  /0.0, 0.86013, 0.51008/ ! [-1-12](021) -> [-a,-b,2c](0,c,b/2)
	data nor4(7,:)  /0.0, 0.86013, -0.51008/ ! [112](0-21) -> [a,b,2c](0,c,-b/2)
	data nor4(8,:)  /0.0, 0.86013, -0.51008/ ! [1-1-2](0-21) -> [a,-b,-2c](0,c,-b/2)         
!
!
!
!
!     -------------------------------------------------------------------------------
!
      end module globalvariables