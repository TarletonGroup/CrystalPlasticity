!     Sept. 24th, 2022
!
!     Everytime UMAT enters here to get temperature dependent properties
!     User has to enter the inputs manually
      module usermaterials
      implicit none
!
!
!
!
      contains
!
!     Material parameters are according to following IDs      
!     material-01: custom material / bcc / phase-1
!     material-02: custom material / fcc / phase-2
!     material-03: custom material / hcp / phase-3
!     material-04: tungsten / bcc / phase-1
!     material-05: copper / fcc / phase-2
!     material-06: carbide / fcc / phase-2
!     material-07: CSMX-4 Nickel alloy / fcc / phase-2 ==> cubicslipystem flag
!     material-08: zirconium / hcp / phase-3
!     material-09: berylium / hcp / phase-3 (not ready yet)
!     material-10: alpha-uranium / alphauranium / phase-4
!     material-11: copper / UKAEA / Vikram Phalke
!     material-12: CuCrZr / UKAEA / Vikram Phalke
!
!
!
!     **********************************************
!     ** MATERIALPARAMETERS sets the material     **
!     ** constants for elasticity and plasticity  **
!     **********************************************
      subroutine materialparam(imat,temperature,
     + iphase,nslip,nscrew,caratio,cubicslip,Cc,
     + gf,G12,v12,alphamat,burgerv,
     + tauc_0,rho_0,rhofor_0,rhosub_0,
     + slipmodel,slipparam,creepmodel,creepparam,
     + hardeningmodel,hardeningparam,
     + irradiationmodel,irradiationparam,
     + sintmat1,sintmat2,hintmat1,hintmat2,
     + backstressparam)
      use errors, only : error
      use userinputs, only : maxnslip, maxnparam
      implicit none
!
!     material id
      integer, intent(in) :: imat
!
!
!     current temperature in Kelvins
      real(8), intent(in) :: temperature
!	  
!     phase id
!     1: BCC, 2: FCC, 3: HCP, 4: alpha-uranium
      integer, intent(out) :: iphase
!
!     number of slip systems
      integer, intent(out) :: nslip
!
!     number of screw systems
      integer, intent(out) :: nscrew
!      
!     c/a ratio for hcp crystals
      real(8), intent(out) :: caratio
!
!     cubic slip for fcc superalloys
      integer, intent(out) :: cubicslip
!
!     elastic stiffness matrix in the crystal reference frame
      real(8), intent(out) :: Cc(6,6)
!
!     shear modulus for Taylor's dislocation law
      real(8), intent(out) :: G12
!
!     Poisson's ratio
      real(8), intent(out) :: v12
!
!     geometric factor for obstacle strength
      real(8), intent(out) :: gf
!
!     thermal eigenstrain to model thermal expansion
      real(8), intent(out) :: alphamat(3,3)
!
!     burgers vectors
      real(8), intent(out) :: burgerv(maxnslip)
!
!     critical resolved shear stress of slip systems
      real(8), intent(out) :: tauc_0(maxnslip)
!
!     initial ssd dislocation density 
      real(8), intent(out) :: rho_0(maxnslip)
!
!     initial forest dislocation density 
      real(8), intent(out) :: rhofor_0
!
!     initial substructure dislocation density 
      real(8), intent(out) :: rhosub_0
!
!     Slip model identifier
!     0: no slip (if only creep is considered)
!     1: sinh law
!     2: double exponent law      
!     3: power law
      integer, intent(out) :: slipmodel
!
!     Slip parameters
!     Array size adjustable
      real(8), intent(out) :: slipparam(maxnparam)
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
!     3: dn/dT / temperature dependence of slip rate sensitivity / [1/K]
!
!
!     Creep model identifier
!     0: no creep
!     1: exponential creep law
!     Default value is set to 0
      integer, intent(out) :: creepmodel
!     Creep parameters
!     Array size adjustable
      real(8), intent(out) :: creepparam(maxnparam)
!
!     Hardening identifier
!     0: Hardening is off (tauc constant)
!     1: Voce type hardening
!     2: Linear hardening
!     3: Kocks-Mecking hardening
!     4: Kocks-Mecking hardening with substructure
!     Default value is set to 0
      integer, intent(out) :: hardeningmodel
!     Hardening parameters
!     Array size adjustable
      real(8), intent(out) :: hardeningparam(maxnparam)
!
!     Irradiation identifier
!     0: Irradiaion is off
!     1: Irradiation is on
!     Default value is set to 0
      integer, intent(out) :: irradiationmodel
!     Irradiation model parameters
!     Array size adjustable
      real(8), intent(out) :: irradiationparam(maxnparam)
!     1: tau_s0: initial solute strength
!     2: gamma_s: saturation value of the cumulative slip
!     3: psi / fraction of mobile density for irradiated material for sinh law
!
!     Interaction matrices
!     Strength interaction between dislocations
      real(8), intent(out) :: sintmat1(maxnslip,maxnslip)
!     Strength interaction dislocation loops related with irradiation
      real(8), intent(out) :: sintmat2(maxnslip,maxnslip)
!     Latent hardening
      real(8), intent(out) :: hintmat1(maxnslip,maxnslip)
!     Hardening interaction matrix between dislocations
      real(8), intent(out) :: hintmat2(maxnslip,maxnslip)
!
!     Slip parameters
!     Array size adjustable
      real(8), intent(out) :: backstressparam(maxnparam)
!
!     burgers vector scalars
      real(8) :: burger1, burger2
!
!     elastic constants scalars
      real(8) :: e1, e2, e3, g13, v13
      real(8) :: g23, v23, v21, v31, v32
!
!     critical resolved shear stress
      real(8) :: xtauc1, xtauc2, xtauc3, xtauc4, xtauc5
!
!     thermal expansion coefficients
      real(8) :: alpha1, alpha2, alpha3
!
!     temperature in celsius
      real(8) :: tcelsius
!
      real(8) :: C11, C12, C44, cst
!
!     local variable for identity martrix
      real(8) :: kdelta(maxnslip,maxnslip), q1, q2
!
      integer :: i, j, k, is, js
!
!
!
!     Set potentially unassigned parameters to zero
!     Slip parameters
      slipparam = 0.
!     Creep parameters
      creepparam = 0.
!     Hardening parameters
      hardeningparam = 0.
!     Irradiation parameters
      irradiationparam = 0.
!     Backstress parameters
      backstressparam = 0.
!
!
!     Interaction matrices
!     Initially set all to zero
      sintmat1=0.
      sintmat2=0.
      hintmat1=0.
      hintmat2=0.
!
!
!
!     Temperature in Celcius
      tcelsius = temperature - 273.15
!
!      
!
!
!
!
!
!
!
!     select crystal type
      select case(imat)
!
!     custom material - bcc (i.e. ferrite)
!     no temperature dependence
      case(1) 
!
!         Slip model
!         sinh law
          slipmodel = 1
!
!
!         Phase id
          iphase = 1
!
!         This can be 12 or 24
          nslip = 12
!
!
!         constant alpha
          slipparam(1) = 0.
!         constant beta
          slipparam(2) = 0.
!         psi - fraction of mobile dislocations
          slipparam(3) = 0.727d-2 
!         rhom0 - mobile dislocation density
          slipparam(4) = 0.035
!         DeltaF - activation energy
          slipparam(5) = 4.646312d-20
!         nu0 - attempt frequency
          slipparam(6) = 1.0d+11
!         gamma0 - multiplier for activation volume
!         =1/sqrt(Psi) in the ref. which was =1/sqrt(1.457e-4)
          slipparam(7) = 0.
!         Activation volume (factor of Burgers vector^3)
          slipparam(8) = 1.
!         Geometric factor (0.1-0.5)
          gf = 0.25
!
!         cubic slip system flag
!         0: inactive / 1: active
          cubicslip = 0 
!
!
!         Screw systems
          nscrew = 4
!
!
!
!         Initialize arrays
          burgerv = 0.
          tauc_0 = 0.
          rho_0 = 0.
!
!         Initial value of SSD density
          rho_0(1:nslip) = 0.
!
!         Initial value of forest density (hardening model-4)
          rhofor_0 = 0.
!
!         Initial value of substructure density (hardening model-4)
          rhosub_0 = 0.
!
!         crss [MPa]
          xtauc1 = 60.
!          xtauc1 = 45.
!
!
!         Burgers vectors [micrometers]
          burger1 = 2.48d-4
!
!
!         thermal expansion coefficients
          alpha1 = 0.
          alpha2 = 0.
          alpha3 = 0.
!
!
!         Example elastic modulus values
!     **********************************************************
!
!         steel-ferrum
!         Cristian Teodosiu, "Elastic Models of crystal Defects" 1982
!         C11=230.1d3
!         C12=134.6d3
!         C44=116.6d3
!
!         steel-ferrite
!         G.V. Kurdjumov, A.G. Khachaturyan, "Nature of axial ratio anomalies
!         of the martensite lattice and mechanism of diffusionless transformation"
!         page 1087
!         C11=233.5d3
!         C12=135.5d3
!         C44=118.0d3
!
!     **********************************************************
!
!         Cubic elastic constants [MPa]
!         Value used for ferrite grains
          C11=233.5d3 
          C12=135.5d3
          C44=118.0d3
!
!         re-calculate elastic constants
          e1 = (C11**2 + C11*C12 - 2.*C12**2)/(C11 + C12)
          v12 = C12/(C11 + C12)
          g12 = C44
!
!
!
!         assign the values based on crystal symmetry          
          e2 = e1
          e3 = e1
          v13 = v12
          v23 = v12
          g13 = g12
          g23 = g12
!
!
!         assign burgers vector scalars
          burgerv(1:nslip) = burger1
!
!         assign crss: same for all slip systems
          tauc_0(1:nslip) = xtauc1
!
!         c/a ratio for hcp crystals
!         dummy output
          caratio=1.
!
!
!         creep model
          creepmodel = 0      
!
!         hardening model
          hardeningmodel = 0  
!
!         irradiation model
          irradiationmodel = 0
!
!
!
!
!     custom material - fcc (i.e. copper)
!     no temperature dependence
      case(2)
!
!         Phase id
          iphase = 2
!
!
!
!         Slip model
!         power law
          slipmodel = 3
!
!         copper
!         Slip model parameters
          slipparam(1) = 1.0d-3
!         rate sensitivity exponent
!          slipparam(2) = 83.333
          slipparam(2) = 20.
!
!!         Inverse slip test parameters
!!         Slip model parameters
!          slipparam(1) = 1.0d-9
!!         rate sensitivity exponent
!          slipparam(2) = 13.         
!
!!         steel
!!         Slip model parameters
!          slipparam(1) = 1.0d-3
!!         rate sensitivity exponent
!          slipparam(2) = 7.143
!
!         Geometric factor (0.1-0.5)
          gf = 0.25
!
!
!         cubic slip system flag
!         0: inactive / 1: active
          cubicslip = 0            
!
!
!         This is 12 for fcc
          nslip = 12
!
!         Screw systems
          nscrew = 6
!
!
!
!         Initialize arrays
          burgerv = 0.
          tauc_0 = 0.
          rho_0 = 0.
!
!         Initial value of SSD density
          rho_0(1:nslip) = 0.
!
!
!         Initial value of forest density (hardening model-4)
          rhofor_0 = 0.
!
!         Initial value of substructure density (hardening model-4)
          rhosub_0 = 0.
!
!         crss [MPa]
!
!!         steel
!          xtauc1 = 32.
!
!         Copper
          xtauc1 = 16.
!
!!         Inverse slip test parameter
!          xtauc1 = 32.
!
!         Burgers vectors [micrometers]
          burger1 = 2.56d-4
!
!
!         thermal expansion coefficients
          alpha1 = 0.
          alpha2 = 0.
          alpha3 = 0.
!
!     Example elastic modulus values
!     **********************************************************
!
!         copper
!         "Texture and Anisotropy", 
!         Cambridge University Press, 1998,  page 300
!         C11=168.0d3
!         C12=121.4d3
!         C44=75.4d3
!
!         "The Mechanics of Crystals and Textured Polycrystals", 
!         Oxford University Press, 1993, page 16
!         C11=166.1d3
!         C12=199.0d3 wrong! correct value: 119.0d3
!         C44=75.6d3
!
!         H.P.R. Frederikse, "Hanbook of Chemistry and Physics" 1995, 12- 38
!         C11=168.3d3
!         C12=122.1d3
!         C44=75.7d3
!
!         aluminum
!         "The Mechanics of Crystals and Textured Polycrystals",
!         Oxford University Press, 1993, page 16 
!         C11=107.3d3
!         C12=60.9d3
!         C44=28.3d3
!
!         H.P.R. Frederikse, "Hanbook of Chemistry and Physics" 1995, 12- 38
!         C11=106.75d3
!         C12=60.41d3
!         C44=28.34d3
!
!         Cristian Teodosiu, "Elastic Models of crystal Defects" 1982
!         C11=106.43d3
!         C12=60.35d3
!         C44=28.21d3
!
!         steel-austenite
!         S. Turteltaub and A.S.J. Suiker, "Transformation-induced plasticit in ferrous alloys", 2005
!         page 1765, Eq. (37)
!         C11=268.5d3
!         C12=156.d3
!         C44=136.d3
!         
!         steel
!         C11=204.6d3
!         C12=137.7d3
!         C44=126.6d3          
!
!     **********************************************************          
!
!         Cubic elastic constants [MPa]
!         Value used for copper
          C11 = 170.d3
          C12 = 124.d3
          C44 = 75.d3
!
!         re-calculate elastic constants 
          e1 = (C11**2 + C11*C12 - 2.*C12**2)/(C11 + C12)
          v12 = C12/(C11 + C12)
          g12 = C44
!
!         assign the values based on crystal symmetry          
          e2 = e1
          e3 = e1
          v13 = v12
          v23 = v12
          g13 = g12
          g23 = g12
!
!
!         assign burgers vector scalars
          burgerv(1:nslip) = burger1

!         assign crss: same for all slip systems
          tauc_0(1:nslip) = xtauc1
!        
!         c/a ratio for hcp crystals
!         dummy output
          caratio = 1.
!
!         creep model
          creepmodel = 0
!
!
!
!         hardening model
          hardeningmodel = 1
!
!!     Kocks-Mecking hardening with substructure evolution
!!     Reference: https://doi.org/10.1016/j.actamat.2010.06.021
!!
!!         hardening model
!          hardeningmodel = 4
!!
!          hardeningparam = 0.
!!         k1 - forest hardening
!          hardeningparam(1) = 40.d-4
!!         k2 - forest annihilation
!          hardeningparam(2) = 1.
!!         q - substructure
!          hardeningparam(7) = 4.
!!         f - substructure
!          hardeningparam(8) = 20.
!!         ksub - substructure
!          hardeningparam(9) = 0.086
!
!
!!         Inverse slip test parameter
!          hardeningmodel = 0
!
!
!         copper          
!         Hardening rate - h0
          hardeningparam(1)=250.
!         Saturation strength for slip - ss
          hardeningparam(2)=190.
!         Hardening exponent - a
          hardeningparam(3)=2.5
!         Latent hardening coefficient - q
          hardeningparam(4)=1.4
!
!
!!         steel
!!         Hardening rate - h0
!          hardeningparam(1)=217.8
!!         Saturation strength for slip - ss
!          hardeningparam(2)=257.
!!         Hardening exponent - a
!          hardeningparam(3)=2.5
!!         Latent hardening coefficient - q
!          hardeningparam(4)=1.1          
!
!
!         irradiation model
          irradiationmodel = 0
!
!
!
!
!         Hardening interactions - latent hardening
          hintmat1 = hardeningparam(4)
          do k = 1, int(nslip/3.)       
	        do i = 1, 3
                  do j = 1, 3
	                hintmat1(3*(k-1)+i, 3*(k-1)+j)=1.
                  enddo
              enddo
          enddo
!!
!         Backstress parameter
          backstressparam(1) = 0.25
!
!     custom material - hcp (i.e. zirconium)
!     no temperature dependence
      case(3) 
!
!         Phase id
          iphase = 3
!
!
!         Slip model
!         sinh law
          slipmodel = 1
!
!         Slip model parameters
!         constant alpha
          slipparam(1) = 0.
!         constant beta
          slipparam(2) = 0.
!         fraction of mobile dislocations
          slipparam(3) = 1.
!         reference mobile dislocations (1/micrometer^2)
          slipparam(4) = 0.01
!         activation energy for slip (J)
          slipparam(5) = 5.127d-20
!         attempt frequency
          slipparam(6) = 1.d11
!         scaling for jump distance
          slipparam(7) = 1.
!         activation volume
          slipparam(8) = 20.93
!
!
!
!         Geometric factor (0.1-0.5)
          gf = 0.22364 
!
!         cubic slip system flag
!         0: inactive / 1: active
          cubicslip = 0          
!
!
!
!         number of slip systems
          nslip = 30
!
!         Screw systems
          nscrew = 9
!
!
!         Initialize arrays
          burgerv = 0.
          tauc_0 = 0.
          rho_0 = 0.
!
!         Initial value of SSD density
          rho_0(1:nslip) = 10.
!
!
!         Initial value of forest density (hardening model-4)
          rhofor_0 = 0.
!
!         Initial value of substructure density (hardening model-4)
          rhosub_0 = 0.
!
!         Burgers vectors [micrometers]
          burger1 = 3.2d-4   
          burger2 = 6.0d-4 ! sqrt(a^2 + c^2)                      
!
!
!         elastic constants [MPa]
          e1 = 98.32d3
          e3 = 123.28d3
          g12 = 32.01d3
          v12 = 0.40
          v13 = 0.24
!
!         crss [MPa]
          xtauc1 = 187.7 ! basal
          xtauc2 = 140.8 ! prismatic
          xtauc3 = 140.8 ! pyramidal
          xtauc4 = 489.8 ! pyramidal-1
          xtauc5 = 2449.1 ! pyramidal-2
!
!         thermal expansion coefficients [1/K]
          alpha1 = 0.
          alpha2 = 0.
          alpha3 = 0.
!
!
!
!         elastic constants based on crystal symmetry
          e2 = e1
          g13 = e3/(1.+v13)/2.
          g23 = g13
          v23 = v13
!
!         assign Burgers vector scalars
          burgerv(1:6) = burger1
          burgerv(7:30) = burger2
!
!         assign crss
          tauc_0(1:3) = xtauc1
          tauc_0(4:6) = xtauc2
          tauc_0(7:12) = xtauc3
          tauc_0(13:24) = xtauc4
          tauc_0(25:30) = xtauc5
!
!         c/a ratio          
          caratio = 1.57
!
!         creep model
          creepmodel = 0
!
!         hardening model
!         linear hardening
          hardeningmodel = 2
!         hardening parameter (1/micrometer^2)
          hardeningparam(1) = 2600.
!
!         irradiation model
          irradiationmodel = 2
!         Irradiation parameters
!         Number of different types of defects
          irradiationparam(1) = 3.
!         Number density of defects (1/micrometer^3)
          irradiationparam(2:4) = 2.033d+4
!         size of defects (micrometer)
          irradiationparam(5:7) = 1.4d-3
!         defect vectors (crystallographic directions)
          irradiationparam(8:10) = (/ 4.0, 6.0, 11.0 /) 
!         Strength interaction matrix coefficients (two independent parameters)
!         when a=0 (a: reaction segment)
          irradiationparam(11) = 1.25
!         when a not equals 0
          irradiationparam(12) = 1.875
!         Hardening (Softening) interaction matrix coefficients (two independent parameters)
!         when a=0 (a: reaction segment)
          irradiationparam(13) = 0.5
!         when a not equals 0
          irradiationparam(14) = 0.          
!
!
!
!
!
!     tungsten - bcc
      case(4) 
!
!         Phase id
          iphase = 1      
!
!
!         Slip model
!         sinh law
          slipmodel = 1
!
!         Burgers vectors [micrometers]
          burger1 = 2.74d-4
!
!         Slip model parameters
!         Partly available in the reference https://doi.org/10.1016/j.ijplas.2018.05.001
!         constant alpha
          slipparam(1) = 0.
!         constant beta
          slipparam(2) = 0.0159
!         psi - fraction of mobile dislocations
          slipparam(3) = 0.727d-2    
!         rhom0 - mobile dislocation density
          slipparam(4) = 0.035
!         DeltaF - activation energy to overcome Pierls barrier
          slipparam(5) = 3.524788e-20
!         nu0 - attempt frequency
          slipparam(6) = 1.0d11
!         gamma0 - multiplier for activation volume 
!         =1/sqrt(Psi) in the ref. which was =1/sqrt(1.457e-4)
          slipparam(7) = 0.
!         AV0
          slipparam(8) = 0.
!
!         Geometric factor (0.1-0.5)
          gf = 0.1        
!
!
!         cubic slip system flag
!         0: inactive / 1: active
          cubicslip = 0             
!
!
!         number of slip systems
          nslip = 12
!
!         Screw systems
          nscrew = 4
!
!
!         Initialize arrays
          burgerv = 0.
          tauc_0 = 0.
          rho_0 = 0.
!         Initial value of SSD density
          rho_0(1:nslip) = 0.     
!
!
!         Initial value of forest density (hardening model-4)
          rhofor_0 = 0.
!
!         Initial value of substructure density (hardening model-4)
          rhosub_0 = 0.
!
!         crss [MPa]
          xtauc1 = 360.0 ! 900 MPa in the reference
!
!         elastic constants [MPa]
          e1 = 421d3
          g12 = 164.4d3
          v12 = 0.28
!
!         thermal expansion coefficients
          alpha1 = 9.5d-6
          alpha2 = alpha1
          alpha3 = 0.5895*alpha1
!
!
!
!
!         elastic constants based on crystal symmetry
          e2 = e1
          e3 = e1
          v13 = v12
          v23 = v12
          g13 = g12
          g23 = g12
!
!         assign burgers vector scalars
          burgerv(1:nslip) = burger1
!
!         assign crss: same for all slip systems
          tauc_0(1:nslip) = xtauc1
!
!         c/a ratio for hcp crystals
!         dummy output
          caratio=1.
!
!         creep model
          creepmodel = 0
!
!         hardening model
          hardeningmodel = 2
          hardeningparam(1) = 1000.
!
!         irradiation model
          irradiationmodel = 1
!         irradiation parameters
!         initial strength
          irradiationparam(1) = 750.
!         solute strength saturation strain
          irradiationparam(2) = 0.025
!         factor for mobile density
          irradiationparam(3) = 3.457d-2
!
!
!
!
!     copper - fcc
      case(5) 
!
!         Phase id
          iphase = 2
!
!         Slip model
!         sinh law
          slipmodel = 1
!
!         Slip model parameters
!         constant alpha
          slipparam(1) = 0.02
!         constant beta
          slipparam(2) = 0.1
!
!         Geometric factor (0.1-0.5)
          gf = 0.25
!
!         cubic slip system flag
!         0: inactive / 1: active
          cubicslip = 0
!
!         number of slip systems
          nslip = 12
!
!         Screw systems
          nscrew = 6
!
!
!         Initialize arrays
          burgerv = 0.
          tauc_0 = 0.
          rho_0 = 0.
!
!         Initial value of SSD density
          rho_0(1:nslip) = 0.
!
!
!         Initial value of forest density (hardening model-4)
          rhofor_0 = 0.
!
!         Initial value of substructure density (hardening model-4)
          rhosub_0 = 0.
!
!         crss [MPa]
          xtauc1 = 20.0
!
!         assign crss: same for all slip systems
          tauc_0(1:nslip) = xtauc1
!
!         Burgers vectors [micrometers]
          burger1 = 2.55d-4
!
!         elastic constants [MPa]
          e1 = 66.69d3
          v12 = 0.4189
          g12 = 75.4d3
!
!
!         elastic constants based on crystal symmetry
          e2 = e1
          e3 = e1
          v13 = v12
          v23 = v12
          g13 = g12
          g23 = g12          
!
!
!         thermal expansion coefficients
          alpha1 = 13.0d-6
          alpha2 = alpha1
          alpha3 = alpha1
!
!
          burgerv(1:nslip)=burger1
!
!         c/a ratio for hcp crystals
!         dummy output
          caratio=1.
!
!
!         creep model
          creepmodel = 0
!
!
!         hardening model
          hardeningmodel = 0
!
!
!         irradiation model
          irradiationmodel = 0
!
!
!
!     carbide - fcc
      case(6) 
!
!         Phase id
          iphase = 2
!
!         Slip model
!         power law
          slipmodel = 3
!
!         Slip model parameters
!         Reference strain rate
          slipparam(1) = 1.0d-3
!         Rate sensitivity exponent
          slipparam(2) = 50
!
!         Geometric factor (0.1-0.5)
          gf = 0.25  
!
!         cubic slip system flag
!         0: inactive / 1: active
          cubicslip = 0             
!
!         number of slip systems
          nslip = 12
!
!         Screw systems
          nscrew = 0
!
!         Initialize arrays
          burgerv = 0.
          tauc_0 = 0.
          rho_0 = 0.
!
!         Initial value of SSD density
          rho_0(1:nslip) = 0.
!
!         Initial value of forest density (hardening model-4)
          rhofor_0 = 0.
!
!         Initial value of substructure density (hardening model-4)
          rhosub_0 = 0.
!
!         crss [MPa]
          xtauc1 = 2300.0
!
!         Burgers vectors [micrometers]
          burger1 = 3.5072d-4
!
!         elastic constants [MPa]
          e1 = 207.0d4
          v12 = 0.28
!
!         thermal expansion coefficients
          alpha1 = 4.5d-6
          alpha2 = alpha1
          alpha3 = alpha1
!
!
!         elastic constants based on crystal symmetry      
          e3 = e1
          e2 = e1
          v13 = v12
          v23 = v12
          g12 = e1/(2.0*(1.0+v12)) 
          g13 = g12
          g23 = g12
!
!         assign Burgers vector scalars
          burgerv(1:nslip) = burger1
!
!         assign crss: same for all slip systems
          tauc_0(1:nslip) = xtauc1
!
!         c/a ratio for hcp crystals
!         dummy output
          caratio=1.
!
!         creep model
          creepmodel = 0
!
!         hardening model
          hardeningmodel = 0
!
!         irradiation model
          irradiationmodel = 0
!
!
!
!
!
!
!     CMSX-4 - fcc + cubic
      case(7)
!
!         Phase id
          iphase = 2
!
!         Slip model
!         Double exponent law
          slipmodel = 2
!
!         Slip model parameters
!         Reference strain rate
          slipparam(1) = 1.0d7 
!         Inner exponent
          slipparam(2) = 0.78
!         Outer exponent
          slipparam(3) = 1.15
!         Activation energy for octahedral slip  (J)
          slipparam(4) = 9.39d-19
!         Activation energy for cubic slip  (J)
          slipparam(5) = 1.17d-18
!
!         Geometric factor (0.1-0.5)
          gf = 0.25          
!
!         cubic slip system flag
!         0: inactive / 1: active
          cubicslip = 0
!              
!         Screw systems
          nscrew = 6
!
!
!
!         number of slip systems
          if (cubicslip == 0) then
              nslip = 12
          elseif (cubicslip == 1) then
              nslip = 18
          else
              call error(3)
          end if
!
!         Initialize arrays
          burgerv = 0.
          tauc_0 = 0.
          rho_0 = 0.
!
!         Initial value of SSD density
          rho_0(1:nslip) = 0.
!
!
!         Initial value of forest density (hardening model-4)
          rhofor_0 = 0.
!
!         Initial value of substructure density (hardening model-4)
          rhosub_0 = 0.
!         
!
!         crss [MPa]
          if (tcelsius <= 850.0) then   !  temperature in celsius
!
              tauc_0=-0.000000001051*tcelsius**4 +
     + 0.000001644382*tcelsius**3 -
     + 0.000738679333*tcelsius**2 + 0.128385617901*tcelsius +
     + 446.547978926622
!
              if (cubicslip == 1) then
!
                  tauc_0(13:18)=-0.000000001077*tcelsius**4 +
     + 0.000001567820*tcelsius**3 -
     + 0.000686532147*tcelsius**2 -
     + 0.074981918833*tcelsius +
     + 571.706771689334
!
              end if
!
          else
!
              tauc_0=-1.1707*tcelsius + 1478.9
!
              if (cubicslip == 1) then
!
                  tauc_0(13:18)=-0.9097*tcelsius + 1183
!
              end if
!
          end if ! end temperature check
!
!         tcelsius dependent stiffness constants [MPa]
          if (tcelsius <= 800.0) then   !  celsius units
!
              c11=-40.841*tcelsius+251300
              c12=-14.269*tcelsius+160965
!
          else
!
              c11=0.111364*tcelsius**2-295.136*tcelsius+382827.0
              c12=-0.000375*tcelsius**3+1.3375*tcelsius**2 -
     + 1537.5*tcelsius+716000
!
          end if ! end temperature check  
!
          g12=-0.00002066*tcelsius**3+0.021718*tcelsius**2 -
     + 38.3179*tcelsius+129864
!
          e1 = (c11-c12)*(c11+2*c12)/(c11+c12)
          v12 = e1*c12/((c11-c12)*(c11+2*c12))
!
!         temperature dependent thermal expansion coefficient
          alpha1 = 9.119d-9*tcelsius +1.0975d-5
!
!         Eralp: Burgers vector was not defined for this
          burger1=2.56d-4
!
!         elastic constants based on crystal symmetry
          e2 = e1
          e3 = e1
          v13 = v12
          v23 = v12
          g13 = g12
          g23 = g12
!
!         assign Burgers vector scalars
          burgerv(1:nslip) = burger1         
!
!         c/a ratio for hcp crystals
!         dummy output
          caratio=1.
!
!         creep model
          creepmodel = 1
!
!         hardening model
          hardeningmodel = 0
!
!         creep parameters
!         reference rate for creep (1/s)
          creepparam(1) = 4.0d+8
!         stress multiplier for creep (1/MPa)
          creepparam(2) = 3.2d-2
!         activation energy for creep (J/mol)
          creepparam(3) = 460000.
!         reference rate for damage (1/s)
          creepparam(4) = 6.0d6
!         stress multiplier for damage (1/MPa)
          creepparam(5) = -5.0d-8
!         activation energy for damage (J/mol)
          creepparam(6) = 340000.
!
!         irradiation model
          irradiationmodel = 0
!
!
!     zirconium - hcp
      case(8)
!
!         Phase id
          iphase = 3
!
!         Slip model
!         sinh law
          slipmodel = 1
!
!         Slip model parameters
!         constant alpha
          slipparam(1) = 0.1
!         constant beta
          slipparam(2) = 0.1
!
!         Geometric factor (0.1-0.5)
          gf = 0.25          
!
!         cubic slip system flag
!         0: inactive / 1: active
          cubicslip = 0
!
!         number of slip systems
          nslip = 12
!
!         Screw systems
          nscrew = 9
!
!
!         Initialize arrays
          burgerv = 0.
          tauc_0 = 0.
          rho_0 = 0.
!
!         Initial value of SSD density
          rho_0(1:nslip) = 0.
!
!         Initial value of forest density (hardening model-4)
          rhofor_0 = 0.
!
!         Initial value of substructure density (hardening model-4)
          rhosub_0 = 0.
!
!
!         Burgers vectors [micrometers]
          burger1 = 2.28d-4   
          burger2 = 4.242d-4 ! sqrt(a^2 + c^2)                      
!
!
!         elastic constants [MPa]
          e1 = 289.38d3
          e3 = 335.17d3
          g12 = 132.80d3
          g13 = 162.50d3
          v12 = 0.09 
          v13 = 0.04
!
!         crss [MPa]
          xtauc1 = 15.2 ! basal
          xtauc2 = 67.7 ! prismatic
          xtauc4 = 2000.0 ! pyramidal
!
!         thermal expansion coefficients [1/K]
          alpha1 = 9.5d-6
          alpha2 = alpha1
          alpha3 = 0.5895*alpha1
!
!
!
!         elastic constants based on crystal symmetry
          e2 = e1
          g23 = g13
          v23 = v13
!
!         assign Burgers vector scalars
          burgerv(1:6) = burger1
          burgerv(7:12) = burger2
!
!         assign crss
          tauc_0(1:3) = xtauc1
          tauc_0(4:6) = xtauc2
          tauc_0(7:12) = xtauc4              
!
!         c/a ratio
          caratio = 1.57
!
!         creep model
          creepmodel = 0
!
!         hardening model
          hardeningmodel = 0
!
!         irradiation model
          irradiationmodel = 0
!
!
!
!
!     Material constants by Alvaro Martinez Pechero
!     Berilyum - hcp
      case(9)           
!
!         Phase id
          iphase = 3
!
!
!         Slip model
!         sinh law
          slipmodel = 1
!
!         Slip model parameters
!         constant alpha
          slipparam(1) = 0.1
!         constant beta
          slipparam(2) = 0.1
!
!         Geometric factor (0.1-0.5)
          gf = 0.25
!
!         cubic slip system flag
!         0: inactive / 1: active
          cubicslip = 0
!
!
!
!         number of slip systems
          nslip = 12
!
!         Screw systems
          nscrew = 3
!
!
!         Initialize arrays
          burgerv = 0.
          tauc_0 = 0.
          rho_0 = 0.
!
!         Initial value of SSD density
          rho_0(1:nslip) = 1.
!
!         Initial value of forest density (hardening model-4)
          rhofor_0 = 0.
!
!
!         Initial value of substructure density (hardening model-4)
          rhosub_0 = 0.
!
!         Burgers vectors [micrometers]
          burger1 = 2.28d-4
          burger2 = 4.242d-4 ! sqrt(a^2 + c^2)
!
!
!         elastic constants [MPa]
          e1 = 289.38d3
          e3 = 335.17d3
          g12 = 132.80d3
          g13 = 162.50d3
          v12 = 0.09
          v13 = 0.04
!
!         crss [MPa]
          xtauc1 = 20.2 ! basal
          xtauc2 = 88.7 ! prismatic
          xtauc4 = 188. ! pyramidal
!
!         thermal expansion coefficients [1/K]
          alpha1 = 0.
          alpha2 = 0.
          alpha3 = 0.
!
!
!
!         elastic constants based on crystal symmetry
          e2 = e1
          g23 = g13
          v23 = v13
!
!         assign Burgers vector scalars
          burgerv(1:6) = burger1
          burgerv(7:12) = burger1
!
!         assign crss
          tauc_0(1:3) = xtauc1
          tauc_0(4:6) = xtauc2
          tauc_0(7:12) = xtauc4
!
!         c/a ratio
          caratio = 1.57
!
!         creep model
          creepmodel = 0
!
!         hardening model
          hardeningmodel = 0
!
!         irradiation model
          irradiationmodel = 0
!
!
!     alpha-uranium - alphauranium
      case(10) 
!
!         Phase id
          iphase = 4
!
!         Slip model
!         power law
          slipmodel = 3
!
!         Slip model parameters
!         reference strain rate
          slipparam(1) = 1.0d-3
!         rate sensitivity exponent
          slipparam(2) = 20.
!
!         Geometric factor (0.1-0.5)
          gf = 0.25          
!
!         cubic slip system flag
!         0: inactive / 1: active
          cubicslip = 0          
!
!         number of slip systems
          nslip = 8
!
!         Screw systems
          nscrew = 0
!
!         Initialize arrays
          burgerv = 0.
          tauc_0 = 0.
          rho_0 = 0.
!
!         Initial value of SSD density
          rho_0(1:nslip) = 0. 
!
!         Initial value of forest density (hardening model-4)
          rhofor_0 = 0.
!
!         Initial value of substructure density (hardening model-4)
          rhosub_0 = 0.
!
!         Burgers vectors [micrometers]
          burgerv(1:2) = 2.85d-4
          burgerv(3:4) = 6.51d-4
          burgerv(5:8) = 11.85d-4
!
!         constant factor tau_0^alpha for crss [MPa]
!         calhoun 2013 values
!         with temperature dependence as in zecevic 2016
!         added after hardening is included
          tauc_0(1) = 24.5
          tauc_0(2) = 85.5
          tauc_0(3) = 166.5
          tauc_0(4) = 166.5
          tauc_0(5) = 235.0
          tauc_0(6) = 235.0
          tauc_0(7) = 235.0
          tauc_0(8) = 235.0
!
!
!
!         elastic moduli [MPa] and poissons ratios
!         see prs literature review by philip earp
!         short crack propagation in uranium, an anisotropic polycrystalline metal
!         temperature dependence according to daniel 1971
          e1 = 203665.987780 * (1.0 - 0.000935*(temperature-293.0))
          e2 = 148588.410104 * (1.0 - 0.000935*(temperature-293.0))
          e3 = 208768.267223 * (1.0 - 0.000935*(temperature-293.0))
          v12 = 0.242363
          v13 = -0.016293
          v23 = 0.387816
          g12 = 74349.442379 * (1.0 - 0.000935*(temperature-293.0))
          g13 = 73421.439060 * (1.0 - 0.000935*(temperature-293.0))
          g23 = 124378.109453 * (1.0 - 0.000935*(temperature-293.0))
!
!         define thermal expansion coefficients as a function of temperature
!         lloyd, barrett, 1966
!         thermal expansion of alpha uranium
!         journal of nuclear materials 18 (1966) 55-59
          alpha1 = 24.22d-6 - 9.83d-9 * temperature + 
     + 46.02d-12 * temperature * temperature
          alpha2 = 3.07d-6 + 3.47d-9 * temperature -
     + 38.45d-12 * temperature * temperature
          alpha3 = 8.72d-6 + 37.04d-9 * temperature +
     + 9.08d-12 * temperature * temperature
!
!
!         c/a ratio for hcp crystals
!         dummy output
          caratio=1.
!
!         creep model
          creepmodel = 0
!
!         hardening model
          hardeningmodel = 0
!
!         irradiation model
          irradiationmodel = 0   
!
!     UKAEA - Vikram Phalke
!     copper - fcc
!     no temperature dependence
      case(11)
!
!         Phase id
          iphase = 2
!
!
!
!         Slip model
!         power law
          slipmodel = 1
!
!         copper
!         Slip model parameters
          slipparam(1) = 2.0d-5
!         rate sensitivity exponent
          slipparam(2) = 1.
!
!
!         Geometric factor (0.1-0.5)
          gf = 0.2
!
!
!         cubic slip system flag
!         0: inactive / 1: active
          cubicslip = 0            
!
!
!         This is 12 for fcc
          nslip = 12
!
!         Screw systems
          nscrew = 6
!
!
!
!         Initialize arrays
          burgerv = 0.
          tauc_0 = 0.
          rho_0 = 0.
!
!         Initial value of SSD density
          rho_0(1:nslip) = 14.98
!
!
!         Initial value of forest density (hardening model-4)
          rhofor_0 = 0.
!
!         Initial value of substructure density (hardening model-4)
          rhosub_0 = 0.
!
!         crss [MPa]
!
!         Copper
          xtauc1 = 1.
!
!
!         Burgers vectors [micrometers]
          burger1 = 2.56d-4
!
!
!         thermal expansion coefficients
          alpha1 = 0.
          alpha2 = 0.
          alpha3 = 0.
!      
!
!         Cubic elastic constants [MPa]
!         Value used for copper
          C11 = 170.d3
          C12 = 124.d3
          C44 = 75.d3
!
!         re-calculate elastic constants 
          e1 = (C11**2 + C11*C12 - 2.*C12**2)/(C11 + C12)
          v12 = C12/(C11 + C12)
          g12 = C44
!
!         assign the values based on crystal symmetry          
          e2 = e1
          e3 = e1
          v13 = v12
          v23 = v12
          g13 = g12
          g23 = g12
!
!
!         assign burgers vector scalars
          burgerv(1:nslip) = burger1

!         assign crss: same for all slip systems
          tauc_0(1:nslip) = xtauc1
!        
!         c/a ratio for hcp crystals
!         dummy output
          caratio = 1.
!
!         creep model
          creepmodel = 0
!
!
!
!         hardening model
          hardeningmodel = 5
!
!
!
!
!
!         copper          
!         Hardening rate - k1
          hardeningparam(1)=0.06
!         Softening rate - k2
          hardeningparam(2)=35.
!         Latent hardening coefficient - q1
          hardeningparam(3)=1.35
!         Latent hardening coefficient - q2
          hardeningparam(4)=1.2
!
!
!!         hardening model
!          hardeningmodel = 3
!!
!          hardeningparam = 0.
!!         k1 - forest hardening
!          hardeningparam(1) = 40.d-4
!!         k2 - forest annihilation
!          hardeningparam(2) = 1.
!
!
!         irradiation model
          irradiationmodel = 0
!
!
!         Define Kronecker delta
          kdelta = 0.
          do is=1,nslip
              kdelta(is,is)=1.
          end do
!    
          q1=hardeningparam(3)
          q2=hardeningparam(4)
          hintmat1=0.
!         Define the hardening interaction matrix
          do is=1,nslip
              do js=1,nslip
                  hintmat1(is,js) = q1 + (1.-q2)*kdelta(is,js)
              end do
          end do
!
!
!
!         Backstress parameter
          backstressparam(1) = 0.
!
!
!     UKAEA - Vikram Phalke
!     CuCrZr - fcc
!     no temperature dependence
      case(12)
!
!         Phase id
          iphase = 2
!
!
!
!         Slip model
!         power law
          slipmodel = 1
!
!         copper
!         Slip model parameters
          slipparam(1) = 2.0d-5
!         rate sensitivity exponent
          slipparam(2) = 1.
!
!
!         Geometric factor (0.1-0.5)
          gf = 0.2
!
!
!         cubic slip system flag
!         0: inactive / 1: active
          cubicslip = 0            
!
!
!         This is 12 for fcc
          nslip = 12
!
!         Screw systems
          nscrew = 6
!
!
!
!         Initialize arrays
          burgerv = 0.
          tauc_0 = 0.
          rho_0 = 0.
!
!         Initial value of SSD density [1/micrometer^2]
          rho_0(1:nslip) = 14.98
!
!
!         Initial value of forest density (hardening model-4)
          rhofor_0 = 0.
!
!         Initial value of substructure density (hardening model-4)
          rhosub_0 = 0.
!
!         crss [MPa]
!
!         CuCrZr
          xtauc1 = 84.6
!
!
!         Burgers vectors [micrometer]
          burger1 = 2.56d-4
!
!
!         thermal expansion coefficients
          alpha1 = 0.
          alpha2 = 0.
          alpha3 = 0.
!      
!
!         Cubic elastic constants [MPa]
!         Value used for copper
          C11 = 170.d3
          C12 = 124.d3
          C44 = 75.d3
!
!         re-calculate elastic constants 
          e1 = (C11**2 + C11*C12 - 2.*C12**2)/(C11 + C12)
          v12 = C12/(C11 + C12)
          g12 = C44
!
!         assign the values based on crystal symmetry          
          e2 = e1
          e3 = e1
          v13 = v12
          v23 = v12
          g13 = g12
          g23 = g12
!
!
!         assign burgers vector scalars
          burgerv(1:nslip) = burger1

!         assign crss: same for all slip systems
          tauc_0(1:nslip) = xtauc1
!        
!         c/a ratio for hcp crystals
!         dummy output
          caratio = 1.
!
!         creep model
          creepmodel = 0
!
!
!
!         hardening model (KM with precipitates)
          hardeningmodel = 5
!
!
!
!
!
!         CuCrZr  
!         Hardening rate - k1
          hardeningparam(1)=0.06
!         Softening rate - k2
          hardeningparam(2)=35.
!         Latent hardening coefficient - q1
          hardeningparam(3)=1.35
!         Latent hardening coefficient - q2
          hardeningparam(4)=1.2
!
!         Parameters related to precipitate strengthening
!         Geometric/Strength factor for precipitates
          hardeningparam(5)=0.08
!         Number density of precipitates [1/micrometer^3]
          hardeningparam(6)=1.76d6
!         Particle diameter [micrometer]
          hardeningparam(7)=3.2d-3
!
!!         hardening model
!          hardeningmodel = 3
!!
!          hardeningparam = 0.
!!         k1 - forest hardening
!          hardeningparam(1) = 40.d-4
!!         k2 - forest annihilation
!          hardeningparam(2) = 1.
!
!
!         irradiation model
          irradiationmodel = 0
!
!
!         Define Kronecker delta
          kdelta = 0.
          do is=1,nslip
              kdelta(is,is)=1.
          end do
!    
          q1=hardeningparam(3)
          q2=hardeningparam(4)
          hintmat1=0.
!         Define the hardening interaction matrix
          do is=1,nslip
              do js=1,nslip
                  hintmat1(is,js) = q1 + (1.-q2)*kdelta(is,js)
              end do
          end do
!
!
!
!         Backstress parameter
          backstressparam(1) = 0.
!
!
!
      case default
!
          call error(1)
!
      end select
!
!     *** set up elastic stiffness matrix in lattice system ***  
!     http://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php#Sect3_2_11
!     however, the voigt notation in abaqus for strain and stress vectors
!     has the following order:
!     stress = (sigma11,sigma22,sigma33,tau12 tau13 tau23)
!     strain = (epsilon11,epsilon22,epsilon33,gamma12,gamma13,gamma23)
!     Calculate the remaining Poisson's ratios
      v21 = e2/e1*v12
      v31 = e3/e1*v13
      v32 = e3/e2*v23
!
!     constant as a factor
      cst = 1./(1.-v12*v21-v23*v32-v31*v13
     + -2.*v21*v32*v13)
!
      Cc=0.
      Cc(1,1) = e1*(1.-v23*v32)*cst
      Cc(2,2) = e2*(1.-v13*v31)*cst
      Cc(3,3) = e3*(1.-v12*v21)*cst
      Cc(1,2) = e1*(v21+v31*v23)*cst
      Cc(1,3) = e1*(v31+v21*v32)*cst
      Cc(2,3) = e2*(v32+v12*v31)*cst
      Cc(4,4) = g12
      Cc(5,5) = g13
      Cc(6,6) = g23
!
!     symmetrize compliance matrix
      Cc(2,1) =  Cc(1,2)
      Cc(3,1) =  Cc(1,3)
      Cc(3,2) =  Cc(2,3)
!
!
!
!     define thermal eigenstrain in the lattice system
      alphamat=0.
      alphamat(1,1) = alpha1 
      alphamat(2,2) = alpha2 
      alphamat(3,3) = alpha3
!
!
!
!
      return
!
      end subroutine materialparam
!
!    
!
!
!
      end module usermaterials