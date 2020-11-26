********************************************
** KMATERIAL sets the material constants  **
**      for elasticity and plasticity     **
********************************************

      SUBROUTINE kMaterialParam(iphase,caratio,compliance,G12,thermat,
     + gammast,burgerv,nSys,tauc,screwplanes,Temperature,
     + tauctwin,nTwin,twinon,nTwinStart,nTwinEnd)

      ! crystal type
      INTEGER,intent(in) :: iphase

      ! number of slip systems
      INTEGER,intent(in) :: nSys

      ! current temperature
      REAL*8,intent(in) :: Temperature

      ! total number of twin systems
      INTEGER,intent(in) :: nTwin

      ! twin systems activation flag
      INTEGER,intent(in) :: twinon

      ! the active twins are the ones in the
      ! interval [nTwinStart,nTwinEnd] in the
      ! twin system file
      INTEGER,intent(in) :: nTwinStart,nTwinEnd

      ! c/a ratio for hcp crystals
      REAL*8,intent(out) :: caratio

      ! elastic compliance matrix in the crystal reference frame 
      REAL*8,intent(out) :: compliance(6,6)

      ! shear modulus for Taylor dislocation law
      REAL*8,intent(out) :: G12

      ! thermal eigenstrain to model thermal expansion
      REAL*8,intent(out) :: thermat(3,3)

      ! prefactor for SSD evolution equation
      REAL*8,intent(out) :: gammast 

      ! Burgers vectors
      REAL*8,intent(out) :: burgerv(nSys)

      ! critical resolved shear stress of slip systems
      REAL*8,intent(out) :: tauc(nSys)

      ! number of screw planes 
      integer,intent(out) :: screwplanes

      ! critical resolved shear stress of twin systems
      REAL*8,intent(out) :: tauctwin(nTwin)

******************************************
** The following parameters must be set **

      ! material name
      ! more materials can be added
      ! by modifying this subroutine
      ! materials available are the following:
      ! 'zirconium', 'alphauranium', 'tungsten', 'copper', 'carbide',
      ! 'olivine'
      character(len=*), parameter :: matname = 'zirconium' 

**       End of parameters to set       **
******************************************

      ! Burgers vector scalars
      REAL*8 :: burger1, burger2

      ! Elastic constants scalars
      REAL*8 :: E1, E2, E3, G13, v12, v13
      REAL*8 :: G23, v23

      ! hcp: basal critical resolved shear stress
      REAL*8 :: XTAUC1

      ! hcp: prismatic critical resolved shear stress
      REAL*8 :: XTAUC2

      ! hcp: prismatic critical resolved shear stress
      REAL*8 :: XTAUC4

      ! thermal expansion coefficients
      REAL*8 :: alpha1, alpha2, alpha3

      INTEGER :: i, j

      ! select crystal type       
      SELECT CASE(iphase)

      case(0) !hcp

        SELECT CASE(matname)

        case('zirconium')

        ! Burgers vectors (um)
        burger1 = 2.28E-4   
        burger2 = 4.242E-4 ! sqrt(a^2 + c^2)                      
        caratio = 1.57

        ! elastic constants (MPa units)
        E1 = 289.38E3
        E3 = 335.17E3
        G12 = 132.80E3
        G13 = 162.50E3
        v12 = 0.09 
        v13 = 0.04

        ! CRSS (MPa units)
        XTAUC1 = 15.2 ! basal
        XTAUC2 = 67.7 ! prismatic
        XTAUC4 = 2000.0 ! pyramidal

        ! thermal expansion coefficients
        alpha1 = 9.5D-6
        alpha2 = alpha1
        alpha3 = 0.5895*alpha1

        ! prefactor for SSD evolution equation
       	gammast = 0.0

        case default
        WRITE(*,*)"Not sure what material"
        END SELECT

      ! elastic constants based on crystal symmetry
      E2 = E1
      G23 = G13
      v23 = v13

      ! assign Burgers vector scalars
      burgerv(1:6) = burger1
      burgerv(7:12) = burger2

      ! assign CRSS
      tauc(1:3) = XTAUC1
      tauc(4:6) = XTAUC2
      tauc(7:12) = XTAUC4

      case(1) !bcc

        SELECT CASE(matname)

        case('tungsten')

        ! CRSS (MPa units)
        XTAUC1 = 360.0

        ! Burgers vectors (um)
        burger1 = 2.74E-4

        ! elastic constants (MPa units)
        E1 = 421E3
        G12 = 164.4E3
        v12 = 0.28

        ! thermal expansion coefficients
        alpha1 = 9.5e-6
        alpha2 = alpha1
        alpha3 = 0.5895*alpha1

        ! prefactor for SSD evolution equation
       	gammast = 0.0

        case default
        WRITE(*,*)"Not sure what material"
        END SELECT

      ! elastic constants based on crystal symmetry
      E2 = E1
      E3 = E1
      v13 = v12
      v23 = v12
      G13 = G12
      G23 = G12
     
      ! assign Burgers vector scalars
      burgerv = burger1

      ! assign CRSS: same for all slip systems
      tauc = XTAUC1
     
      case(2) !fcc

        SELECT CASE(matname)

        case('copper')

        ! CRSS (MPa units)
        XTAUC1 = 20.0

        ! Burgers vectors (um)
        burger1 = 2.55E-4

        ! elastic constants (MPa units)
        E1 = 66.69E3
        v12 = 0.4189
        G12 = 75.4E3

        ! thermal expansion coefficients
        alpha1 = 13.0e-6
        alpha2 = alpha1
        alpha3 = alpha1

        ! prefactor for SSD evolution equation
        gammast = 0.0

        case default
        WRITE(*,*)"Not sure what material"
        END SELECT

      screwplanes = 2
      
      ! elastic constants based on crystal symmetry
      E2 = E1
      E3 = E1
      v13 = v12
      v23 = v12
      G13 = G12
      G23 = G12
      
      ! assign Burgers vector scalars
      burgerv = burger

      ! assign CRSS: same for all slip systems
      tauc = XTAUC1
      
      case(3) ! carbide

        SELECT CASE(matname)

        case('carbide')

        ! CRSS (MPa units)
        XTAUC1 = 2300.0

        ! Burgers vectors (um)
        burger = 3.5072e-4

        ! elastic constants (MPa units)
        E1 = 207.0E+4
        v12 = 0.28

        ! thermal expansion coefficients
        alpha1 = 4.5e-6
        alpha2 = alpha1
        alpha3 = alpha1

        case default
        WRITE(*,*)"Not sure what material"
        END SELECT
	  
      ! elastic constants based on crystal symmetry      
      E3 = E1
      v13 = v12
      G12 = E1/(2.0*(1.0+v12)) 
      G13 = G12
      G23 = G12
     
      ! assign Burgers vector scalars
      burgerv = burger

      ! assign CRSS: same for all slip systems
      tauc = XTAUC1
      
      case(4) ! olivine

        SELECT CASE(matname)

        case('olivine')

        ! CRSS (MPa units)
        XTAUC1 = 2000.0

        ! Burgers vectors (um)
        burger1 = 2.74E-4

        ! elastic constants (MPa units)
        E1 =  421E3
        v12 = 0.28
        G12 = 164.4E3

        ! thermal expansion coefficients
        alpha1 = 4.5e-6
        alpha2 = alpha1
        alpha3 = alpha1

        case default
        WRITE(*,*)"Not sure what material"
        END SELECT

      ! elastic constants based on crystal symmetry       
      E2 = E1
      E3 = E1
      v13 = v12
      v23 = v12
      G13 = G12
      G23 = G12

      ! assign Burgers vector scalars
      burgerv = burger1

      ! assign CRSS
      tauc(1:7) = XTAUC1*(/2.0,1.0,1.0,1.5,4.0,4.0,1.0/)

      case(5) !alpha-Uranium

        SELECT CASE(matname)

        case('alphauranium')

        ! Burgers vectors (micrometre)
        burgerv(1:2) = 2.85e-4
        burgerv(3:4) = 6.51e-4
        burgerv(5:8) = 11.85e-4

        ! Constant factor tau_0^alpha for CRSS (MPa)
        ! Calhoun 2013 values
        ! with temperature dependence as in Zecevic 2016
        ! added after hardening is included
        tauc(1) = 24.5
        tauc(2) = 85.5
        tauc(3) = 166.5
        tauc(4) = 166.5
        tauc(5) = 235.0
        tauc(6) = 235.0
        tauc(7) = 235.0
        tauc(8) = 235.0

        ! Daniel 1971, Figure 9
        ! softest value (MPa)
        if (twinon == 1) then
          tauctwin(nTwinStart) = 6.25
          tauctwin(nTwinEnd) = 6.25
        end if

        ! elastic moduli (MPa) and Poissons ratios
        ! see PRS Literature Review by Philip Earp
        ! Short Crack Propagation in Uranium, an Anisotropic Polycrystalline Metal
        ! temperature dependence according to Daniel 1971
        E1 = 203665.987780 * (1.0 - 0.000935*(Temperature-293.0))
        E2 = 148588.410104 * (1.0 - 0.000935*(Temperature-293.0))
        E3 = 208768.267223 * (1.0 - 0.000935*(Temperature-293.0))
        v12 = 0.242363
        v13 = -0.016293
        v23 = 0.387816
        G12 = 74349.442379 * (1.0 - 0.000935*(Temperature-293.0))
        G13 = 73421.439060 * (1.0 - 0.000935*(Temperature-293.0))
        G23 = 124378.109453 * (1.0 - 0.000935*(Temperature-293.0))

        ! define thermal expansion coefficients as a function of temperature
        ! Lloyd, Barrett, 1966
        ! Thermal expansion of alpha Uranium
        ! Journal of Nuclear Materials 18 (1966) 55-59
        alpha1 = 24.22e-6 - 9.83e-9 * Temperature + 46.02e-12 * Temperature * Temperature
        alpha2 = 3.07e-6 + 3.47e-9 * Temperature - 38.45e-12 * Temperature * Temperature
        alpha3 = 8.72e-6 + 37.04e-9 * Temperature + 9.08e-12 * Temperature * Temperature

        case default
        WRITE(*,*)"Not sure what material"
        END SELECT

      case default
      WRITE(*,*)"Not sure what crystal type"
      END SELECT

C     *** SET UP ELASTIC STIFFNESS MATRIX IN LATTICE SYSTEM ***  
C     notation as in: https://en.wikipedia.org/wiki/Hooke%27s_law 
C     However, the Voigt notation in Abaqus for strain and stress vectors
C     has the following order:
C     stress = (sigma11,sigma22,sigma33,tau12 tau13 tau23)
C     strain = (epsilon11,epsilon22,epsilon33,gamma12,gamma13,gamma23)

      compliance(1,1:3) = (/1./E1,-v12/E1,-v13/E1/)
      compliance(2,2:3) =         (/1./E2,-v23/E2/)
      compliance(3,3:3) =                 (/1./E3/)
      compliance(4,4:4) =                       (/1./G12/)
      compliance(5,5:5) =                       (/1./G13/)
      compliance(6,6:6) =                       (/1./G23/)

      ! symmetrize compliance matrix
      DO i=2,6
         DO j=1,i-1
            compliance(i,j)=compliance(j,i)
         END DO
      END DO

      ! define thermal eigenstrain in the lattice system
      thermat(1,1) = alpha1 
      thermat(2,2) = alpha2 
      thermat(3,3) = alpha3

      RETURN

      END

