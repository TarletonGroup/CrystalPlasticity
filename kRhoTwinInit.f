************************************
** Initialize dislocation density **
**     and twin phase field       **
************************************

      SUBROUTINE kRhoTwinInit(nTwin,M,NSTATV,STATEV,twinon,nTwinStart,nTwinEnd,
     + iphase,COORDS,nSys)

      INCLUDE 'ABA_PARAM.INC'

      ! dimension of the space
      INTEGER,intent(in) :: M

      ! total number of twin systems
      INTEGER,intent(in) :: nTwin

      ! twin activation flag
      INTEGER,intent(in) :: twinon

      ! the active twins are the ones in the
      ! interval [nTwinStart,nTwinEnd] in the
      ! twin system file
      INTEGER,intent(in) :: nTwinStart
      INTEGER,intent(in) :: nTwinEnd

      ! crystal type
      INTEGER,intent(in) :: iphase

      ! coordinates of this IP
      REAL*8,intent(in) :: COORDS(3)

      ! number of slip systems
      INTEGER,intent(in) :: nSys

      REAL*8,intent(inout) :: STATEV(NSTATV)

      ! initial random twin phase field
      real*8 :: TwinInitRand

      ! temporary twin directions and normals
      real*8 :: ttwindir(M)
      real*8 :: ttwinnor(M)

      ! rotation matrix needed to
      ! rotate twin directions and normals
      real*8 :: gmatinv(M,M)

      ! twin directions and normals
      real*8 :: dtwindir(nTwin,M)
      real*8 :: dtwinnor(nTwin,M)

      integer :: i, j

      ! twin directions are necessary
      ! to determine the Schmid tensor
      ! associated with the initial twins
      include 'xTwinDirAlphaUranium.f'

      ! twin normals are necessary
      ! to determine the position of an IP with
      ! respect to a formed twin
      include 'xTwinNormAlphaUranium.f'

      ! get rotation matrix from state variable
      DO i=1,3
        DO j=1,3
          gmatinv(i,j) = STATEV(j+(i-1)*3)
        END DO
      END DO

      if (twinon == 1) then ! twin active

        DO i=nTwinStart,nTwinEnd ! cycle over active twin systems

	  ! initial random twin volume fraction
          CALL RANDOM_NUMBER(TwinInitRand)
          TwinInitRand = 0.001*TwinInitRand

          ! twin of 2nd system at 45 degrees
          if (i == nTwinEnd) then
            if (abs(COORDS(1)-COORDS(2)) < 1.414) then
              TwinInitRand = 0.49
            end if
          end if
         
          ! rotate twin direction and normal
          ! to find the Schmid tensor of the twin in the
          ! sample reference frame, use matmul(M,v) function
          DO j=1,3
            ttwindir(j) = dtwindir(i,j) ! already normalized
            ttwinnor(j) = dtwinnor(i,j) ! already normalized
          END DO
          ttwindir = matmul(gmatinv,ttwindir)
          ttwinnor = matmul(gmatinv,ttwinnor)

	  ! initial plastic deformation
          ! found by tensor product of ttwindir and ttwinnor
          STATEV(81) = STATEV(81) + 0.299*ttwindir(1)*ttwinnor(1)*TwinInitRand
          STATEV(82) = STATEV(82) + 0.299*ttwindir(1)*ttwinnor(2)*TwinInitRand
          STATEV(83) = STATEV(83) + 0.299*ttwindir(1)*ttwinnor(3)*TwinInitRand
          STATEV(84) = STATEV(84) + 0.299*ttwindir(2)*ttwinnor(1)*TwinInitRand
          STATEV(85) = STATEV(85) + 0.299*ttwindir(2)*ttwinnor(2)*TwinInitRand
          STATEV(86) = STATEV(86) + 0.299*ttwindir(2)*ttwinnor(3)*TwinInitRand
          STATEV(87) = STATEV(87) + 0.299*ttwindir(3)*ttwinnor(1)*TwinInitRand
          STATEV(88) = STATEV(88) + 0.299*ttwindir(3)*ttwinnor(2)*TwinInitRand
          STATEV(89) = STATEV(89) + 0.299*ttwindir(3)*ttwinnor(3)*TwinInitRand

          ! random initial twin volume fraction
          STATEV(106+i) = TwinInitRand

        END DO ! end cycle over active twin systems

      end if ! twin active

! see: Grilli, Cocks, Tarleton  
! A phase field model for the growth and characteristic thickness of deformation-induced twins
! Journal of the Mechanics and Physics of Solids, Volume 143, October 2020, 104061
      if (iphase == 5) then

        ! initial substructure dislocation density, alpha-Uranium
        STATEV(65) = 4.0 
        
        ! initial forest dislocation densities
        do i=1,nSys
          STATEV(56+i) = 4.0
        end do

      end if

      RETURN

      END

