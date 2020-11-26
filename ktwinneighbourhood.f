C Non-local twin model
C find neigbourhood of an integration point
C and calculate the integral of the twin volume fraction
C on that neighbourhood

C see: Grilli, Cocks, Tarleton  
C A phase field model for the growth and characteristic thickness of deformation-induced twins
C Journal of the Mechanics and Physics of Solids, Volume 143, October 2020, 104061

C calculation of TwinUpDownEl, TwinUpDownIP
C at the second increment when kgausscoords
C have been assigned for all elements and IPs
      SUBROUTINE kfindneighbourhood(noel,npt,COORDS,STATEV,NSTATV,
     + Twin1UpDownEl,Twin1UpDownIP,Twin2UpDownEl,Twin2UpDownIP,
     + nTwinStart,nTwinEnd,nTwin)

      integer,intent(in) :: NSTATV

      real*8,intent(in) :: COORDS(3)
      real*8,intent(in) :: STATEV(NSTATV)

      integer,intent(in) :: noel,npt

      integer,intent(in) :: nTwinStart
      integer,intent(in) :: nTwinEnd

      integer,intent(in) :: nTwin

      include 'mycommon.f'

      integer,intent(inout) :: Twin1UpDownEl(ArrayNUpDown), Twin1UpDownIP(ArrayNUpDown)
      integer,intent(inout) :: Twin2UpDownEl(ArrayNUpDown), Twin2UpDownIP(ArrayNUpDown)

      ! temporary vector from current IP
      ! to neighbouring IP and its norm
      real*8 :: NeighbourConnection(3)
      real*8 :: NeighbourConnectionAbs

      ! temporary gauss point coordinates
      real*8 :: gausscoords(3)

      ! rotation matrix needed to
      ! rotate twin directions and normals
      real*8 :: gmatinv(3,3)

      ! temporary twin normal
      real*8 :: ttwinnor(3)

      ! twin normals
      real*8 :: dtwinnor(nTwin,3)

      ! temporary projection of the distance
      ! from neighbouring IP on the twin plane normal
      real*8 :: NeighbourProjection

      ! temporary distance from neighbouring IP
      ! to the axis of the cylinder passing through
      ! the current IP and parallel to the twin plane normal
      real*8 :: NeighbourDistCyl

      ! temporary indices of the neighbouring IPs
      ! two discrete twin systems per grain
      ! but indices can be arbitrary based on nTwinStart and nTwinEnd
      integer :: NeighbourUpDownIndex(nTwin)

      ! flag telling that NUpDown has been exceeded
      ! not enough space in the NUpDown variable
      ! for neighbouring points
      integer :: FlagNUpDown

      ! temporary indices of SMAIntArray
      ! converting noel, npt and NeighbourUpDownIndex
      ! into the SMA array index
      integer :: ArrayIndexUpDown

      integer :: elem, ip, i, j, twinindex

      ! twin normals are necessary
      ! to determine the position of an IP with
      ! respect to a formed twin
      include 'xTwinNormAlphaUranium.f'

        NeighbourUpDownIndex(1:nTwin) = 0
        FlagNUpDown = 0
! loop over all elements and IPs

        do elem=1,nElements
          do ip=1,nintpts

            NeighbourConnection(1:3) = 0.0
            NeighbourConnectionAbs = 0.0

            if (elem /= noel .or. ip /= npt) then ! check that IP is not the current one
                
              if (kGrainIndex(noel,npt) == kGrainIndex(elem,ip)) then ! check that IP belongs to the same grain

                ! assign IP coordinate to a temporary variable
                gausscoords(1) = kgausscoords(elem,ip,1)
                gausscoords(2) = kgausscoords(elem,ip,2)
                gausscoords(3) = kgausscoords(elem,ip,3)

                ! check that x coordinate is close to
                ! the current element
                if (abs(COORDS(1) - gausscoords(1)) < 10.001) then

                  ! check that y coordinate is close to
                  ! the current element
                  if (abs(COORDS(2) - gausscoords(2)) < 10.001) then

                    ! check that z coordinate is close to
                    ! the current element
                    if (abs(COORDS(3) - gausscoords(3)) < 1.001) then

                      ! calculate vector from current IP to
                      ! neighbouring IP
                      NeighbourConnection(1) = gausscoords(1) - COORDS(1)
                      NeighbourConnection(2) = gausscoords(2) - COORDS(2)
                      NeighbourConnection(3) = gausscoords(3) - COORDS(3)
                      NeighbourConnectionAbs = dot_product(NeighbourConnection(1:3),NeighbourConnection(1:3))
                      NeighbourConnectionAbs = sqrt(NeighbourConnectionAbs)

                      ! assign rotation matrix to gmatinv
                      DO i=1,3
                        DO j=1,3
                          gmatinv(i,j) = STATEV(j+(i-1)*3)
                        END DO
                      END DO

                      DO twinindex=nTwinStart,nTwinEnd ! cycle over twin systems

                        ! dtwinnor is assumed normalized
                        ! rotate it to the sample reference frame
                        DO i=1,3
                          ttwinnor(i) = dtwinnor(twinindex,i) ! already normalized
                        END DO
                        ttwinnor = matmul(gmatinv,ttwinnor)

                        ! calculate projection of NeighbourConnection
                        ! on the twin plane normal
                        NeighbourProjection = abs(dot_product(NeighbourConnection(1:3),ttwinnor(1:3)))

                        NeighbourDistCyl = NeighbourConnectionAbs*NeighbourConnectionAbs
                        NeighbourDistCyl = NeighbourDistCyl - NeighbourProjection*NeighbourProjection
                        NeighbourDistCyl = sqrt(NeighbourDistCyl)

                        ! check if the neighbouring point
                        ! is up/down (inside the cylinder or not)
                        if (NeighbourDistCyl < 1.001 .and. NeighbourProjection < 5.001) then ! up/down
                          NeighbourUpDownIndex(twinindex) = NeighbourUpDownIndex(twinindex) + 1
                          if (NeighbourUpDownIndex(twinindex) <= NUpDown) then
      ArrayIndexUpDown = (noel-1)*nintpts*NUpDown+(npt-1)*NUpDown+NeighbourUpDownIndex(twinindex)
                            if (twinindex == nTwinStart) then
                              Twin1UpDownEl(ArrayIndexUpDown) = elem
                              Twin1UpDownIP(ArrayIndexUpDown) = ip
                            end if
                            if (twinindex == nTwinEnd) then
                              Twin2UpDownEl(ArrayIndexUpDown) = elem
                              Twin2UpDownIP(ArrayIndexUpDown) = ip
                            end if
                          else
                            FlagNUpDown = 1 ! NUpDown needs to be increased
                          end if
                        end if ! check up/down

                      END DO ! cycle over twin systems

                    end if ! check z coordinate

                  end if ! check y coordinate

                end if ! check x coordinate
              
              end if ! check that IP belongs to the same grain
              
            end if ! check that IP is not the current one

          end do ! loop over all IPs
        end do ! loop over all elements

        ! assign number of Neighbours to common block
        kNoNeighbours(noel,npt,1) = NeighbourUpDownIndex(nTwinStart)
        kNoNeighbours(noel,npt,2) = NeighbourUpDownIndex(nTwinEnd)
        
        if (FlagNUpDown == 1) then
          if (maxval(NeighbourUpDownIndex) > NUpDownExceeded) then
            NUpDownExceeded = maxval(NeighbourUpDownIndex)
            write(*,*) "WARNING: number of up/down neighbouring IPs exceeded"
            write(*,*) "Max value reached by NeighbourUpDownIndex(nTwinStart:nTwinEnd) is:"
            write(*,*) "maxval(NeighbourUpDownIndex)"
            write(*,*) maxval(NeighbourUpDownIndex)
          end if
        end if

      RETURN

      END

C calculate integral of the twin volume fraction
C in the up/down region

      SUBROUTINE ktwinvolfracintegral(noel,npt,TwinIntegral,
     + Twin1UpDownEl,Twin1UpDownIP,Twin2UpDownEl,Twin2UpDownIP,
     + nTwinStart,nTwinEnd,nTwin)

      integer,intent(in) :: nTwinStart
      integer,intent(in) :: nTwinEnd

      integer,intent(in) :: nTwin

      real*8,intent(inout) :: TwinIntegral(nTwin)

      integer,intent(in) :: noel,npt

      include 'mycommon.f'

      integer,intent(in) :: Twin1UpDownEl(ArrayNUpDown), Twin1UpDownIP(ArrayNUpDown)
      integer,intent(in) :: Twin2UpDownEl(ArrayNUpDown), Twin2UpDownIP(ArrayNUpDown)

      integer :: I, twinindex

      ! temporary indices of SMAIntArray
      ! converting noel, npt and NeighbourUpDownIndex
      ! into the SMA array index
      integer :: ArrayIndexUpDown

      ! temporary element index
      ! to find neighbouring IPs
      integer :: noelTwin
      ! temporary IP index
      ! to find neighbouring IPs
      integer :: nptTwin

      DO twinindex=nTwinStart,nTwinEnd ! cycle over twin systems

        DO I=1,NUpDown ! cycle over up/down neighbours
          ArrayIndexUpDown = (noel-1)*nintpts*NUpDown+(npt-1)*NUpDown+I
          if (twinindex == nTwinStart) then
            noelTwin = Twin1UpDownEl(ArrayIndexUpDown)
            nptTwin = Twin1UpDownIP(ArrayIndexUpDown)
          end if
          if (twinindex == nTwinEnd) then
            noelTwin = Twin2UpDownEl(ArrayIndexUpDown)
            nptTwin = Twin2UpDownIP(ArrayIndexUpDown)
          end if

          if (noelTwin /= 0 .and. nptTwin /= 0) then ! this IP was not assigned: e.g. element close to the boundary
            if (twinindex == nTwinStart) then
              TwinIntegral(twinindex) = TwinIntegral(twinindex) + kTwinVolFrac(noelTwin,nptTwin,1)
            end if
            if (twinindex == nTwinEnd) then
              TwinIntegral(twinindex) = TwinIntegral(twinindex) + kTwinVolFrac(noelTwin,nptTwin,2)
            end if
          end if
        END DO ! cycle over up/down neighbours

      END DO ! end cycle over twin systems

      ! average over neighbours
      if (kNoNeighbours(noel,npt,1) == 0) then ! kNoNeighbours is not assigned at the first increment
        TwinIntegral(nTwinStart) = 0.0
      else
        TwinIntegral(nTwinStart) = TwinIntegral(nTwinStart) / kNoNeighbours(noel,npt,1)
      end if
      if (kNoNeighbours(noel,npt,2) == 0) then ! kNoNeighbours is not assigned at the first increment
        TwinIntegral(nTwinEnd) = 0.0
      else
        TwinIntegral(nTwinEnd) = TwinIntegral(nTwinEnd) / kNoNeighbours(noel,npt,2)
      end if

      RETURN

      END






