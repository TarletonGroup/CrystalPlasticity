! Nicolo Grilli
! University of Oxford
! AWE project 2020
! 7th April 2020

! calculate integral of the twin volume fraction
! in a cylinder parallel to the twin plane normal
! two twin systems

! use allocatable arrays to avoid passing
! the 2 GB limit for static arrays in the stack
! as common blocks do
! use the command SMALocalIntArrayCreateInit(ID,SIZE,0)

!Use this to avoid editing multiple common blocks independently
		integer, parameter :: nElements = 8035
		integer, parameter :: nintpts = 8
		integer, parameter :: nsdv  = 125
		! number of cohesive elements
		integer, parameter :: nElemUel = 1
		! maximum number of IPs in the cylinder
		integer, parameter :: NUpDown = 3564
		! NUpDown * nElements * nintpts
        integer, parameter :: ArrayNUpDown = 229093920
        real*8 :: kFp, kgausscoords, kcurlFp, kTwinVolFrac
        real*8 :: kSigma0
        real*8 :: kDeltaEff
        integer :: kNoNeighbours

        ! grain index for each element and IP
        integer :: kGrainIndex, NUpDownExceeded

      COMMON/UMPS/kFp(nElements, nintpts, 9),
     +    kgausscoords(nElements,nintpts,3),
     +    kcurlFp(nElements, nintpts, 9),
     +    kTwinVolFrac(nElements,nintpts,2), ! twin volume fraction
     +    kNoNeighbours(nElements,nintpts,2), ! number of up/down neighbours
     +    kGrainIndex(nElements,nintpts),
     +    kSigma0(nElements,nintpts), ! output the max stress to failure
     +    kDeltaEff(nElements,nintpts), ! output the max effective opening
     +    NUpDownExceeded ! keep track of the max value of the neighbouring IP reached



