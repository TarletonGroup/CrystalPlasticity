********************************************************
** KCRSS calculates the CRSS of slip and twin systems **
********************************************************

      SUBROUTINE kCRSS(iphase,tauc,nSys,G12,burgerv,gndtot,irradiate,
     + tauSolute,gndcut,rhofor,rhosub,Temperature,homogtwin,
     + nTwinStart,nTwinEnd,twinvolfrac,tauctwin,nTwin,TwinIntegral,
     + twinvolfractotal,twinon)

      INCLUDE 'ABA_PARAM.INC'

      ! crystal type
      INTEGER,intent(in) :: iphase

      ! number of slip systems
      INTEGER,intent(in) :: nSys

      ! total number of twin systems
      INTEGER,intent(in) :: nTwin

      ! shear modulus for Taylor dislocation law
      REAL*8,intent(in) :: G12

      ! Burgers vectors
      REAL*8,intent(in) :: burgerv(nSys)

      ! scalar total GND density
      REAL*8,intent(in) :: gndtot

      ! activate irradiation effect
      INTEGER,intent(in) :: irradiate

      ! increase in tauc due to solute force
      REAL*8,intent(in) :: tauSolute

      ! GND density (immobile)
      REAL*8,intent(in) :: gndcut(nSys)

      ! forest dislocation density
      REAL*8,intent(in) :: rhofor(nSys)

      ! substructure dislocation density
      REAL*8,intent(in) :: rhosub 

      ! Current temperature
      REAL*8,intent(in) :: Temperature

      ! homogenize twin model
      INTEGER,intent(in) :: homogtwin

      ! the active twins are the ones in the
      ! interval [nTwinStart,nTwinEnd] in the
      ! twin system file
      INTEGER,intent(in) :: nTwinStart,nTwinEnd
	  
	  ! twin systems activation flag
      INTEGER,intent(in) :: twinon

      ! twin volume fraction
      REAL*8,intent(in) :: twinvolfrac(nTwin)

      ! average of the twin volume fraction
      ! over the neighbourhood
      ! two twin systems
      REAL*8,intent(in) :: TwinIntegral(nTwin)

      ! total twin volume fraction 
      REAL*8,intent(in) :: twinvolfractotal

      ! critical resolved shear stress of slip systems
      REAL*8,intent(inout) :: tauc(nSys)

      ! critical resolved shear stress of twin systems
      REAL*8,intent(inout) :: tauctwin(nTwin)

      INTEGER :: i

      ! check crystal type
      if (iphase == 1) then

          ! Taylor dislocation law
          tauc = tauc + 0.0065*G12*(burgerv(1))*sqrt(gndtot)
         
          if (irradiate == 1) then
              tauc = tauc + tauSolute
          end if

      else if (iphase == 2) then

          ! Taylor dislocation law
          tauc = tauc + 0.32*G12*(burgerv(1))*sqrt(gndcut)

      else if (iphase == 5) then 

          ! alpha-Uranium model with forest and substructure dislocations
          ! R.J. McCabe, L. Capolungo, P.E. Marshall, C.M. Cady, C.N. Tomé
          ! Deformation of wrought uranium: Experiments and modeling
          ! Acta Materialia 58 (2010) 5447–5459
          tauc(1) = tauc(1) + 19.066 * sqrt(rhofor(1)) + 1.8218 * sqrt(rhosub) * log(1.0 / (burgerv(1) * sqrt(rhosub)))
          tauc(2) = tauc(2) + 18.832 * sqrt(rhofor(2)) + 1.7995 * sqrt(rhosub) * log(1.0 / (burgerv(2) * sqrt(rhosub)))
          tauc(3) = tauc(3) + 54.052 * sqrt(rhofor(3)) + 5.1650 * sqrt(rhosub) * log(1.0 / (burgerv(3) * sqrt(rhosub)))
          tauc(4) = tauc(4) + 54.052 * sqrt(rhofor(4)) + 5.1650 * sqrt(rhosub) * log(1.0 / (burgerv(4) * sqrt(rhosub)))
          tauc(5) = tauc(5) + 123.357 * sqrt(rhofor(5)) + 11.7875 * sqrt(rhosub) * log(1.0 / (burgerv(5) * sqrt(rhosub)))
          tauc(6) = tauc(6) + 123.357 * sqrt(rhofor(6)) + 11.7875 * sqrt(rhosub) * log(1.0 / (burgerv(6) * sqrt(rhosub)))
          tauc(7) = tauc(7) + 123.357 * sqrt(rhofor(7)) + 11.7875 * sqrt(rhosub) * log(1.0 / (burgerv(7) * sqrt(rhosub)))
          tauc(8) = tauc(8) + 123.357 * sqrt(rhofor(8)) + 11.7875 * sqrt(rhosub) * log(1.0 / (burgerv(8) * sqrt(rhosub)))

          ! Zecevic 2016 temperature dependence
          tauc(1) = tauc(1) * exp(-(Temperature-293.0)/140.0)
          tauc(2) = tauc(2) * exp(-(Temperature-293.0)/140.0)
          tauc(3) = tauc(3) * exp(-(Temperature-293.0)/140.0)
          tauc(4) = tauc(4) * exp(-(Temperature-293.0)/140.0)
          tauc(5) = tauc(5) * exp(-(Temperature-293.0)/140.0)
          tauc(6) = tauc(6) * exp(-(Temperature-293.0)/140.0)
          tauc(7) = tauc(7) * exp(-(Temperature-293.0)/140.0)
          tauc(8) = tauc(8) * exp(-(Temperature-293.0)/140.0)

          ! Daniel, Lesage, 1971 minimum value as minima
          tauc(1) = max(tauc(1),4.0)
          tauc(2) = max(tauc(2),4.0)
          tauc(3) = max(tauc(3),4.0)
          tauc(4) = max(tauc(4),4.0)
          tauc(5) = max(tauc(5),4.0)
          tauc(6) = max(tauc(6),4.0)
          tauc(7) = max(tauc(7),4.0)
          tauc(8) = max(tauc(8),4.0)

          if (twinon == 1) then ! twin active

            if (homogtwin == 1) then ! homogenized twin model

	      ! add cross hardening of one twin system on the other
              DO i=nTwinStart,nTwinEnd
                tauctwin(i) = tauctwin(i) + 96.79*twinvolfractotal
              END DO

            else ! discrete twin model

              ! add hardening in the nucleation stage
              ! 50% is the critical twin volume fraction at which the
              ! softest value is reached
              DO i=nTwinStart,nTwinEnd
                if (twinvolfrac(i) < 0.5) then
                  tauctwin(i) = tauctwin(i) + 37.5*(0.5-twinvolfrac(i))
                  tauctwin(i) = tauctwin(i) + 2000.0*TwinIntegral(i)
                else ! twinvolfrac(i) > 0.5
                  tauctwin(i) = tauctwin(i) + 37.5*(twinvolfrac(i)-0.5)
                  tauctwin(i) = tauctwin(i) + 2000.0*TwinIntegral(i)
                end if
              END DO

              ! local interaction between twin systems
              ! only when threshold is passed
              if (twinvolfrac(nTwinEnd) > 0.5) then
                tauctwin(nTwinStart) = tauctwin(nTwinStart) + 200.0*twinvolfrac(nTwinEnd)
              end if
              if (twinvolfrac(nTwinStart) > 0.5) then
                tauctwin(nTwinEnd) = tauctwin(nTwinEnd) + 200.0*twinvolfrac(nTwinStart)
              end if

            end if ! choice of twin model (homogeneous/discrete)

          end if ! twin active

      end if ! check crystal type

      RETURN

      END
