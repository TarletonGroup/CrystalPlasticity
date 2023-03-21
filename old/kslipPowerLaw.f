** Nicolo Grilli
** 28th January 2019
** alpha-Uranium
** AWE project
**
** including twinning
** NO rotation of slip systems due to twinning
** but slip is allowed inside twins
**
** for full twin developement a characteristic
** time is introduced if twin volume fraction > 0.5
** This leads the twin volume fraction to reach 1.0
** after about the characteristic time
**

      subroutine kslipPowerLaw(xNorm,xDir,xTwinNorm,xTwinDir,
     + tau,tautwin,signtau,signtautwin,tauc,tauctwin,
     + burgerv,dtime,nSys,nTwin,iphase,irradiate,
     + gndcut,gndtot,rhossd,twinvolfrac,twinvolfractotal,
     + Lp,tmat,gammaDot,gammaTwinDot,twinon,
     + nTwinStart,nTwinEnd)

      implicit none
      integer,intent(in) :: nSys,nTwin,iphase,irradiate,twinon
      integer,intent(in) :: nTwinStart,nTwinEnd
      real*8, intent(in) :: dtime, signtau(nSys),signtautwin(nTwin), gndcut(nSys), gndtot,rhossd
      real*8, intent(in) :: xNorm(nSys,3),xDir(nSys,3),tau(nSys),tauc(nSys),burgerv(nSys)
      real*8, intent(in) :: xTwinNorm(nTwin,3),xTwinDir(nTwin,3),tautwin(nTwin),tauctwin(nTwin)
      real*8, intent(in) :: twinvolfrac(nTwin), twinvolfractotal ! twin volume fraction

      real*8, intent(inout) :: gammaTwinDot(nTwin)
      real*8, intent(out) :: Lp(3,3),tmat(6,6), gammaDot(nSys)
      
      integer :: i, j
      real*8  :: alpha,beta,result1, rhom,dF,f,T,k,gamma0,b,psi,V,
     + SNij(3,3),sni(6),nsi(6),SNNS(6,6),NSij(3,3),result4(6,6),
     + tempNorm(3), tempDir(3)

      ! alpha-Uranium shear strain rate law
      real*8 :: gamma0dot, gamma0dottwin ! slip rate prefactor
      real*8 :: nsliprate ! slip rate exponent
      real*8 :: xtau ! ratio: resolved shear stress / CRSS
      
C
C  *** CALCULATE THE DERIVATIVE OF PLASTIC STRAIN INCREMENT WITH
C   RESPECT TO THE STRESS DEFINED AS tmat***
C     
      

      gamma0dot = 1.0e-3
      gamma0dottwin = 1.0e-3
      nsliprate = 20.0

      
      !ne = microns, stress = MPa, F = microN, therefore E = pJ

C
      tmat=0.; Lp = 0.;result4=0.

      Do I=1,nSys ! shear due to slip without twinning

         xtau = tau(I) / tauc(I)

         ! 0.5 for alpha-Uranium model due to the exponent 20 in the shear rate law	
         ! slip is allowed inside twins
         if (xtau >= 0.5) THEN

              gammaDot(I)=gamma0dot*(xtau**nsliprate)*signtau(I) ! signed
         
              tempNorm = xNorm(I,:); 
              tempDir = xDir(I,:)
              SNij = spread(tempDir,2,3)*spread(tempNorm,1,3) ! b_i tensor product n_j
              NSij = spread(tempNorm,2,3)*spread(tempDir,1,3) ! n_i b_j
              CALL KGMATVEC6(SNij,sni)         
              CALL KGMATVEC6(NSij,nsi)
              ! nsi and sni are the same thing
              SNNS = spread(sni,2,6)*spread(nsi,1,6)
              result1 = xtau**(nsliprate-1.0)
          
              ! result4 = dt dL_{ij}/dsigma_{kl} with indices in Voigt notation
              ! such that 4,5,6 components contributes as (b_i n_j + b_j n_i)
              ! slip inside twin is allowed
              result4=result4+(gamma0dot*nsliprate/tauc(I))*dtime*result1*SNNS

              ! update plastic velocity gradient
              ! slip is allowed inside twins
              Lp = Lp + gammaDot(I)*SNij

          else
              gammaDot(I)=0.0
          end if
C
      END DO

      if (twinon == 1) then ! twin active

        Do I=nTwinStart,nTwinEnd ! shear due to twinning

          xtau = tautwin(I) / tauctwin(I)

          ! 0.5 for alpha-Uranium model due to the exponent 20 in the shear rate law
          ! if twins have already occupied the whole volume
          ! no additional plastic strain from twinning is possible
          if ((xtau >= 0.5 .or. twinvolfrac(I) > 0.5) .and. twinvolfrac(I) < 1.0) THEN

              gammaTwinDot(I)=gamma0dottwin*(xtau**nsliprate)*signtautwin(I) ! signtautwin always positive

              if (gammaTwinDot(I) > 1000.0*gamma0dottwin) then
                gammaTwinDot(I) = 1000.0*gamma0dottwin
              end if

              ! add contribution if activation threshold 0.5 has been passed
              ! here characteristic time is 1s
              if (twinvolfrac(I) > 0.5) then
                gammaTwinDot(I)=gammaTwinDot(I)+0.299*(1.0-twinvolfrac(I))
              end if

              tempNorm = xTwinNorm(I,:);
              tempDir = xTwinDir(I,:)
              SNij = spread(tempDir,2,3)*spread(tempNorm,1,3) ! b_i tensor product n_j
              NSij = spread(tempNorm,2,3)*spread(tempDir,1,3) ! n_i b_j
              CALL KGMATVEC6(SNij,sni)         
              CALL KGMATVEC6(NSij,nsi)
              ! nsi and sni are the same thing
              SNNS = spread(sni,2,6)*spread(nsi,1,6)
              if (gammaTwinDot(I) > 999.0*gamma0dottwin) then
                  result1 = 707.9
              else
                  result1 = xtau**(nsliprate-1.0)
              end if

              ! result4 = dt dL_{ij}/dsigma_{kl} with indices in Voigt notation
              ! such that 4,5,6 components contributes as (b_i n_j + b_j n_i)
              result4 = result4 + (gamma0dottwin*nsliprate/tauctwin(I))*dtime*result1*SNNS

              ! update plastic velocity gradient
              Lp = Lp + gammaTwinDot(I)*SNij

          else
              gammaTwinDot(I)=0.0
          end if

        END DO

      end if ! twin active

      tmat = 0.5*(result4+transpose(result4))
             
      return
      end 
