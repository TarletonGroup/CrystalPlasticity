
      subroutine kslip6ET(xNorm,xDir,tau,signtau,tauc,burgerv,dtime,
     + nSys,iphase,irradiate, gndcut,gndtot,rhossd,Lp,tmat,gammaDot)
!Slip Rule with constant coefficients for simplicity

         
      implicit none
      integer,intent(in):: nSys,iphase, irradiate 
      real*8, intent(in) :: dtime, signtau(nSys), gndcut(nSys), gndtot,rhossd
      real*8, intent(in) :: xNorm(nSys,3),xDir(nSys,3),tau(nSys),tauc(nSys),burgerv(nSys)

      real*8, intent(out) :: Lp(3,3),tmat(6,6), gammaDot(nSys)
      
      integer :: i
      real*8  :: alpha,beta,result1, rhom,dF,f,T,k,gamma0,b,psi,V,
     + SNij(3,3),sni(6),nsi(6),SNNS(6,6),NSij(3,3),result4(6,6),
     + tempNorm(3), tempDir(3)
      
C
C  *** CALCULATE THE DERIVATIVE OF PLASTIC STRAIN INCREMENT WITH 
C   RESPECT TO THE STRESS DEFINED AS tmat***
C     
      if (iphase == 1) then
          
          rhom =  0.035 !mobile dislocation density 
          dF = 3.4559E-8 !! Energy barrier microN.microns = pJ          
          f = 1.0E11 ! attempt frequency /s
          T = 293.0
          k = 1.381E-11 ! pJ / K
          gamma0 = 8.33E-6 ! some reference strain which appears in eg http://dx.doi.org/10.1016/j.ijsolstr.2015.02.023
          b = burgerv(1)                    
             
          if (irradiate == 1) then              
              psi = 3.457e-2   
          else
              psi = 0.727e-2  
          end if          
          
          alpha = rhom*b*b*f*exp(-dF/(k*T)) ! prefactor /s
          
          V = b*b/sqrt(psi*rhossd) ! activation volume /micron^3
          
          beta = gamma0*V/(k*T) ! rate sensitivity 1/MPa
          
      elseif (iphase == 2) then

          alpha = 0.02
          beta = 0.1
          
      endif
          
      !ne = microns, stress = MPa, F = microN, therefore E = pJ
    
C
      tmat=0.; Lp = 0.;result4=0.
	     
      !rhogndold=sum(gndold)
      Do I=1,nSys

         if (tau(I) >= tauc(I)) THEN				
                                                 
             gammaDot(I)=alpha*sinh(beta*signtau(I)*(tau(I) - tauc(I)) )
         
              tempNorm = xNorm(I,:); 
              tempDir = xDir(I,:)
              SNij = spread(tempDir,2,3)*spread(tempNorm,1,3)
              NSij = spread(tempNorm,2,3)*spread(tempDir,1,3)
              CALL KGMATVEC6(SNij,sni)         
              CALL KGMATVEC6(NSij,nsi) 
              SNNS = spread(sni,2,6)*spread(nsi,1,6)
              result1 = cosh(beta*signtau(I)*(tau(I) - tauc(I)))
          
              result4 = result4 + alpha*beta*dtime*result1*SNNS          
              Lp = Lp + gammaDot(I)*SNij
          else
              gammaDot(I)=0.0
          end if
C
      END DO
      
      tmat = 0.5*(result4+transpose(result4))
             
      return
      end 
