!Slip Rule Re-development
      subroutine kslip5ET(xNorm,xDir,tau,signtau,tauc,burgerv,rhossd,gndtot,
     +     gndall,gndcut,gndmob,dtime,ne,ns,iphase,Lp,tmat)
                 

         
      implicit none
      integer,intent(in):: ne,ns,iphase 
      real*8, intent(in) :: rhossd,gndtot,dtime, gndall(ne+ns), gndcut(ne), gndmob(ne), signtau(ne)
      real*8, intent(in) :: xNorm(ne,3),xDir(ne,3),tau(ne),tauc(ne),burgerv(ne)

      real*8, intent(out) :: Lp(3,3),tmat(6,6)
      
      integer :: i
      real*8  :: alpha,beta,result1,
     + xlambdap,xvol,rhom,rhom0,psi,K,T,dF,f,
     + SNij(3,3),sni(6),nsi(6),SNNS(6,6),NSij(3,3),result4(6,6),
     + tempNorm(3), tempDir(3),gammaDot(ne), gamm0, rhognd
      
C
C  *** CALCULATE THE DERIVATIVE OF PLASTIC STRAIN INCREMENT WITH 
C   RESPECT TO THE STRESS DEFINED AS tmat***
C      
      if (iphase == 0) then
          !Slava Be    
          psi = 1.0
          rhom0 = 0.003 !Slava Be
          dF = 3.4559E-20 *1E12 ! microN.microns = pJ
          f = 1e11
          gamm0 = 5e-2 
          
      elseif (iphase == 2) then
          psi = 1E-7 !1.457e-4
          rhom0 = 5.0      
          dF = 0.0 !  !3.4559E-20 *1E12 ! microN.microns = pJ
          f = 50E+11 !1e11
          gamm0 = 6e-4  !8.33E-6 
      endif
      
      K = 1.381E-23 *1E12 ! pJ / K
      T = 293.0 
      !ne = microns, stress = MPa, F = microN, therefore E = pJ
    
C
      tmat=0.; Lp = 0.;result4=0.
	     
      !rhogndold=sum(gndold)
      Do I=1,ne

         if (tau(I) >= tauc(I)) THEN				
             
             rhom = psi*(rhom0 + gndmob(I))
             rhognd = gndtot !gndcut(I)
            
             xlambdap = 1.0/sqrt((rhognd+rhossd)) !overall     
             xvol = xlambdap*burgerv(I)*burgerv(I)
         
             beta = gamm0*xvol/(K*T) 
          
             alpha = rhom*burgerv(I)*burgerv(I)*f*exp(-dF/(K*T))
          
             gammaDot(I)=alpha*sinh(beta*signtau(I)*(tau(I) - tauc(I)) )
         
              tempNorm = xNorm(I,:); tempDir = xDir(I,:)
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
      end subroutine kslip5ET
