C   Nicolo Grilli
C   University of Bristol
C   Christos Skamniotis
C   University of Oxford
C   11 Novembre 2021 
C
C   Double exponent slip rule in:
C   Zhengxuan Fan & Serge Kruch (2020) A comparison of different crystal
C   plasticity finite-element models on the simulation of nickel alloys, 
C   Materials at High Temperatures, 37:5, 328-339, DOI: 10.1080/09603409.2020.1801951

      subroutine kslipDoubleExponent(xNorm,xDir,tau,signtau,tauc,burgerv,dtime,
     + nSys,iphase,CurrentTemperature,Backstress,Lp,tmat,gammaDot)
         
      implicit none
	  
	  ! number of slip system
      integer, intent(in):: nSys
	  
	  ! phase
	  integer, intent(in):: iphase

      ! slip directions and normals	  
      real*8, intent(in) :: xNorm(nSys,3),xDir(nSys,3)
	  
	  ! resolved shear stress and critical resolved shear stress
	  ! and sign of the resolved shear stress
	  real*8, intent(in) :: tau(nSys), tauc(nSys), signtau(nSys)
	  
	  ! Burgers vectors
	  real*8, intent(in) :: burgerv(nSys)
	  
	  ! Temperature
	  real*8, intent(in) :: CurrentTemperature
	  
	  ! time step
      real*8, intent(in) :: dtime
	  
	  ! Backstress state variable
	  real*8, intent(in) :: Backstress(nSys)
	  
	  ! plastic velocity gradient
      real*8, intent(out) :: Lp(3,3)
	  
	  ! and its derivative with respect to the stress
	  real*8, intent(out) :: tmat(6,6)
	  
	  ! plastic strain rate on each slip system
	  real*8, intent(out) :: gammaDot(nSys)
	 

******************************************
** The following parameters must be set **

      ! Free energy variation (J/mol)
      real*8, parameter :: dF = 286000.0
	  
	  ! Boltzmann constant (J/K)
	  real*8, parameter :: kB = 1.38e-23

**       End of parameters to set       **
******************************************

      integer :: i
      real*8  :: alpha,beta,result1, rhom,,f,T,k,gamma0,b,psi,V,
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