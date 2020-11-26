!Slip Rule Re-development
      subroutine kslip0(xNorm,xDir,tau,tauc,caratio,dtime,nSys,r,iphase,
     + xlp,tmat)
      implicit none
      integer,intent(in):: nSys,iphase
      real*8,intent(in) :: dtime,r,caratio
      real*8,intent(in) :: xNorm(nSys,3),xDir(nSys,3),tau(nSys),tauc(nSys)
      real*8,intent(out) :: xlp(3,3),tmat(6,6)
      
      integer :: i
      real*8  :: xalpha,xbeta,result1,
     + xsnt(3,3),xsnv(6),xnsv(6),xsnnst(6,6),xnst(3,3),result4(6,6),
     + tempNorm(3), tempDir(3),gammaDot(nSys)
      
C
C  *** CALCULATE THE DERIVATIVE OF PLASTIC STRAIN INCREMENT WITH 
C   RESPECT TO THE STRESS DEFINED AS tmat***
C
      tmat=0.; xlp = 0.;result4=0.
	     Do I=1,nSys
         
         xalpha = 0.1; xbeta = 0.1
!         xalpha = 5.0e-7; xbeta = 5.0e-7
         !c/a ratio for HCP! !same burgers magnitude for first 12
         if(iphase==0 .and. i > 12) then 
            xalpha = caratio*caratio*xalpha
            xbeta = caratio*caratio*xbeta
         end if
C
         if (tau(I) >= tauc(I)) THEN
          gammaDot(I)=xalpha*sinh(xbeta*abs(tau(I)-r-tauc(I)))
         
          tempNorm = xNorm(I,:); tempDir = xDir(I,:)
          xsnt = spread(tempDir,2,3)*spread(tempNorm,1,3)
          xnst = spread(tempNorm,2,3)*spread(tempDir,1,3)
          CALL KGMATVEC6(xsnt,xsnv)         
          CALL KGMATVEC6(xnst,xnsv) 
          xsnnst = spread(xsnv,2,6)*spread(xnsv,1,6)
          result1 = cosh(xbeta*abs(tau(I)-r-tauc(I)))
          
          result4 = result4 + xalpha*xbeta*dtime*result1*xsnnst          
          xlp = xlp + gammaDot(I)*xsnt
         else
            gammaDot(I)=0.0
         end if
C
        END DO
      
      tmat = 0.5*(result4+transpose(result4))
             
      return
      end subroutine kslip0
