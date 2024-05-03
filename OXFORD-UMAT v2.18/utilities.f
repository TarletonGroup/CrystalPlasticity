! *************************************************
! *************************************************
! *          UTILITY SUBROUTINES                  *
! *************************************************
! *************************************************
! Updated on Sept 23rd, 2022
! Reorganized by Eralp Demir
! Converged into a module file
! Function trace is written as a subroutine
!
!
      module utilities
      implicit none
      contains
!
!
! *************************************************
! *           TRACE OF A 3X3 MATRIX               *
! *************************************************
      subroutine trace3x3(a,aii)
      implicit none
      real(8), intent(in) :: a(3,3)
      real(8), intent(out) :: aii
!      
      aii = a(1,1)+a(2,2)+a(3,3)
!
      return
      end subroutine trace3x3
!
!
!
! *************************************************
! *           TRACE OF A2X2 MATRIX               *
! *************************************************
      subroutine trace2x2(a,aii)
      implicit none
      real(8), intent(in) :: a(2,2)
      real(8), intent(out) :: aii
!
      aii = a(1,1)+a(2,2)
!
      return
      end subroutine trace2x2
!
!
! *************************************************
! *      TRANSFER 3X3 MATRIX TO 9X1 COLUMN VECTOR *
! *************************************************
      subroutine vecmat9(dvin,dmout)
      implicit none
      real(8), intent(in) :: dvin(9)
      real(8), intent(out) :: dmout(3,3)
      integer :: i
!
      dmout(1,1) = dvin(1)
      dmout(1,2) = dvin(2)
      dmout(1,3) = dvin(3)
!
      dmout(2,1) = dvin(4)
      dmout(2,2) = dvin(5)
      dmout(2,3) = dvin(6)      
!
      dmout(3,1) = dvin(7)
      dmout(3,2) = dvin(8)
      dmout(3,3) = dvin(9)
!
      return
      end subroutine vecmat9
!
!
! *************************************************
! *      TRANSFER 3X3 MATRIX TO 9X1 COLUMN VECTOR *  
! *************************************************
      subroutine matvec9(dmin,dvout)
      implicit none
      real(8), intent(in) :: dmin(3,3)
      real(8), intent(out) :: dvout(9)
      integer :: i
!
      dvout(1) = dmin(1,1)
      dvout(2) = dmin(1,2)
      dvout(3) = dmin(1,3)
!
      dvout(4) = dmin(2,1)
      dvout(5) = dmin(2,2)
      dvout(6) = dmin(2,3)
!
      dvout(7) = dmin(3,1)
      dvout(8) = dmin(3,2)
      dvout(9) = dmin(3,3)
!
      return
      end subroutine matvec9
!
!
! *************************************************
! *      TRANSFER 3X3 MATRIX TO 6X1 COLUMN VECTOR *  
! *************************************************
      subroutine matvec6(dmin,dvout)
      implicit none
      real(8), intent(in) :: dmin(3,3)
      real(8), intent(out) :: dvout(6)
      integer :: i
!
      do i=1,3
          dvout(i)=dmin(i,i)
      end do
!
      dvout(4) = (dmin(1,2)+dmin(2,1))/2.
      dvout(5) = (dmin(1,3)+dmin(3,1))/2.
      dvout(6) = (dmin(2,3)+dmin(3,2))/2.
!
      return
      end subroutine matvec6
!
!
! *************************************************
! *   TRANSFER 6X1 COLUMN VECTOR TO 3X3 MATRIX    *
! *************************************************
!
      subroutine vecmat6(dvin,dmout)
      implicit none
      real(8), intent(in) :: dvin(6)
      real(8), intent(out) :: dmout(3,3)
      integer :: i
!
      do i=1,3
            dmout(i,i) = dvin(i)
      end do
!
      dmout(1,2) = dvin(4)
      dmout(2,1) = dvin(4)
      dmout(1,3) = dvin(5)
      dmout(3,1) = dvin(5)
      dmout(3,2) = dvin(6)
      dmout(2,3) = dvin(6)
!
      return
      end subroutine vecmat6
!
!
! *************************************************
! *   VECTOR PRODUCT OF 3X1 WITH 3X1 GIVING 3X1   *
! *************************************************
      subroutine vecprod(dvin1,dvin2,dvout)
      implicit none
      real(8), intent(in) :: dvin1(3), dvin2(3)
      real(8), intent(out) :: dvout(3)
!
      dvout(1)=dvin1(2)*dvin2(3)-dvin1(3)*dvin2(2)
      dvout(2)=dvin1(3)*dvin2(1)-dvin1(1)*dvin2(3)
      dvout(3)=dvin1(1)*dvin2(2)-dvin1(2)*dvin2(1)
!
      return
      end subroutine vecprod
!
!
! *************************************************
! *   DOT PRODUCT OF 3X1 WITH 3X1 GIVING 1X1      *
! *************************************************
      subroutine dotprod(dvin1,dvin2,dvout)
      implicit none
      real(8), intent(in) :: dvin1(3), dvin2(3)
      real(8), intent(out) :: dvout(3)
!
      dvout = dvin1(1)*dvin2(1)+dvin1(2)*dvin2(2)+dvin1(3)*dvin2(3)
      dvout = abs(dvout)
!
      return
      end subroutine dotprod
!
!
! ****************************************************
! * TRANSFER GENERAL 3X3 MATRIX TO 6X1 COLUMN VECTOR *
! ****************************************************
!     Off-diagonal terms are doubled!
!     This is valid for shear conversion only!
      subroutine gmatvec6(dmin,dvout)
      implicit none
      real(8), intent(in) :: dmin(3,3)
      real(8), intent(out) :: dvout(6)
      integer :: i
!
      do i=1,3
          dvout(i)=dmin(i,i)
      end do
!
      dvout(4) = dmin(1,2)+dmin(2,1)
      dvout(5) = dmin(3,1)+dmin(1,3)
      dvout(6) = dmin(2,3)+dmin(3,2)            
!
      return
      end subroutine gmatvec6
!
! *************************************************
! *        INVERSE OF A MATRIX  WITH LAPACK       *
! *************************************************
! Checked inverse against python's numpy.linalg.inv
      subroutine lapinverse(xmatin,m,info,xmatout)
      implicit none
!
      integer,intent(in) :: m
      real(8),intent(in):: xmatin(m,m) !abaqus won't allow xmatin(:,:)
!
      integer,parameter :: lwork = 64
      real(8),parameter :: zero=1.0d-12
      integer :: i,j
!
      integer,intent(out) :: info
      real(8),intent(out) :: xmatout(m,m)
      integer ::  ipiv(m)
      real(8)  :: a(m,m)
      real(8) :: work(lwork)
!
!
! https://software.intel.com/en-us/mkl-developer-reference-fortran-getri#626EB2AE-CA6A-4233-A6FA-04F54EF7A6E6
!
!      EXTERNAL DGETRI, DGETRF
!
!
!
      xmatout = 0.
!     ED: I have added to run this
!     The next line shall be removed
      info=1
!
      a = xmatin !don't input array xmatin
!
!      call DGETRF( m, m, a, m, ipiv, info )
!
      if(info == 0) then
!          call DGETRI( m, a, m, ipiv, work, lwork, info )
          !write(*,*) "work(1) == min lwork needed", work(1)
      else
          xmatout = 0.
!          write(*,*)"dgetrf, illegal value at = ",-info,". no inverse"         
      end if
!
      do i=1,m;
          do j=1,m;
              if (abs(a(i,j))<= zero) a(i,j) = 0. 
          end do 
      end do
      xmatout = a
!
!
!      
!
      return
      end subroutine lapinverse
!
!
!
!
! ****************************************************
! *        INVERSE OF A MATRIX  WITHOUT LAPACK       *
! ****************************************************
!
      subroutine nolapinverse(ain,c,n)
!     ============================================================
!     Inverse matrix
!     Method: Based on Doolittle LU factorization for Ax=b
!     Alex G. December 2009
!     -----------------------------------------------------------
!     input ...
!     a(n,n) - array of coefficients for matrix A
!     n      - dimension
!     output ...
!     c(n,n) - inverse matrix of A
!
!     ===========================================================
      implicit none
!
      integer, intent(in) :: n
      real(8), intent(out) :: c(n,n)
      real(8), intent(in) :: ain(n,n)
      real(8) :: L(n,n), U(n,n), b(n), d(n), x(n)
      real(8) :: coeff, a(n,n)
      integer :: i, j, k
!
!     step 0: initialization for matrices L and U and b
!     Fortran 90/95 allows such operations on matrices
      L=0.
      U=0.
      b=0.
      a=ain
!
!     step 1: forward elimination
      do k=1,n-1
          do i=k+1,n
              coeff=a(i,k)/a(k,k)
              L(i,k) = coeff
              do j=k+1,n
                  a(i,j) = a(i,j)-coeff*a(k,j)
              end do
          end do
      end do
!
!     Step 2: prepare L and U matrices 
!     L matrix is a matrix of the elimination coefficient
!     + the diagonal elements are 1.0
      do i=1,n
          L(i,i) = 1.0
      end do  
!     U matrix is the upper triangular part of A
      do j=1,n
          do i=1,j
              U(i,j) = a(i,j)
          end do
      end do
!
!     Step 3: compute columns of the inverse matrix C
      do k=1,n
          b(k)=1.0
          d(1) = b(1)
!     Step 3a: Solve Ld=b using the forward substitution
          do i=2,n
              d(i)=b(i)
              do j=1,i-1
                  d(i) = d(i) - L(i,j)*d(j)
              end do
          end do
!         Step 3b: Solve Ux=d using the back substitution
          x(n)=d(n)/U(n,n)
          do i = n-1,1,-1
              x(i) = d(i)
              do j=n,i+1,-1
                  x(i)=x(i)-U(i,j)*x(j)
              end do
              x(i) = x(i)/u(i,i)
          end do
!         Step 3c: fill the solutions x(n) into column k of C
          do i=1,n
              c(i,k) = x(i)
          end do
          b(k)=0.0
      end do
!
      end subroutine nolapinverse
!
!
!
!
!
! *************************************************
! *          SUBROUTINE MATRIX SQUARE ROOT        *
! *************************************************
!
      subroutine msqrt(a,b)
      implicit none
      real(8), intent(inout) :: a(3,3)
      real(8), intent(out) :: b(3,3)
      integer :: i
      real(8) :: diag(3,3), q(3,3), d(3),
     + qtrans(3,3), res(3,3)
!
!
      diag=0.; b=0.; q=0.;res=0.;d=0.
!             
      call jacobi(a,3,d,q) 
      call eigsrt(d,q,3)
!
      do i=1,3
         if (d(i) .ge. 0) then
            diag(i,i)=sqrt(d(i))
         else
            write (6,*) 'the matrix is not positive definite'
         end if
      end do
!
      res = matmul(q,diag)
      b = matmul(res,transpose(q))
!
      return
      end subroutine msqrt
!
!
!
! *************************************************
! * SUBROUTINE MATRIX EIGENVALUES AND EIGENVECTORS*
! *************************************************
!
      subroutine jacobi(a,n,d,v)
      implicit none
      integer, intent(in) :: n
      real(8), intent(inout) :: a(n,n)
      real(8), intent(out) :: d(n), v(n,n)
!
      integer, parameter :: nmax=500
      integer :: i, j, ip, iq, nrot
      real(8) :: b(nmax), z(nmax), dial(n,n),
     + sm, tresh, g, h, t, theta, c, s, ta
!
!
      do ip=1,n
         do iq=1,n
            v(ip,iq)=0.
         end do
         v(ip,ip)=1.
      end do
!
      do ip=1,n
         b(ip)=a(ip,ip)
         d(ip)=b(ip)
         z(ip)=0.
      end do
!
      nrot=0
      do i=1,50
         sm=0.
         do ip=1,n-1
            do iq=ip+1,n
               sm=sm+abs(a(ip,iq))
            end do
         end do
!
         if (sm .eq. 0.) return
         if (i .lt. 4) then
           tresh=0.2*sm/n**2
         else
           tresh=0.
         end if
!
         do ip=1,n-1
           do iq=ip+1,n
              g=100.*abs(a(ip,iq))
              if ((i .gt. 4) .and. (abs(d(ip))+g .eq. abs(d(ip)))
     + .and. (abs(d(ip))+g .eq. abs(d(iq)))) then
                a(ip,iq)=0.
              else if (abs(a(ip,iq)) .gt. tresh) then
                h=d(iq)-d(ip)
                if (abs(h)+g .eq. abs(h)) then
                   t=a(ip,iq)/h 
                else
                   theta=0.5*h/a(ip,iq)
                   t=1./(abs(theta)+sqrt(1.+theta**2))
                   if (theta .lt. 0) t=-t
                end if
                c=1./sqrt(1.+t**2)
                s=t*c
                ta=s/(1.+c)
                h=t*a(ip,iq)
                z(ip)=z(ip)-h
                z(iq)=z(iq)+h
                d(ip)=d(ip)-h
                d(iq)=d(iq)+h
                a(ip,iq)=0.
!
                do j=1,ip-1
                  g=a(j,ip)
                  h=a(j,iq)
                  a(j,ip)=g-s*(h+g*ta)
                  a(j,iq)=h+s*(g-h*ta)
                end do
!
                do j=ip+1,iq-1
                 g=a(ip,j)
                 h=a(j,iq)
                 a(ip,j)=g-s*(h+g*ta)
                 a(j,iq)=h+s*(g-h*ta)
                end do
!
                do j=iq+1,n
                 g=a(ip,j)
                 h=a(iq,j)
                 a(ip,j)=g-s*(h+g*ta)
                 a(iq,j)=h+s*(g-h*ta)
                end do
!
                do j=1,n
                  g=v(j,ip)
                  h=v(j,iq)
                  v(j,ip)=g-s*(h+g*ta)
                  v(j,iq)=h+s*(g-h*ta)
                end do
!
                nrot=nrot+1
              end if
           end do
         end do
!
         do ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip)
            z(ip)=0.
         end do
      end do
!      
!
      return
      end subroutine jacobi
!
!
! *************************************************
! *          SUBROUTINE SORT EIGENVALUES          *
! *************************************************
!
      subroutine eigsrt(d,v,n)
      implicit none
      integer, intent(in) :: n
      real(8), intent(inout) :: d(n)
      real(8), intent(inout) :: v(n,n)
!
      integer :: i, j, k
      real(8) :: p
!
      do i=1,n-1
        k=i
        p=d(i)
        do j=i+1,n
          if(d(j).ge.p)then
            k=j
            p=d(j)
          endif
        end do
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
          end do
        endif
      end do
!
      return
      end subroutine eigsrt
!
!
! *************************************************
! *      THE DETERMINANT OF A 3X3 MATRIX          *
! *************************************************
      subroutine deter3x3(dmin,d)
      implicit none
      real(8), intent(in) :: dmin(3,3)
      real(8), intent(out) :: d
!
      d = 0.
      d = dmin(1,1)*dmin(2,2)*dmin(3,3) + 
     + dmin(1,2)*dmin(2,3)*dmin(3,1) + 
     + dmin(2,1)*dmin(3,2)*dmin(1,3) -
     + dmin(1,3)*dmin(2,2)*dmin(3,1) -
     + dmin(1,1)*dmin(2,3)*dmin(3,2) -
     + dmin(1,2)*dmin(2,1)*dmin(3,3)
!
      return
      end subroutine deter3x3
!
!
! *************************************************
! *      THE DETERMINANT OF A 2X2 MATRIX          *
! *************************************************
      subroutine deter2x2(dmin,d)
      implicit none
      real(8), intent(in) :: dmin(2,2)
      real(8), intent(out) :: d
!
      d=0.
      d=dmin(1,1)*dmin(2,2)-dmin(1,2)*dmin(2,1)
!
      return
      end subroutine deter2x2
!
! *************************************************
! *      Build 4th order rotation matrix (tsigma) *
! *************************************************
!     valid for rotating symmetric tensors
      subroutine rotord4sig(xrot,tsigma) 
      implicit none
      real(8), intent(in) :: xrot(3,3)
      real(8), intent(out) :: tsigma(6,6)
!
      tsigma(1,1) = xrot(1,1)*xrot(1,1)
      tsigma(1,2) = xrot(1,2)*xrot(1,2)
      tsigma(1,3) = xrot(1,3)*xrot(1,3)
      tsigma(1,4) = 2.*xrot(1,1)*xrot(1,2)
      tsigma(1,6) = 2.*xrot(1,2)*xrot(1,3)
      tsigma(1,5) = 2.*xrot(1,3)*xrot(1,1)
!
      tsigma(2,1) = xrot(2,1)*xrot(2,1)
      tsigma(2,2) = xrot(2,2)*xrot(2,2)
      tsigma(2,3) = xrot(2,3)*xrot(2,3)
      tsigma(2,4) = 2.*xrot(2,1)*xrot(2,2)
      tsigma(2,6) = 2.*xrot(2,2)*xrot(2,3)
      tsigma(2,5) = 2.*xrot(2,3)*xrot(2,1)
!
      tsigma(3,1) = xrot(3,1)*xrot(3,1)
      tsigma(3,2) = xrot(3,2)*xrot(3,2)
      tsigma(3,3) = xrot(3,3)*xrot(3,3)
      tsigma(3,4) = 2.*xrot(3,1)*xrot(3,2)
      tsigma(3,6) = 2.*xrot(3,2)*xrot(3,3)
      tsigma(3,5) = 2.*xrot(3,3)*xrot(3,1)
!
      tsigma(4,1) = xrot(1,1)*xrot(2,1)
      tsigma(4,2) = xrot(1,2)*xrot(2,2)
      tsigma(4,3) = xrot(1,3)*xrot(2,3)
      tsigma(4,4) = xrot(1,1)*xrot(2,2) + xrot(1,2)*xrot(2,1)
      tsigma(4,6) = xrot(1,2)*xrot(2,3) + xrot(2,2)*xrot(1,3) 
      tsigma(4,5) = xrot(1,3)*xrot(2,1) + xrot(2,3)*xrot(1,1)
!
      tsigma(6,1) = xrot(2,1)*xrot(3,1)
      tsigma(6,2) = xrot(2,2)*xrot(3,2)
      tsigma(6,3) = xrot(2,3)*xrot(3,3)
      tsigma(6,4) = xrot(2,1)*xrot(3,2) + xrot(3,1)*xrot(2,2)
      tsigma(6,6) = xrot(2,2)*xrot(3,3) + xrot(2,3)*xrot(3,2) 
      tsigma(6,5) = xrot(2,3)*xrot(3,1) + xrot(3,3)*xrot(2,1)
!
      tsigma(5,1) = xrot(3,1)*xrot(1,1)
      tsigma(5,2) = xrot(3,2)*xrot(1,2)
      tsigma(5,3) = xrot(3,3)*xrot(1,3)
      tsigma(5,4) = xrot(3,1)*xrot(1,2) + xrot(1,1)*xrot(3,2)
      tsigma(5,6) = xrot(3,2)*xrot(1,3) + xrot(1,2)*xrot(3,3) 
      tsigma(5,5) = xrot(3,3)*xrot(1,1) + xrot(3,1)*xrot(1,3)
!      
      return
      end subroutine rotord4sig
!
!
! *************************************************
! *      Build 4th order rotation matrix (tstran) *
! *************************************************
!     valid for symmetric tensors
      subroutine rotord4str(xrot,tstran)
      implicit none
      real(8), intent(in) :: xrot(3,3)
      real(8), intent(out) :: tstran(6,6)
!  
      tstran(1,1) = xrot(1,1)*xrot(1,1)
      tstran(1,2) = xrot(1,2)*xrot(1,2)
      tstran(1,3) = xrot(1,3)*xrot(1,3)
      tstran(1,4) = xrot(1,1)*xrot(1,2)
      tstran(1,6) = xrot(1,2)*xrot(1,3)
      tstran(1,5) = xrot(1,3)*xrot(1,1)
!
      tstran(2,1) = xrot(2,1)*xrot(2,1)
      tstran(2,2) = xrot(2,2)*xrot(2,2)
      tstran(2,3) = xrot(2,3)*xrot(2,3)
      tstran(2,4) = xrot(2,1)*xrot(2,2)
      tstran(2,6) = xrot(2,2)*xrot(2,3)
      tstran(2,5) = xrot(2,3)*xrot(2,1)
!
      tstran(3,1) = xrot(3,1)*xrot(3,1)
      tstran(3,2) = xrot(3,2)*xrot(3,2)
      tstran(3,3) = xrot(3,3)*xrot(3,3)
      tstran(3,4) = xrot(3,1)*xrot(3,2)
      tstran(3,6) = xrot(3,2)*xrot(3,3)
      tstran(3,5) = xrot(3,3)*xrot(3,1)
!
      tstran(4,1) = 2.*xrot(1,1)*xrot(2,1)
      tstran(4,2) = 2.*xrot(1,2)*xrot(2,2)
      tstran(4,3) = 2.*xrot(1,3)*xrot(2,3)
      tstran(4,4) = xrot(1,1)*xrot(2,2) + xrot(1,2)*xrot(2,1)
      tstran(4,6) = xrot(1,2)*xrot(2,3) + xrot(2,2)*xrot(1,3) 
      tstran(4,5) = xrot(1,3)*xrot(2,1) + xrot(2,3)*xrot(1,1)
!
      tstran(6,1) = 2.*xrot(2,1)*xrot(3,1)
      tstran(6,2) = 2.*xrot(2,2)*xrot(3,2)
      tstran(6,3) = 2.*xrot(2,3)*xrot(3,3)
      tstran(6,4) = xrot(2,1)*xrot(3,2) + xrot(3,1)*xrot(2,2)
      tstran(6,6) = xrot(2,2)*xrot(3,3) + xrot(2,3)*xrot(3,2) 
      tstran(6,5) = xrot(2,3)*xrot(3,1) + xrot(3,3)*xrot(2,1)
!
      tstran(5,1) = 2.*xrot(3,1)*xrot(1,1)
      tstran(5,2) = 2.*xrot(3,2)*xrot(1,2)
      tstran(5,3) = 2.*xrot(3,3)*xrot(1,3)
      tstran(5,4) = xrot(3,1)*xrot(1,2) + xrot(1,1)*xrot(3,2)
      tstran(5,6) = xrot(3,2)*xrot(1,3) + xrot(1,2)*xrot(3,3) 
      tstran(5,5) = xrot(3,3)*xrot(1,1) + xrot(3,1)*xrot(1,3)
!
      return
      end subroutine rotord4str
!
!
! ***************************************************************
! *   Build 4th order lower (and upper) tensor product          *
! *   NB: When contracted from 4th to the format used for       *
! *   C-matrix, it turns out that the upper and lower products  *
! *   are the same!                                             *
! ***************************************************************
!     valid for symmetric tensors
      subroutine ltprod(a,b,c)
      implicit none
      real(8), intent(in) :: a(3,3), b(3,3)
      real(8), intent(out) :: c(6,6)
      integer :: i, j
!
      c = 0.
      c(1,1) = 2.*a(1,1)*b(1,1)
      c(1,2) = 2.*a(1,2)*b(1,2)
      c(1,3) = 2.*a(1,3)*b(1,3)
      c(1,4) = a(1,1)*b(1,2)+a(1,2)*b(1,1)
      c(1,5) = a(1,1)*b(1,3)+a(1,3)*b(1,1)
      c(1,6) = a(1,2)*b(1,3)+a(1,3)*b(1,2)
!
      c(2,2) = 2.*a(2,2)*b(2,2)
      c(2,3) = 2.*a(2,3)*b(2,3)
      c(2,4) = a(2,1)*b(2,2)+a(2,2)*b(2,1)
      c(2,5) = a(2,1)*b(2,3)+a(2,3)*b(2,1)
      c(2,6) = a(2,2)*b(2,3)+a(2,3)*b(2,2)
!
      c(3,3) = 2.*a(3,3)*b(3,3)
      c(3,4) = a(3,1)*b(3,2)+a(3,2)*b(3,1)
      c(3,5) = a(3,1)*b(3,3)+a(3,3)*b(3,1)
      c(3,6) = a(3,2)*b(3,3)+a(3,3)*b(3,2)
!
      c(4,4) = a(1,1)*b(2,2)+a(1,2)*b(2,1)
      c(4,5) = a(1,1)*b(2,3)+a(1,3)*b(2,1)
      c(4,6) = a(1,2)*b(2,3)+a(1,3)*b(2,2)
!
      c(5,5) = a(1,1)*b(3,3)+a(1,3)*b(3,1)
      c(5,6) = a(1,2)*b(3,3)+a(1,3)*b(3,2)
!
      c(6,6) = a(2,2)*b(3,3)+a(2,3)*b(3,2)
!
      do i=2,6
         do j=1,i-1
            c(i,j)=c(j,i)
         end do
      end do
!
      c = 0.5*c
!
      return
      end subroutine ltprod
!
!
! **************************************************
! *      Build 4th order tensor product (kronecker)*
! **************************************************
!     valid for symmetric tensors
      subroutine tprod(a,b,c)
      implicit none
      real(8), intent(in) :: a(3,3), b(3,3)
      real(8), intent(out) :: c(6,6)
      integer :: i, j
!
      c = 0.
      c(1,1) = 2.*a(1,1)*b(1,1)
      c(1,2) = 2.*a(1,1)*b(2,2)
      c(1,3) = 2.*a(1,1)*b(3,3)
      c(1,4) = a(1,1)*b(1,2)+a(1,1)*b(2,1)
      c(1,5) = a(1,1)*b(1,3)+a(1,1)*b(3,1)
      c(1,6) = a(1,1)*b(2,3)+a(1,1)*b(3,2)
!
      c(2,2) = 2.*a(2,2)*b(2,2)
      c(2,3) = 2.*a(2,2)*b(3,3)
      c(2,4) = a(2,2)*b(1,2)+a(2,2)*b(2,1)
      c(2,5) = a(2,2)*b(1,3)+a(2,2)*b(3,1)
      c(2,6) = a(2,2)*b(2,3)+a(2,2)*b(3,2)
!
      c(3,3) = 2.*a(3,3)*b(3,3)
      c(3,4) = a(3,3)*b(1,2)+a(3,3)*b(2,1)
      c(3,5) = a(3,3)*b(1,3)+a(3,3)*b(3,1)
      c(3,6) = a(3,3)*b(2,3)+a(3,3)*b(3,2)
!
      c(4,4) = a(1,2)*b(1,2)+a(1,2)*b(2,1)
      c(4,5) = a(1,2)*b(1,3)+a(1,2)*b(3,1)
      c(4,6) = a(1,2)*b(2,3)+a(1,2)*b(3,2)
!
      c(5,5) = a(1,3)*b(1,3)+a(1,3)*b(3,1)
      c(5,6) = a(1,3)*b(2,3)+a(1,3)*b(3,2)
! 
      c(6,6) = a(2,3)*b(2,3)+a(2,3)*b(3,2)
!
      do i=2,6
         do j=1,i-1
            c(i,j)=c(j,i)
         end do
      end do
!
      c = 0.5*c
!  
      return
      end subroutine tprod
!
!
!
! **************************************
! *         MULTIPLY 3x3 MATRICES      *
! **************************************
      subroutine mmult(dm1,dm2,dm)
      implicit none
      real(8), intent(in) :: dm1(3,3), dm2(3,3)
      real(8), intent(out) :: dm(3,3)
      integer :: i, j, k
      real(8) :: x
!
!
      do i=1,3
          do j=1,3
              x=0.0
              do k=1,3
                  x=x+dm1(i,k)*dm2(k,j)
              end do
              dm(i,j)=x
          end do
      end do
!
      return
      end subroutine mmult
!
!
!	This subroutine inverts a 3x3 matrix
!	INPUT:	Matrix								---	A(3,3)
!	OUTPUT:	Invereted matrix, determinant		---	invA(3,3),det
	subroutine inv3x3(A,invA,det)
      use globalvariables, only: smallnum
	implicit none
      real(8), intent(in)  :: A(3,3)
      real(8), intent(out) :: invA(3,3), det
	integer :: i,j
!
!
!	First calculate the determinant
	call deter3x3(A,det)
!	If the determinant is greater than certain value
	if (abs(det) < smallnum) then
		invA=0.0d+0
	else
		invA(1,1)=((A(2,2)*A(3,3))-(A(2,3)*A(3,2)))/det
		invA(2,1)=-((A(2,1)*A(3,3))-(A(2,3)*A(3,1)))/det
		invA(3,1)=((A(2,1)*A(3,2))-(A(2,2)*A(3,1)))/det
		invA(1,2)=-((A(1,2)*A(3,3))-(A(1,3)*A(3,2)))/det
		invA(2,2)=((A(1,1)*A(3,3))-(A(1,3)*A(3,1)))/det
		invA(3,2)=-((A(1,1)*A(3,2))-(A(1,2)*A(3,1)))/det
		invA(1,3)=((A(1,2)*A(2,3))-(A(1,3)*A(2,2)))/det
		invA(2,3)=-((A(1,1)*A(2,3))-(A(2,1)*A(1,3)))/det
		invA(3,3)=((A(1,1)*A(2,2))-(A(1,2)*A(2,1)))/det
	endif
	return
      end subroutine inv3x3
!
!
!	This subroutine inverts a 2x2 matrix
!	INPUT:	Matrix								---	A(2,2)
!	OUTPUT:	Invereted matrix, determinant		---	invA(2,2),det
	subroutine inv2x2(A,invA,det)
      use globalvariables, only: smallnum
	implicit none
      real(8), intent(in)  :: A(2,2)
      real(8), intent(out) :: invA(2,2), det
	integer :: i,j
!
!
!	First calculate the determinant
	call deter2x2(A,det)
!	If the determinant is greater than certain value
	if (abs(det) < smallnum) then
		invA=0.0d+0
	else
		invA(1,1) = A(2,2)/det
		invA(1,2) =-A(1,2)/det
		invA(2,1) =-A(2,1)/det
          invA(2,2) = A(1,1)/det
	endif
	return
	end subroutine inv2x2
!
!	Euler angles to crystal orientation matrix
!	INPUT:	Angles(deg)			---	ang(3)
!	OUTPUT:	Orientation matrix	---	R(3,3)
!     USES:     Number pi           --- pi
	subroutine Euler2ori(Euler,R)
	use globalvariables, only : pi
	implicit none
	real(8), intent(in) :: Euler(3)
      real(8), intent(out) :: R(3,3)
	real(8) :: phi1, phi2, Phi
!
!     convert to radians
	phi1=Euler(1)*pi/180.
      Phi=Euler(2)*pi/180.
	phi2=Euler(3)*pi/180.
!
!
      R=0.
      R(1,1)=(cos(phi1)*cos(phi2))-(sin(phi1)*sin(phi2)*cos(Phi))
      R(2,1)=-(cos(phi1)*sin(phi2))-(sin(phi1)*cos(phi2)*cos(Phi))
      R(3,1)=sin(phi1)*sin(Phi)
      R(1,2)=(sin(phi1)*cos(phi2))+(cos(phi1)*sin(phi2)*cos(Phi))
      R(2,2)=-(sin(phi1)*sin(phi2))+(cos(phi1)*cos(phi2)*cos(Phi))
      R(3,2)=-cos(phi1)*sin(Phi)
      R(1,3)=sin(phi2)*sin(Phi)
      R(2,3)=cos(phi2)*sin(Phi)
      R(3,3)=cos(Phi)
!
	return
      end subroutine Euler2ori
!
!
!	Orientation matrix from Euler angles
!	INPUT:	Orientation matrix	---	R(3)
!	OUTPUT:	Bunge angles(deg)   ---	ang(3)
!     USES:     Number pi         --- pi
	subroutine ori2Euler(R,Euler)
      use globalvariables, only : pi
	implicit none
      real(8), intent(in) :: R(3,3)
      real(8), intent(out) :: Euler(3)
	real(8) :: phi1, phi2, Phi
!
	if (R(3,3).eq.1.) then
		Phi=0.
		phi1=atan2(R(1,2),R(1,1))
		phi2=0.
	else
		Phi=acos(R(3,3))
		phi1=atan2(R(3,1)/sin(Phi),-R(3,2)/sin(Phi))
		phi2=atan2(R(1,3)/sin(Phi),R(2,3)/sin(Phi))
      endif
      Euler(1)=phi1
      Euler(2)=Phi
      Euler(3)=phi2
!     convert to degrees
	Euler=Euler*180./pi
!
	return
      end subroutine ori2Euler   
!
!
!     This subroutine calculates inverse of a square matrix
!     by singular value decomposition
      subroutine SVDinverse(A,n,invA,err)
      use globalvariables, only: smallnum
      implicit none
      integer n
      real(8), intent(in) :: A(n,n)
      real(8), intent(out) :: invA(n,n)
      integer, intent(out) :: err
!     local variables
      real(8) :: w(n), V(n,n), U(n,n)
      real(8) :: invS(n,n), tol
      integer :: i
!
!
!     tolerance value
      tol = sqrt(smallnum)
!
      U = A
!     Singular Value Decomposition
      call svdcmp(U,n,n,n,n,w,V,err)
!
!     subroutine returns
!     U = A
!     V = V
!     S(i,i) = w(i)
!
!     Calculate the inverse
!     A = U * S * V^T
!     invA = V * inv(S) * U^T
!     inv(S) is the pseudo inverse 1/w(i)
!
!     If there is no error
      if (err==0) then
!    
          invS = 0.
          do i = 1, n
              if (abs(w(i)) > tol) then
                  invS(i,i) = 1. / w(i)
              end if
          end do
!         Corrected by Alvaro
          invA = matmul(matmul(V,invS),transpose(U))
!     In case of an error
      else
!
          V = 0.
          invA = 0.
!
!
      end if
!
      end subroutine SVDinverse
!
!
!     Generalized inverse by singular value decomposition
!     Written by Alvaro Martinez 08-03-2023
!
      subroutine SVDgeninverse(A,n,m,invA,err)
      use globalvariables, only: smallnum
      implicit none
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(8), intent(in) :: A(n,m)
      real(8), intent(out) :: invA(m,n)
      integer, intent(out) :: err
!     local variables
      real(8) :: w(m), V(m,m), U(n,m)
      real(8) :: invS(m,m), tol
      integer :: i
!
!
!     tolerance value
      tol = sqrt(smallnum)
!
      U = A
!     Singular Value Decomposition
      call svdcmp(U,n,m,n,m,w,V,err)
!     subroutine returns
!     U = A
!     V = V
!     S(i,i) = w(i)
!
!     Calculate the inverse
!     A = U * S * V^T
!     invA = V * inv(S) * U^T
!     inv(S) is the pseudo inverse 1/w(i)
!
!     If there is no error
      if (err==0) then
          invS = 0.
          do i = 1, m
              if (abs(w(i)) > tol) then
                  invS(i,i) = 1. / w(i)
              end if
          end do
!
          invA = matmul(matmul(V,invS),transpose(U))
!     In case of an error
      else
!
          V = 0.
          invA = 0.
!
!
      end if
!
      end subroutine SVDgeninverse
!
!
!
!     Codes for singular value decomposition
!     Numerical Recipies in F77
!     https://websites.pmc.ucsc.edu/~fnimmo/eart290c_17/NumericalRecipesinF77.pdf
!  
!     SVDcmp subroutine
!     Given a matrix (1:m,1:n) with physical dimensions mp by np,
!     this routine computes its singular value decomposition,
!     A = U W VT.  The matrix U replaces A on output.  The diagonal
!     matrix of singular values W is output as a vector w(1:n)
!     The matrix V (not the transpose VT) is the output as V(1:n,1:n) 
!
      subroutine svdcmp(A,m,n,mp,np,w,V,err)
      use errors, only: error
      implicit none
      integer, intent(in) :: m,mp,n,np
      real(8), intent(inout) :: A(mp,np)
      real(8), intent(out) :: V(np,np), w(np)
      integer, intent(out) :: err
      integer, parameter :: nmax=1000
!     uses pythag
      integer i,its,j,jj,k,l,nm
      real(8) anorm,c,f,g,h,s,scale,x,y,z,rv1(nmax),pyt
!     initialize error flag
      err=0
!
      g=0.0
      scale=0.0
      anorm=0.0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0
        s=0.0
        scale=0.0
        if(i.le.m)then
          do 11 k=i,m
            scale=scale+abs(A(k,i))
11        continue
          if(scale.ne.0.0)then
            do 12 k=i,m
              A(k,i)=A(k,i)/scale
              s=s+A(k,i)*A(k,i)
12          continue
            f=A(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            A(i,i)=f-g
            do 15 j=l,n
              s=0.0
              do 13 k=i,m
                s=s+A(k,i)*A(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                A(k,j)=A(k,j)+f*A(k,i)
14            continue
15          continue
            do 16 k=i,m
              A(k,i)=scale*A(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0
        s=0.0
        scale=0.0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scale=scale+abs(A(i,k))
17        continue
          if(scale.ne.0.0)then
            do 18 k=l,n
              A(i,k)=A(i,k)/scale
              s=s+A(i,k)*A(i,k)
18          continue
            f=A(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            A(i,l)=f-g
            do 19 k=l,n
              rv1(k)=A(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0
              do 21 k=l,n
                s=s+A(j,k)*A(i,k)
21            continue
              do 22 k=l,n
                A(j,k)=A(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n
              A(i,k)=scale*A(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.0)then
            do 26 j=l,n
              V(j,i)=(A(i,j)/A(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0
              do 27 k=l,n
                s=s+A(i,k)*V(k,j)
27            continue
              do 28 k=l,n
                V(k,j)=V(k,j)+s*V(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            V(i,j)=0.0
            V(j,i)=0.0
31        continue
        endif
        V(i,i)=1.0
        g=rv1(i)
        l=i
32    continue
      do 39 i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do 33 j=l,n
          A(i,j)=0.0
33      continue
        if(g.ne.0.0)then
          g=1.0/g
          do 36 j=l,n
            s=0.0
            do 34 k=l,m
              s=s+A(k,i)*A(k,j)
34          continue
            f=(s/A(i,i))*g
            do 35 k=i,m
              A(k,j)=A(k,j)+f*A(k,i)
35          continue
36        continue
          do 37 j=i,m
            A(j,i)=A(j,i)*g
37        continue
        else
          do 38 j= i,m
            A(j,i)=0.0
38        continue
        endif
        A(i,i)=A(i,i)+1.0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0
          s=1.0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            call pythag(f,g,h)
            w(i)=h
            h=1.0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=A(j,nm)
              z=A(j,i)
              A(j,nm)=(y*c)+(z*s)
              A(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=w(k)
          if(l.eq.k)then
            if(z.lt.0.0)then
              w(k)=-z
              do 44 j=1,n
                V(j,k)=-V(j,k)
44            continue
            endif
            goto 3
          endif
          if(its.eq.30) then !pause 'no convergence in svdcmp'
              err=1
          endif
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          call pythag(f,1.0d+0,g)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            call pythag(f,h,z)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=V(jj,j)
              z=V(jj,i)
              V(jj,j)= (x*c)+(z*s)
              V(jj,i)=-(x*s)+(z*c)
45          continue
            call pythag(f,h,z)
            w(j)=z
            if(z.ne.0.0)then
              z=1.0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=A(jj,j)
              z=A(jj,i)
              A(jj,j)= (y*c)+(z*s)
              A(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      end subroutine svdcmp
!     (C) Copr. 1986-92 Numerical Recipes Software
!
!
!
!
!
!
      subroutine pythag(a,b,pyt)
      implicit none
      real(8), intent(in):: a,b
      real(8), intent(out):: pyt
      real(8) absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pyt=absa*sqrt(1.+(absb/absa)**2)
      else
        if(absb.eq.0.)then
          pyt=0.
        else
          pyt=absb*sqrt(1.+(absa/absb)**2)
        endif
      endif
      return
      end subroutine pythag
!     (C) Copr. 1986-92 Numerical Recipes Software
!
!
!
!
      end module utilities
