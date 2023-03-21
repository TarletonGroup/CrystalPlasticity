!Slip Rule Re-development
      subroutine ksvd(A,m,n,B)
      implicit none
      integer,parameter :: LWORK=5315
      real*8,parameter  :: zero = 1.0e-7
      
      integer,intent(in):: m,n
      real*8,intent(in) :: A(m,n)
      real*8,intent(out) :: B(n,m)
            
      real*8 :: U(m,n),VT(m,n),diagm(n,m),qt1(m,m),WORK(LWORK)
      integer,dimension(:),allocatable  :: IWORK
      real*8,dimension(:),allocatable :: S
      
      integer :: i,j,ialloc,INFO, LDA, LDU, LDVT
      real*8 :: xbig
      
      CHARACTER*1                  :: JOBZ
      character (len=*),parameter :: fmt2="(24(' ',(I2,1X)))",
     + fmt3="(9(' ',(F8.5,1X)))",fmt6= "(24(' ',(ES11.3,1X)))"
     
      EXTERNAL DGESDD

      !Allocate dimensions to key arrays
      !============================
      allocate(S(min(m,n)),IWORK(8*min(m,n)),STAT=ialloc)
      S=0.; B=0.; U=0.;VT=0.;diagm=0.;qt1=0.; IWORK=0
      
      !Find SVD of tau* (with LAPACK)
      !============================
      JOBZ = 'A'; LDA = m; LDU = m; LDVT = n
      
      call DGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK,
     $                   LWORK, IWORK, INFO )
      
      xbig = maxval(S)
      do i=1,ubound(S,1) 
!         if(abs(S(i)) <= zero) S(i) = 0.
         if(abs(S(i)) <= 0.0001*xbig) S(i) = 0. !i.e. difference of four orders of magnitude
      end do
      
      write(6,*)"singular values,xbig",xbig; write(6,fmt6) S
      
      !Find psuedoinverse of A
      !============================
      do i=1,ubound(S,1); if(S(i) /= 0.) diagm(i,i) = 1./S(i); end do
      qt1 = matmul(diagm,transpose(U))
      B = matmul(transpose(VT),qt1)  
             
      return
      end subroutine ksvd
