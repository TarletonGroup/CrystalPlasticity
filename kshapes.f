! *************************************************
! *       shape functions and derivatives         *
! * 20-noded reduced integration (2x2x2) element  *
! *************************************************
      subroutine kshapes(yp,yq,yr,xnat,xn,dndloc) !20-noded element
      
      implicit none      
      integer :: i
      integer,parameter:: nnodes = 20
      real*8,intent(in):: yp,yq,yr      
      real*8,intent(in):: xnat(nnodes,3)
      real*8,intent(out):: xn(nnodes),dndloc(nnodes,3)
            
      real(kind=8):: gn(nnodes),gnr(nnodes),gns(nnodes),gnt(nnodes)
     
      real(kind=8):: ri,si,ti,gr,gs,gt,dgr,dgs,dgt,xgauss
      
      character(len=*),parameter :: fmt20 = "(' ',20(F4.2,1X))"

!     Zero output arrays.
      xn = 0.; dndloc = 0.
      
      do i = 1,nnodes
       ri = xnat(i,1);si = xnat(i,2);ti = xnat(i,3)       
       !---------------------------
       if(ri == 1. .or. ri==-1.) then
         gr = 0.5*(1.+ri*yp)
         dgr = 0.5*ri
       else
         gr = (1.-yp*yp)
         dgr = -2.0*yp  
       end if         
       !---------------------------
       if(si == 1. .or. si==-1.) then
         gs = 0.5*(1.+si*yq)
         dgs = 0.5*si
       else
         gs = (1.-yq*yq)
         dgs = -2.0*yq  
       end if         
       !---------------------------
       if(ti == 1. .or. ti==-1.) then
         gt = 0.5*(1.+ti*yr)
         dgt = 0.5*ti
       else
         gt = (1.-yr*yr) 
         dgt = -2.0*yr 
       end if  
       !---------------------------
       gn(i) = gr*gs*gt
       !---------------------------
       gnr(i) = dgr*gs*gt
       gns(i) = gr*dgs*gt
       gnt(i) = gr*gs*dgt
      end do
      
!      write(6,*)"gn in shape functions"; write(6,fmt20)gn
      
!     Build the shape functions
      do i = 9,20
        xn(i) = gn(i)
      end do

!     Uel node numbering order      
      xn(1) = gn(1) - (gn(9)+gn(12)+gn(17))/2.
      xn(2) = gn(2) - (gn(9)+gn(10)+gn(18))/2.
      xn(3) = gn(3) - (gn(10)+gn(11)+gn(19))/2.
      xn(4) = gn(4) - (gn(11)+gn(12)+gn(20))/2.
      xn(5) = gn(5) - (gn(13)+gn(16)+gn(17))/2.
      xn(6) = gn(6) - (gn(13)+gn(14)+gn(18))/2.
      xn(7) = gn(7) - (gn(14)+gn(15)+gn(19))/2.
      xn(8) = gn(8) - (gn(15)+gn(16)+gn(20))/2.      
      
!     Now for derivatives!
      do i = 9,20
        dndloc(i,1) = gnr(i)
        dndloc(i,2) = gns(i)
        dndloc(i,3) = gnt(i)
      end do
      dndloc(1,1) = gnr(1) - (gnr(9)+gnr(12)+gnr(17))/2.
      dndloc(2,1) = gnr(2) - (gnr(9)+gnr(10)+gnr(18))/2.
      dndloc(3,1) = gnr(3) - (gnr(10)+gnr(11)+gnr(19))/2.
      dndloc(4,1) = gnr(4) - (gnr(11)+gnr(12)+gnr(20))/2.
      dndloc(5,1) = gnr(5) - (gnr(13)+gnr(16)+gnr(17))/2.
      dndloc(6,1) = gnr(6) - (gnr(13)+gnr(14)+gnr(18))/2.
      dndloc(7,1) = gnr(7) - (gnr(14)+gnr(15)+gnr(19))/2.
      dndloc(8,1) = gnr(8) - (gnr(15)+gnr(16)+gnr(20))/2.  
      
      dndloc(1,2) = gns(1) - (gns(9)+gns(12)+gns(17))/2.
      dndloc(2,2) = gns(2) - (gns(9)+gns(10)+gns(18))/2.
      dndloc(3,2) = gns(3) - (gns(10)+gns(11)+gns(19))/2.
      dndloc(4,2) = gns(4) - (gns(11)+gns(12)+gns(20))/2.
      dndloc(5,2) = gns(5) - (gns(13)+gns(16)+gns(17))/2.
      dndloc(6,2) = gns(6) - (gns(13)+gns(14)+gns(18))/2.
      dndloc(7,2) = gns(7) - (gns(14)+gns(15)+gns(19))/2.
      dndloc(8,2) = gns(8) - (gns(15)+gns(16)+gns(20))/2. 
      
      dndloc(1,3) = gnt(1) - (gnt(9)+gnt(12)+gnt(17))/2.
      dndloc(2,3) = gnt(2) - (gnt(9)+gnt(10)+gnt(18))/2.
      dndloc(3,3) = gnt(3) - (gnt(10)+gnt(11)+gnt(19))/2.
      dndloc(4,3) = gnt(4) - (gnt(11)+gnt(12)+gnt(20))/2.
      dndloc(5,3) = gnt(5) - (gnt(13)+gnt(16)+gnt(17))/2.
      dndloc(6,3) = gnt(6) - (gnt(13)+gnt(14)+gnt(18))/2.
      dndloc(7,3) = gnt(7) - (gnt(14)+gnt(15)+gnt(19))/2.
      dndloc(8,3) = gnt(8) - (gnt(15)+gnt(16)+gnt(20))/2.   

      return
      end subroutine

! *************************************************
! *       shape functions and derivatives         *
! * 8-noded full integration (2x2x2) element      *
! *************************************************
      subroutine kshapes8(yp,yq,yr,xnat,xn,dndloc) !8-noded element.
      
      implicit none      
      integer :: i
      integer,parameter::nnodes = 8
      real*8,intent(in):: yp,yq,yr      
      real*8,intent(in):: xnat(nnodes,3)
      real*8,intent(out):: xn(nnodes),dndloc(nnodes,3)
            
      real(kind=8):: gn(nnodes),gnr(nnodes),gns(nnodes),gnt(nnodes)
     
      real(kind=8):: ri,si,ti,gr,gs,gt,dgr,dgs,dgt,xgauss
      
      character(len=*),parameter :: fmt20 = "(' ',20(F4.2,1X))"

!     Zero output arrays.
      xn = 0.; dndloc = 0.
      
      do i = 1,nnodes !20
       ri = xnat(i,1);si = xnat(i,2);ti = xnat(i,3)       
       !---------------------------
       if(ri == 1. .or. ri==-1.) then
         gr = 0.5*(1.+ri*yp)
         dgr = 0.5*ri
       else
         gr = (1.-yp*yp)
         dgr = -2.0*yp  
       end if         
       !---------------------------
       if(si == 1. .or. si==-1.) then
         gs = 0.5*(1.+si*yq)
         dgs = 0.5*si
       else
         gs = (1.-yq*yq)
         dgs = -2.0*yq  
       end if         
       !---------------------------
       if(ti == 1. .or. ti==-1.) then
         gt = 0.5*(1.+ti*yr)
         dgt = 0.5*ti
       else
         gt = (1.-yr*yr) 
         dgt = -2.0*yr 
       end if  
       !---------------------------
       gn(i) = gr*gs*gt
       !---------------------------
       gnr(i) = dgr*gs*gt
       gns(i) = gr*dgs*gt
       gnt(i) = gr*gs*dgt
      end do

!     Uel node numbering order      
      xn(1) = gn(1)
      xn(2) = gn(2)
      xn(3) = gn(3)
      xn(4) = gn(4)
      xn(5) = gn(5)
      xn(6) = gn(6)
      xn(7) = gn(7)
      xn(8) = gn(8)
      
!     Now for derivatives!
      dndloc(1,1) = gnr(1)
      dndloc(2,1) = gnr(2)
      dndloc(3,1) = gnr(3)
      dndloc(4,1) = gnr(4)
      dndloc(5,1) = gnr(5)
      dndloc(6,1) = gnr(6)
      dndloc(7,1) = gnr(7)
      dndloc(8,1) = gnr(8)
      
      dndloc(1,2) = gns(1)
      dndloc(2,2) = gns(2)
      dndloc(3,2) = gns(3)
      dndloc(4,2) = gns(4)
      dndloc(5,2) = gns(5)
      dndloc(6,2) = gns(6)
      dndloc(7,2) = gns(7)
      dndloc(8,2) = gns(8)
      
      dndloc(1,3) = gnt(1)
      dndloc(2,3) = gnt(2)
      dndloc(3,3) = gnt(3)
      dndloc(4,3) = gnt(4)
      dndloc(5,3) = gnt(5)
      dndloc(6,3) = gnt(6)
      dndloc(7,3) = gnt(7)
      dndloc(8,3) = gnt(8)

      return
      end subroutine      
