      subroutine kdirns(xrot,TwinRot,iphase,L,nTwin,
     + ddir1,dnor1,dtwindir1,dtwinnor1,
     + caratio,cubicslip)
      
      implicit none
      integer :: i,j,k,slip
      integer,parameter :: m = 3
      integer,intent(in):: iphase,L,nTwin
      
      real*8,intent(in):: xrot(m,m), TwinRot(nTwin,m,m), caratio
	  
	  ! activate cubic slip systems for CMSX-4 alloy
      INTEGER,intent(in) :: cubicslip
	  
      real*8,intent(out) :: ddir1(L,m),dnor1(L,m)
	  real*8,intent(out) :: dtwindir1(nTwin,m),dtwinnor1(nTwin,m)

      real(kind=8) :: ddir(L,m),dnor(L,m)
	  real(kind=8) :: dtwindir(nTwin,m),dtwinnor(nTwin,m)
      real(kind=8) :: tdir(m),tnor(m),ttwindir(m),ttwinnor(m)
	  real(kind=8) :: tdir1(m),tnor1(m),ttwindir1(m),ttwinnor1(m)
     
      real(kind=8):: xdirmag,xnormag,xtwindirmag,xtwinnormag

         ddir1 = 0.0; dnor1 = 0.0 
      
      if(iphase .eq. 0) then !hcp
         include 'xDir0ET.f'
         include 'xNorm0ET.f'         
      else if(iphase .eq. 1) then !bcc (24/48)
         include 'xDir1.f'
         include 'xNorm1.f'   
      else if(iphase .eq. 2) then !fcc (12)
         if (cubicslip == 0) then ! cubic slip
            include 'xDir2.f'
            include 'xNorm2.f'
         else
            include 'xDir2c.f'
            include 'xNorm2c.f' 
         end if 
      else if (iphase .eq. 4) then ! olivine (7)
          include 'xDir4.f'
          include 'xNorm4.f'
      else if (iphase .eq. 5) then ! alpha-Uranium (8)
          include 'xDirAlphaUranium.f'
          include 'xNormAlphaUranium.f'
          include 'xTwinDirAlphaUranium.f'
          include 'xTwinNormAlphaUranium.f'
      end if     

      do k=1,L ! rotate slip directions (without twinning)
        do i=1,m
        tdir(i) = ddir(k,i)
        tnor(i) = dnor(k,i)
        end do
        
        if(iphase .eq. 0) then !hcp ET added
          tdir(3) =  tdir(3)*caratio
          tnor(3) =  tnor(3)/caratio
        end if  
        
        tdir1 = matmul(xrot,tdir)
        tnor1 = matmul(xrot,tnor)
        
        xdirmag = sqrt(dot_product(tdir1,tdir1))
        xnormag = sqrt(dot_product(tnor1,tnor1))
        
        do i=1,m
        ddir1(k,i) = tdir1(i)/xdirmag
        dnor1(k,i) = tnor1(i)/xnormag
        end do
      end do

      if (iphase .eq. 5) then ! alpha-Uranium (2 twin systems)

        do k=1,nTwin ! rotate twin direction and normal with gmatinv
          do i=1,m
            ttwindir(i) = dtwindir(k,i)
            ttwinnor(i) = dtwinnor(k,i)
          end do

          ttwindir1 = matmul(xrot,ttwindir) ! rotation with gmatinv
          ttwinnor1 = matmul(xrot,ttwinnor)

          xtwindirmag = sqrt(dot_product(ttwindir1,ttwindir1))
          xtwinnormag = sqrt(dot_product(ttwinnor1,ttwinnor1))

          do i=1,m ! assign value to output
            dtwindir1(k,i) = ttwindir1(i)/xtwindirmag
            dtwinnor1(k,i) = ttwinnor1(i)/xtwinnormag
          end do
        end do

      end if
      
      return
      end 

