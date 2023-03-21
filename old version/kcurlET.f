C     STRAIN GRADIENTS
C *******************************************************************************
C *                     Compute Curl of a 2nd order tensor                      *
C *              As a vector, curl of row written in corresponding column       *
C *                          Full integration                                   *
C *******************************************************************************
       SUBROUTINE kcurlET(curlFp,Fp,xnat,gauss,gausscoords) 

       INCLUDE 'ABA_PARAM.INC'

      integer, parameter:: nnodes=8


      !scalars

       real*8,intent(out) :: curlFP(8,9)
      !arrays gauss is gauss coordinates of isoparametric element in (s1,s2,s3) space

      real*8,intent(in) :: xnat(nnodes,3),gauss(nnodes,3),
     + gausscoords(3,nnodes), Fp(8,9)

      real*8 :: J(3,3),Jinv(3,3),dNds(nnodes,3),dndx(nnodes,3),
     + fnode(3,nnodes), fmat1(3,3),dmout(3,3),N(nnodes), z(3)
     
      !Extrapolation arrays
      real*8 :: z11i(nnodes),z11n(nnodes),z12i(nnodes),z12n(nnodes),
     + z13i(nnodes),z13n(nnodes),z21i(nnodes),z21n(nnodes),z22i(nnodes),
     + z22n(nnodes),z23i(nnodes),z23n(nnodes),z31i(nnodes),z31n(nnodes),
     + z32i(nnodes),z32n(nnodes),z33i(nnodes),z33n(nnodes),
     + Nmat(nnodes,nnodes),xnmatI(nnodes,nnodes)

      real(kind=8):: ri,si,ti,gr,gs,gt,dgr,dgs,dgt,xgauss 
      
      fnode = 0.;fmat1 = 0.;Nmat=0.;xnmatI=0.
      z11i=0.;z11n=0.;z12i=0.;z12n=0.;z13i=0.;z13n=0.
      z21i=0.;z21n=0.;z22i=0.;z22n=0.;z23i=0.;z23n=0.
      z31i=0.;z31n=0.;z32i=0.;z32n=0.;z33i=0.;z33n=0.
C
C   USE EIGHT GAUSS POINT FOR GRADIENTS
C
C    LOOP OVER EIGHT INTEGRATION POINTS

!    Evaluate curlT at the integration points, and simultaneously populate Nmat which will be later inverted for extrapolation purposes.   

      do kint2 = 1,8

C    SPECIFY z - INTEGRATION POINT
    
      z = gauss(kint2,1:3)

C    SHAPE FUNCTIONS AND DERIVATIVES  
      do i = 1,nnodes 
       ri = xnat(i,1)
       si = xnat(i,2)
       ti = xnat(i,3)       
       !---------------------------
       gr = 0.5*(1. + ri*z(1))
       dgr = 0.5*ri
       !---------------------------       
       gs = 0.5*(1.+si*z(2))
       dgs = 0.5*si       
       !---------------------------       
       gt = 0.5*(1.+ti*z(3))
       dgt = 0.5*ti
       !---------------------------
       N(i) = gr*gs*gt
       !---------------------------
       dNds(i, 1) = dgr*gs*gt ! dnds1(i)
       dNds(i, 2) = gr*dgs*gt ! dnds2(i)
       dNds(i, 3) = gr*gs*dgt ! dnds3(i)
      end do       
      
      Nmat(kint2,:) = N    
C   
C     SET UP JACOBIAN           
C  
      J = matmul(gausscoords,dNds)

C    AND ITS INVERSE
C
      call KDETER(J,det)
      
      if (abs(det) <= 1.0e-6 .or. det /= det) then !last part true if det=NaN
         dmout = 0.0
      else  
         call lapinverse(J,3,info,Jinv)
         
         if(info /= 0) then 
             write(6,*) "inverse failure: J in kcurl"             
          end if
C
      dndx = matmul(dNds,Jinv) 
C
C    DETERMINE first column of curlf: Read row, to determine column.
C
C   
      fnode(1:3, 1:8) = transpose(Fp(1:8,1:3))
      fmat1 = matmul(fnode,dndx)
      
      
      !Curlfp at integeration points
      z11i(kint2) = fmat1(3,2) - fmat1(2,3)
      z21i(kint2) = fmat1(1,3) - fmat1(3,1)
      z31i(kint2) = fmat1(2,1) - fmat1(1,2)
            
C
C    DETERMINE second column of curlf
C
C
      fnode(1:3, 1:8) = transpose(Fp(1:8,4:6))
      fmat1 = matmul(fnode,dndx)
     
      
      !Curlfp at integeration points
      z12i(kint2) = fmat1(3,2) - fmat1(2,3)
      z22i(kint2) = fmat1(1,3) - fmat1(3,1)
      z32i(kint2) = fmat1(2,1) - fmat1(1,2)             
C
C    DETERMINE third column of curlf
C
C
      fnode(1:3, 1:8) = transpose(Fp(1:8,7:9))
      fmat1 = matmul(fnode,dndx)
    

      !Curlfp at integeration points
      z13i(kint2) = fmat1(3,2) - fmat1(2,3)
      z23i(kint2) = fmat1(1,3) - fmat1(3,1)
      z33i(kint2) = fmat1(2,1) - fmat1(1,2)         
C
      end if

      
      end do !kint2
      
      !All integration points done. Extrapolation begins
      
      call lapinverse(Nmat,nnodes,info2,xnmatI)
      if(info2 /= 0) then
          write(6,*) "inverse failure: xnmat in kcurl"
      end if      
      
      z11n = matmul(xnmatI,z11i)
      z21n = matmul(xnmatI,z21i)
      z31n = matmul(xnmatI,z31i)
      
      z12n = matmul(xnmatI,z12i)
      z22n = matmul(xnmatI,z22i)
      z32n = matmul(xnmatI,z32i)
      
      z13n = matmul(xnmatI,z13i)
      z23n = matmul(xnmatI,z23i)
      z33n = matmul(xnmatI,z33i)


      !The storage is done by row.
      do kint=1,nnodes
          curlFp(kint,1) = z11n(kint)
          curlFp(kint,2) = z12n(kint)
          curlFp(kint,3) = z13n(kint)
          
          curlFp(kint,4) = z21n(kint)
          curlFp(kint,5) = z22n(kint)
          curlFp(kint,6) = z23n(kint)
          
          curlFp(kint,7) = z31n(kint)
          curlFp(kint,8) = z32n(kint)
          curlFp(kint,9) = z33n(kint)
          
      end do !kint 
C     
      RETURN
      END
C
C
