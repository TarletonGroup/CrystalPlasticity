!     Jan. 01st, 2023
!     Eralp Demir      
!
!      
!     ABAQUS Analysis User's Manual (6.10)
!     2D-Plane strain elements
!     CPE4   	4-node bilinear         4-integration points
!     CPE6   	6-node quadratic        3-integration points
!     CPE8   	8-node biquadratic      9-integration points
!     CPE8R  	8-node biquadratic      4-integration points (reduced integration)
!
!
!     2D-Plane stress elements
!     CPS4   	4-node bilinear         4-integration points
!     CPS6   	6-node quadratic        3-integration points
!     CPS8   	8-node biquadratic      9-integration points
!     CPS8R  	8-node biquadratic      4-integration points (reduced integration)
!   
!
!     3D - element types
!     C3D6    6-node linear triangular prism      4-integration points
!     C3D8    8-node linear brick                 8-integration points
!     C3D10   10-node quadratic tetrahedron       4-integration points
!     C3D15   15-node quadratic triangular prism  9-integration points
!     C3D20   20-node quadratic brick             27-integration points
!     C3D20R  20-node quadratic brick             8-integration points (reduced integration)
!
      module meshprop
      implicit none
      contains
!
!     Finite element properties
!     based on element type
      subroutine feprop
      use errors, only : error
      use globalvariables, only : eltyp, numel, 
     + numtens, numdim, numpt, nnpel,
     + Nmat, invNmat, dNmat, dNmatc, 
     + ipweights, calculategradient
      use userinputs, only : gndlinear, gndmodel
      use utilities, only: nolapinverse
!
!
      implicit none
!
!     Gaussian point coordinates
      real(8) :: a, b
!     Parametric coordinates
      real(8) :: g, h, r
!
!     Separate variables are used for each element type
!     in order to avoid using allocatables!
!
!     Element type: CPE4 or CPS4 or CPE8R or CPS8R
!     GP coordinates (numpt x numdim)
      real(8) :: GPcoords_CP4(4,2)
!     Shape functions (nnpel)
      real(8) :: N_CP4(4)
!     Shape function derivatives (numdim x nnpel)
      real(8) :: dN_CP4(2,4)
!     Shape function matrix (numpt x nnpel)
      real(8) :: Nmat_CP4(4,4)
!     Inverse of the shape function matrix (nnpel x numpt)
      real(8) :: invNmat_CP4(4,4)
!     Shape function derivatives (numdim x nnpel x numpt)
      real(8) :: dNmat_CP4(4,2,4)
!     Shape function derivatives at the center (numdim x nnpel)    
      real(8) :: dNmatc_CP4(2,4)
!
!
!     Element type: CPE6 or CPS6
!     GP coordinates (numpt x numdim)
      real(8) :: GPcoords_CP6(3,2)
!     Shape functions (nnpel)
!     Linear version
      real(8) :: N_CP3(3)
!     Quadratic version
      real(8) :: N_CP6(6)
!     Shape function derivatives (numdim x nnpel)
!     Linear version
      real(8) :: dN_CP3(2,3)
!     Quadratic version
      real(8) :: dN_CP6(2,6)
!     Shape function matrix
!     Linear version (numpt x nnpel)
      real(8) :: Nmat_CP3(3,3)
!     Quadratic version (numpt+numpt_sub x nnpel)
      real(8) :: Nmat_CP6(6,6)
!     Inverse of the shape function matrix (nnpel x numpt)
!     Linear version
      real(8) :: invNmat_CP3(3,3)
!     Quadratic version
      real(8) :: invNmat_CP6(6,3)
!     Shape function derivatives (numdim x nnpel x numpt)
!     Linear version
      real(8) :: dNmat_CP3(3,2,3)
!     Quadratic version
      real(8) :: dNmat_CP6(3,2,6)
!     Shape function derivatives at the center (numdim x nnpel)
!     Linear version
      real(8) :: dNmatc_CP3(2,3)
!     Quadratic version      
      real(8) :: dNmatc_CP6(2,6)
!     Sub element properties
!     Global coordinates (numpt_sub x numdim)
      real(8) :: GPcoords_sub_CP6(3,2)
!     Local coordinates (numpt_sub x numdim)
      real(8) :: GPcoords_sub_CP3(3,2)
!     Sub element mapping for additional IPs (numpt+numpt_sub x numpt)
      real(8) :: IPmap_CP6(6,3)
!     Inverse of the square matrix (nnpel x numpt + numpt_sub)
      real(8) :: coeff_CP6(6,6)
!
!
!     Element type: CPE8 or CPS8
!     GP coordinates (numpt x numdim)
      real(8) :: GPcoords_CP8(9,2)
!     Shape functions (nnpel)
      real(8) :: N_CP8(8)
!     Shape function derivatives (numdim x nnpel)
      real(8) :: dN_CP8(2,8)
!     Shape function matrix (numpt x nnpel)
      real(8) :: Nmat_CP8(9,8)
!     inverted square matrix (nnpel x nnpel)
      real(8) :: coeff_CP8(8,8)
!     Inverse of the shape function matrix (nnpel x numpt)
      real(8) :: invNmat_CP8(8,9)
!     Shape function derivatives (numdim x nnpel x numpt)
      real(8) :: dNmat_CP8(9,2,8)
!     Shape function derivatives at the center (numdim x nnpel)
      real(8) :: dNmatc_CP8(2,8)
!     Square matrix (nnpel x nnpel)
      real(8) :: NTN_CP8(8,8)
!
!     Shape function matrix (numpt x nnpel)
      real(8) :: Nmat_CP8lin(9,4)
!     Inverse of the shape function matrix (nnpel x numpt)
      real(8) :: invNmat_CP8lin(4,9)
!     Shape function derivatives (numpt x numdim x nnpel)
      real(8) :: dNmat_CP8lin(9,2,4)
!     Shape function derivatives at the center (numdim x nnpel)
      real(8) :: dNmatc_CP8lin(2,4)
!     Square matrix (nnpel x nnpel)
      real(8) :: NTN_CP8lin(4,4)
!     inverted square matrix (nnpel x nnpel)
      real(8) :: coeff_CP8lin(4,4)
!
!
!     Element type: C3D8 or C3D20R
!     GP coordinates (numpt x numdim)
      real(8) :: GPcoords_C3D8(8,3)
!     Shape functions (nnpel)
      real(8) :: N_C3D8(8)
!     Shape function derivatives (numdim x nnpel)
      real(8) :: dN_C3D8(3,8)
!     Shape function matrix (numpt x nnpel)
      real(8) :: Nmat_C3D8(8,8)
!     Inverse of the shape function matrix (nnpel x numpt)
      real(8) :: invNmat_C3D8(8,8)
!     Shape function derivatives (numdim x nnpel x numpt)
      real(8) :: dNmat_C3D8(8,3,8)
!     Shape function derivatives at the center (numdim x nnpel)
      real(8) :: dNmatc_C3D8(3,8)
!
!
!     Element type: C3D10
!     GP coordinates (numpt x numdim)
      real(8) :: GPcoords_C3D10(4,3)
!     Shape functions (nnpel)
!     Linear version
      real(8) :: N_C3D4(4)
!     Quadratic version
      real(8) :: N_C3D10(10)
!     Shape function derivatives (numdim x nnpel)
!     Linear version
      real(8) :: dN_C3D4(3,4)
!     Quadratic version
      real(8) :: dN_C3D10(3,10)
!     Shape function matrix
!     Linear version (numpt x nnpel)
      real(8) :: Nmat_C3D4(4,4)
!     Linear version (numpt_sub x nnpel)
      real(8) :: Nmat_sub_C3D4(6,4)
!     Quadratic version (numpt+numpt_sub x nnpel)
      real(8) :: Nmat_C3D10(10,10)
!     Inverse of the shape function matrix (nnpel x numpt)
!     Linear version
      real(8) :: invNmat_C3D4(4,4)
!     Quadratic version
      real(8) :: invNmat_C3D10(10,4)
!     Shape function derivatives (numdim x nnpel x numpt)
!     Linear version
      real(8) :: dNmat_C3D4(4,3,4)
!     Quadratic version
      real(8) :: dNmat_C3D10(4,3,10)
!     Shape function derivatives at the center (numdim x nnpel)
!     Linear version
      real(8) :: dNmatc_C3D4(3,4)
!     Quadratic version
      real(8) :: dNmatc_C3D10(3,10)
!     Sub element properties
!     Global coordinates (numpt_sub x numdim)
      real(8) :: GPcoords_sub_C3D10(6,3)
!     Local coordinates (numpt_sub x numdim)
      real(8) :: GPcoords_sub_C3D4(6,3)
!     Sub element mapping for additional IPs (numpt+numpt_sub x numpt)
      real(8) :: IPmap_C3D10(10,4)
!     Inverse of the square matrix (nnpel x numpt + numpt_sub)
      real(8) :: coeff_C3D10(10,10)
!
!
!
!     Element type: C3D15
!     GP coordinates (numpt x numdim)
      real(8) :: GPcoords_C3D15(9,3)
!     Shape functions (nnpel)
!     Linear version
      real(8) :: N_C3D6(6)
!     Quadratic version
      real(8) :: N_C3D15(15)
!     Shape function derivatives (numdim x nnpel)
!     Linear version
      real(8) :: dN_C3D6(3,6)
!     Quadratic version
      real(8) :: dN_C3D15(3,15)
!     Shape function matrix
!     Linear version (numpt x nnpel_sub)
      real(8) :: Nmat_C3D6(9,6)
!     Quadratic version (numpt x nnpel_sub)
      real(8) :: Nmat_sub_C3D6(6,6)
!     Quadratic version (numpt+numpt_sub x nnpel)
      real(8) :: Nmat_C3D15(15,15)
!     Coefficient matrix
!     Linear version (nnpel_sub x numpt)
      real(8) :: coeff_C3D6(6,6)
!     Linear version (nnpel_sub x numpt)
      real(8) :: invNmat_C3D6(6,9)
!     Quadratic version (nnpel x numpt)
      real(8) :: invNmat_C3D15(15,9)
!     Shape function derivatives (numdim x nnpel x numpt)
!     Linear version
      real(8) :: dNmat_C3D6(9,3,6)
!     Quadratic version
      real(8) :: dNmat_C3D15(9,3,15)
!     Shape function derivatives at the center (numdim x nnpel)
!     Linear version
      real(8) :: dNmatc_C3D6(3,6)
!     Quadratic version
      real(8) :: dNmatc_C3D15(3,15)
!     Sub element properties
!     Global coordinates (numpt_sub x numdim)
      real(8) :: GPcoords_sub_C3D15(6,3)
!     Local coordinates (numpt_sub x numdim)
      real(8) :: GPcoords_sub_C3D6(6,3)
!     Sub element mapping for additional IPs (numpt+numpt_sub x numpt)
      real(8) :: IPmap_C3D15(15,9)
!     Inverse of the square matrix (nnpel x numpt + numpt_sub)
      real(8) :: coeff_C3D15(15,15)
!     Square matrix (nnpel x nnpel)
      real(8) :: NTN_C3D6(6,6)  
!
!
!
!
!     Element type: C3D20
!     GP coordinates (numpt x numdim)
      real(8) :: GPcoords_C3D20(27,3)
!     Shape functions (nnpel)
      real(8) :: N_C3D20(20)
!     Shape function derivatives (numdim x nnpel)
      real(8) :: dN_C3D20(3,20)
!     Shape function matrix (numpt x nnpel)
      real(8) :: Nmat_C3D20(27,20)
!     inverted square matrix (nnpel x nnpel)
      real(8) :: coeff_C3D20(20,20)
!     Inverse of the shape function matrix (nnpel x numpt)
      real(8) :: invNmat_C3D20(20,27)
!     Shape function derivatives (numdim x nnpel x numpt)
      real(8) :: dNmat_C3D20(27,3,20)
!     Shape function derivatives at the center (numdim x nnpel)
      real(8) :: dNmatc_C3D20(3,20)
!     Square matrix (nnpel x nnpel)
      real(8) :: NTN_C3D20(20,20)
!
!     Shape function matrix (numpt x nnpel)
      real(8) :: Nmat_C3D20lin(27,8)
!     Inverse of the shape function matrix (nnpel x numpt)
      real(8) :: invNmat_C3D20lin(8,27)
!     Shape function derivatives (numpt x numdim x nnpel)
      real(8) :: dNmat_C3D20lin(27,3,8)
!     Shape function derivatives at the center (numdim x nnpel)
      real(8) :: dNmatc_C3D20lin(3,8)
!     Square matrix (nnpel x nnpel)
      real(8) :: NTN_C3D20lin(8,8)
!     inverted square matrix (nnpel x nnpel)
      real(8) :: coeff_C3D20lin(8,8)
!
!
!     Sub element properties
      integer :: nnpel_sub, numpt_sub
!     Parametric Gauss point coordinates
      integer :: ipt
!     Dummy integers
      integer :: i, j
!
!
!     Identify the element type
!     Based on 
!     - number of integration points per element: numpt
!     - dimension of the problem: ntens
      call findelementtype(numtens,numpt,eltyp)
!
!
!
!
!!     Output the element type
!      write(*,*) 'ABAQUS element type: ', eltyp
!
      select case(eltyp)
!
!     Element type: CPE3 or CPS3 or CPE6R or CPS6R
!     Use the same interpolation functions for the reduced elements
      case('CPE3', 'CPS3', 'CPE6R', 'CPS6R')
!
!         GND calculation is not possible for this element type
          if (gndmodel>0) then
              call error(7)
          end if
!
!         Number of dimensions of the problem
          numdim = 2
!
!
!         Single integration point
          calculategradient = 0
!
!         Number of nodes per element
          nnpel = 3
!
!         IP integration weights
          allocate(ipweights(numpt))
          ipweights=0.
!
!         Number of nodes per the sub element
          nnpel_sub = 0
!
!         Additional Gauss points for quadratic interpolation
          numpt_sub = 0
!
!     Element type: CPE4R or CPS4R
!     Use the same interpolation functions for the reduced elements
      case('CPE4R', 'CPS4R')
!
!         GND calculation is not possible for this element type
          if (gndmodel>0) then
              call error(7)
          end if
!
!         Number of dimensions of the problem
          numdim = 2
!
!
!         Single integration point
          calculategradient = 0
!
!         Number of nodes per element
          nnpel = 4
!
!         IP integration weights
          allocate(ipweights(numpt))
          ipweights=0.
!
!         Number of nodes per the sub element
          nnpel_sub = 0
!
!         Additional Gauss points for quadratic interpolation
          numpt_sub = 0
!
!     Element type: CPE4 or CPS4 or CPE8R or CPS8R
!     Use the same interpolation functions for the reduced elements
      case('CPE4', 'CPS4', 'CPE8R', 'CPS8R')
!
!         Number of dimensions of the problem
          numdim = 2
!
!
!
!         Number of nodes per element
          nnpel = 4
!
!         IP integration weights
          allocate(ipweights(numpt))
          ipweights=1.
!
!
!         Number of nodes per the sub element
          nnpel_sub = 0
!
!         Additional Gauss points for quadratic interpolation
          numpt_sub = 0    
!
!         Gauss points
          GPcoords_CP4 = 0.
          GPcoords_CP4(1,:) = (/ -1., -1. /)
          GPcoords_CP4(2,:) = (/  1., -1. /)
          GPcoords_CP4(3,:) = (/ -1.,  1. /)
          GPcoords_CP4(4,:) = (/  1.,  1. /)
          GPcoords_CP4 = GPcoords_CP4/sqrt(3.)
!
          write(*,*) 'Linear interpolation for will be used GNDs'
!
!         Compute shape functions and their derivatives at the integration points
          do ipt = 1, numpt
!
!
!             Gauss point coordinates
              g = GPcoords_CP4(ipt,1)
              h = GPcoords_CP4(ipt,2)
!
!             Calculate shape functions and their derivatives
              call CP4_N_dN(nnpel,numdim,g,h,N_CP4,dN_CP4)
!
!             Shape function mapping
              Nmat_CP4(ipt,1:nnpel) = N_CP4
!
!             Shape function derivative
              dNmat_CP4(ipt,1:numdim,1:nnpel) = dN_CP4
!
          end do          
!
!
!         Compute inverse of the interpolation matrix
          call nolapinverse(Nmat_CP4,invNmat_CP4,4)
!
!
!         Shape function derivatives at the element center
          g = 0.
          h = 0.
!
!         Calculate shape functions and their derivatives
          call CP4_N_dN(nnpel,numdim,g,h,N_CP4,dN_CP4)
!
!         Assign the center value
          dNmatc_CP4 = dN_CP4
!
!         Allocate global arrays
          allocate(Nmat(numpt,nnpel))
          Nmat=0.
          allocate(dNmat(numpt,numdim,nnpel))
          dNmat=0.
          allocate(invNmat(nnpel,numpt))
          invNmat=0.
          allocate(dNmatc(numdim,nnpel))
          dNmatc=0.
!
!         Assign global arrays
          Nmat(:,:) = Nmat_CP4(:,:)
          dNmat(:,:,:) = dNmat_CP4(:,:,:)
          invNmat(:,:) = invNmat_CP4(:,:)
          dNmatc(:,:) = dNmatc_CP4(:,:)
!
!
!     Element type: CPE6 or CPS6
      case('CPE6', 'CPS6')
!
!         Number of dimensions of the problem
          numdim = 2
!
!
!
!
!         IP integration weights
          allocate(ipweights(numpt))
          ipweights = 1./6.
!
!
!
!
!         Gauss points
          GPcoords_CP6 = 0.
          GPcoords_CP6(1,:) = (/ 1./6., 1./6. /)
          GPcoords_CP6(2,:) = (/ 2./3., 1./6. /)
          GPcoords_CP6(3,:) = (/ 1./6., 2./3. /)
!
!
!
!
!         Based on the selected option
!         Linear interpolation functions
          if (gndlinear==1) then
!
              write(*,*) 'Linear interpolation for will be used GNDs'
!
!             Number of nodes per element (linear element)
!             Assumed as a linear element
              nnpel = 3
!
!             Number of nodes per the sub element
              nnpel_sub = 0
!
!             Additional Gauss points for quadratic interpolation
              numpt_sub = 0
!
!
!             Use CP3 (linear) element function
!
!             Compute shape functions and their derivatives at the integration points          
              do ipt = 1, numpt
!
!                 Gauss point coordinates
                  g = GPcoords_CP6(ipt,1)
                  h = GPcoords_CP6(ipt,2)
!
!                 Calculate shape functions and their derivatives
                  call CP3_N_dN(nnpel,numdim,g,h,N_CP3,dN_CP3)
!
!                 Shape function mapping
                  Nmat_CP3(ipt,1:nnpel) = N_CP3
!
!                 Shape function derivative
                  dNmat_CP3(ipt,1:numdim,1:nnpel) = dN_CP3
!
              end do          
!
!
!             Compute inverse of the interpolation matrix
              call nolapinverse(Nmat_CP3,invNmat_CP3,3)
!
!
!             Shape function derivatives at the element center
              g = 1./3.
              h = 1./3.
!
!             Calculate shape functions and their derivatives
              call CP3_N_dN(nnpel,numdim,g,h,N_CP3,dN_CP3)      
!
!             Assign the center value
              dNmatc_CP3 = dN_CP3
!
!
!
!             Allocate global arrays
              allocate(Nmat(numpt,nnpel))
              Nmat=0.
              allocate(dNmat(numpt,numdim,nnpel))
              dNmat=0.
              allocate(invNmat(nnpel,numpt))
              invNmat=0.
              allocate(dNmatc(numdim,nnpel))
              dNmatc=0.
!
!             Assign global arrays
              Nmat(:,:) = Nmat_CP3(:,:)
              dNmat(:,:,:) = dNmat_CP3(:,:,:)
              invNmat(:,:) = invNmat_CP3(:,:)
              dNmatc(:,:) = dNmatc_CP3(:,:)
!
!
!
!
!         Quadratic interpolation functions
          else
!
              write(*,*) 'Quadratic interpolation will be used for GNDs'
!
!             Number of nodes per element
              nnpel = 6
!
!             Number of nodes per the sub element
              nnpel_sub = 3
!
!             Additional Gauss points for quadratic interpolation
              numpt_sub = nnpel-numpt
!
!             Calculated using the linear sub-element
!             Global coordinates (in its quadratic form)
              GPcoords_sub_CP6 = 0.
              GPcoords_sub_CP6(1,:) = (/ 5./12., 1./6. /)
              GPcoords_sub_CP6(2,:) = (/ 5./12., 5./12. /)
              GPcoords_sub_CP6(3,:) = (/ 1./6.,  5./12. /)
!
!             Local coordinates (in its linear form)
              GPcoords_sub_CP3 = 0.
              GPcoords_sub_CP3(1,:) = (/ 1./2., 0. /)
              GPcoords_sub_CP3(2,:) = (/ 1./2., 1./2. /)
              GPcoords_sub_CP3(3,:) = (/ 0., 1./2. /)
!
!
!
!             Use CP6 (quadratic) element function
!
!             Compute shape functions and their derivatives at the integration points
              do ipt = 1, numpt
!
!                 Gauss point coordinates
                  g = GPcoords_CP6(ipt,1)
                  h = GPcoords_CP6(ipt,2)
!
!                 Calculate shape functions and their derivatives
                  call CP6_N_dN(nnpel,numdim,g,h,N_CP6,dN_CP6)
!
!                 Shape function mapping
                  Nmat_CP6(ipt,1:nnpel) = N_CP6
!
!                 Shape function derivative
                  dNmat_CP6(ipt,1:numdim,1:nnpel) = dN_CP6
!
              end do
!
!
!
!
!
!
!             Compute shape functions and their derivatives at the integration points              
              do ipt = 1, numpt_sub
!
!                 Gauss point coordinates - global
                  g = GPcoords_sub_CP6(ipt,1)
                  h = GPcoords_sub_CP6(ipt,2)
!
!                 Calculate shape functions and their derivatives
                  call CP6_N_dN(nnpel,numdim,g,h,N_CP6,dN_CP6)
!
!                 Shape function mapping
                  Nmat_CP6(numpt+ipt,1:nnpel) = N_CP6
!
!
!
!                 Gauss point coordinates - local
                  g = GPcoords_sub_CP3(ipt,1)
                  h = GPcoords_sub_CP3(ipt,2)
!
!                 Calculate shape functions and their derivatives
                  call CP3_N_dN(nnpel_sub,numdim,g,h,N_CP3,dN_CP3)
!
!                 Shape function mapping
                  Nmat_CP3(ipt,1:nnpel_sub) = N_CP3
!
!
!
              end do
!
!
              IPmap_CP6=0.
!             Identity matrix
              do i=1,numpt
                  IPmap_CP6(i,i)=1.
              end do
              IPmap_CP6(numpt+1:numpt+numpt_sub,1:nnpel_sub) =
     + Nmat_CP3
!
!
!
!
!             Compute inverse of the interpolation matrix
              call nolapinverse(Nmat_CP6,coeff_CP6,6)
!
!
              invNmat_CP6 = matmul(coeff_CP6,IPmap_CP6)
!
!
!             Shape function derivatives at the element center
              g = 1./3.
              h = 1./3.
!
!             Calculate shape functions and their derivatives
              call CP6_N_dN(nnpel,numdim,g,h,N_CP6,dN_CP6)
!
!             Assign the center value
              dNmatc_CP6 = dN_CP6
!
!
!             Allocate global arrays
              allocate(Nmat(numpt+numpt_sub,nnpel))
              Nmat=0.
              allocate(dNmat(numpt,numdim,nnpel))
              dNmat=0.
              allocate(invNmat(nnpel,numpt))
              invNmat=0.
              allocate(dNmatc(numdim,nnpel))
              dNmatc=0.
!
!             Assign global arrays
              Nmat(:,:) = Nmat_CP6(:,:)
              dNmat(:,:,:) = dNmat_CP6(:,:,:)
              invNmat(:,:) = invNmat_CP6(:,:)
              dNmatc(:,:) = dNmatc_CP6(:,:)
!
!
!
!
!
          endif
!
!
!
!
!
!     Element type: CPE8 or CPS8
      case('CPE8', 'CPS8')
!
!         Number of dimensions of the problem
          numdim = 2
!
!
!
!         Integration weights
          allocate(ipweights(numpt))
          ipweights=0.
!
          ipweights(1) = 0.308641975308642
          ipweights(2) = 0.493827160493827
          ipweights(3) = 0.308641975308642
          ipweights(4) = 0.493827160493827
          ipweights(5) = 0.790123456790124
          ipweights(6) = 0.493827160493827
          ipweights(7) = 0.308641975308642
          ipweights(8) = 0.493827160493827
          ipweights(9) = 0.308641975308642
!
!         Number of nodes per the sub element
          nnpel_sub = 0
!
!         Additional Gauss points for quadratic interpolation
          numpt_sub = 0
!
!         Gauss points
          GPcoords_CP8 = 0.
          GPcoords_CP8(1,:) = (/ -1., -1. /)
          GPcoords_CP8(2,:) = (/  0., -1. /)
          GPcoords_CP8(3,:) = (/  1., -1. /)
          GPcoords_CP8(4,:) = (/ -1.,  0. /)
          GPcoords_CP8(5,:) = (/  0.,  0. /)
          GPcoords_CP8(6,:) = (/  1.,  0. /)
          GPcoords_CP8(7,:) = (/ -1.,  1. /)
          GPcoords_CP8(8,:) = (/  0.,  1. /)
          GPcoords_CP8(9,:) = (/  1.,  1. /)
          GPcoords_CP8 = sqrt(0.6) * GPcoords_CP8
!
!         Based on the selected option
!         Linear interpolation functions
          if (gndlinear==1) then
!
              write(*,*) 'Linear interpolation for will be used GNDs'
!
!             Number of nodes per element (linear element)
!             Assumed as a linear element
              nnpel = 4
!
!
!             Use CP4 (linear) element function
!             
!             Compute shape functions and their derivatives at the integration points              
              do ipt = 1, numpt
!
!                 Gauss point coordinates
                  g = GPcoords_CP8(ipt,1)
                  h = GPcoords_CP8(ipt,2)
!
!                 Calculate shape functions and their derivatives
                  call CP4_N_dN(nnpel,numdim,g,h,N_CP4,dN_CP4)
!
!                 Shape function mapping
                  Nmat_CP8lin(ipt,1:nnpel) = N_CP4
!
!                 Shape function derivative
                  dNmat_CP8lin(ipt,1:numdim,1:nnpel) = dN_CP4
!
              end do
!
!
!             Make Nmat a square matrix
              NTN_CP8lin = matmul(transpose(Nmat_CP8lin),Nmat_CP8lin)
!
!
!             Compute inverse of the interpolation matrix
              call nolapinverse(NTN_CP8lin,coeff_CP8lin,4)
!
!             Overall mapping for mapping from IPs to nodes
              invNmat_CP8lin=matmul(coeff_CP8lin,transpose(Nmat_CP8lin))
!
!
!             Shape function derivatives at the element center
              g = 0.
              h = 0.
!
!             Calculate shape functions and their derivatives
              call CP4_N_dN(nnpel,numdim,g,h,N_CP4,dN_CP4)
!
!             Assign the center value
              dNmatc_CP8lin = dN_CP4
!
!
!
!             Allocate global arrays
              allocate(Nmat(numpt,nnpel))
              Nmat=0.
              allocate(dNmat(numpt,numdim,nnpel))
              dNmat=0.
              allocate(invNmat(nnpel,numpt))
              invNmat=0.
              allocate(dNmatc(numdim,nnpel))
              dNmatc=0.
!
!             Assign global arrays
              Nmat(:,:) = Nmat_CP8lin(:,:)
              dNmat(:,:,:) = dNmat_CP8lin(:,:,:)
              invNmat(:,:) = invNmat_CP8lin(:,:)
              dNmatc(:,:) = dNmatc_CP8lin(:,:)
!
!
!
!
          else
!
!
!
!
!             Number of nodes per element
              nnpel = 8
!
!
              write(*,*) 'Quadratic interpolation will be used for GNDs'
!
!
!             Compute shape functions and their derivatives at the integration points
              do ipt = 1, numpt
!
!                 Gauss point coordinates
                  g = GPcoords_CP8(ipt,1)
                  h = GPcoords_CP8(ipt,2)
!
!                 Calculate shape functions and their derivatives
                  call CP8_N_dN(nnpel,numdim,g,h,N_CP8,dN_CP8)
!
!                 Shape function mapping
                  Nmat_CP8(ipt,1:nnpel) = N_CP8
!
!                 Shape function derivative
                  dNmat_CP8(ipt,1:numdim,1:nnpel) = dN_CP8
!
              end do
!
!             Make Nmat a square matrix
              NTN_CP8 = matmul(transpose(Nmat_CP8),Nmat_CP8)
!
!
!             Compute inverse of the interpolation matrix
              call nolapinverse(NTN_CP8,coeff_CP8,8)
!
!             Overall mapping for mapping from IPs to nodes
              invNmat_CP8 = matmul(coeff_CP8,transpose(Nmat_CP8))
!
!
!             Shape function derivatives at the element center
              g = 0.
              h = 0.
!
!             Calculate shape functions and their derivatives
              call CP8_N_dN(nnpel,numdim,g,h,N_CP8,dN_CP8)
!
!             Assign the center value
              dNmatc_CP8 = dN_CP8
!
!
!             Allocate global arrays
              allocate(Nmat(numpt,nnpel))
              Nmat=0.
              allocate(dNmat(numpt,numdim,nnpel))
              dNmat=0.
              allocate(invNmat(nnpel,numpt))
              invNmat=0.
              allocate(dNmatc(numdim,nnpel))
              dNmatc=0.
!
!             Assign global arrays
              Nmat(:,:) = Nmat_CP8(:,:)
              dNmat(:,:,:) = dNmat_CP8(:,:,:)
              invNmat(:,:) = invNmat_CP8(:,:)
              dNmatc(:,:) = dNmatc_CP8(:,:)
!
!
          end if
!
!
      case('C3D4')
!
!         GND calculation is not possible for this element type
          if (gndmodel>0) then
              call error(7)
          end if
!
!
!
!         Number of dimensions of the problem
          numdim = 3
!
!
!         Single integration point
          calculategradient = 0
!
!         Number of nodes per element      
          nnpel = 4
!
!         IP integration weights
          allocate(ipweights(numpt))
          ipweights=0.
!
!         Number of nodes per the sub element
          nnpel_sub = 0
!
!         Additional Gauss points for quadratic interpolation
          numpt_sub = 0
!
!
      case('C3D6')
!
!         GND calculation is not possible for this element type
          if (gndmodel>0) then
              call error(7)
          end if
!
!
!
!         Number of dimensions of the problem
          numdim = 3
!
!
!         Single integration point
          calculategradient = 0
!
!         Number of nodes per element
          nnpel = 6
!
!         IP integration weights
          allocate(ipweights(numpt))
          ipweights=0.
!
!         Number of nodes per the sub element
          nnpel_sub = 0
!
!         Additional Gauss points for quadratic interpolation
          numpt_sub = 0
!
!
!
!
      case('C3D8R')
!
!         GND calculation is not possible for this element type
          if (gndmodel>0) then
              call error(7)
          end if
!
!
!
!         Number of dimensions of the problem
          numdim = 3
!
!
!         Single integration point
          calculategradient = 0
!
!         Number of nodes per element      
          nnpel = 8
!
!         IP integration weights
          allocate(ipweights(numpt))
          ipweights=0.
!
!         Number of nodes per the sub element
          nnpel_sub = 0
!
!         Additional Gauss points for quadratic interpolation
          numpt_sub = 0
!
!
!
!
!     Element type: C3D8 or C3D20R
!     Use the same interpolation functions for the reduced element
      case('C3D8','C3D20R')
!
!         Number of dimensions of the problem
          numdim = 3
!
!
!
!         Number of nodes per element  
          nnpel = 8
!
!         IP integration weights
          allocate(ipweights(numpt))
          ipweights = 1.
!
!
!         Number of nodes per the sub element
          nnpel_sub = 0
!
!         Additional Gauss points for quadratic interpolation
          numpt_sub = 0
!         
!
!         Gauss points
          GPcoords_C3D8 = 0.
          GPcoords_C3D8(1,:) = (/ -1., -1., -1. /)
          GPcoords_C3D8(2,:) = (/  1., -1., -1. /)
          GPcoords_C3D8(3,:) = (/ -1.,  1., -1. /)
          GPcoords_C3D8(4,:) = (/  1.,  1., -1. /)
          GPcoords_C3D8(5,:) = (/ -1., -1.,  1. /)
          GPcoords_C3D8(6,:) = (/  1., -1.,  1. /)
          GPcoords_C3D8(7,:) = (/ -1.,  1.,  1. /)
          GPcoords_C3D8(8,:) = (/  1.,  1.,  1. /)
          GPcoords_C3D8 = GPcoords_C3D8/sqrt(3.)
!
!
          write(*,*) 'Linear interpolation for will be used GNDs'
!
!
!
!
!         Compute shape functions and their derivatives at the integration points
          do ipt = 1, numpt
!
!             Gauss point coordinates
              g = GPcoords_C3D8(ipt,1)
              h = GPcoords_C3D8(ipt,2)
              r = GPcoords_C3D8(ipt,3)
!
!             
!             Calculate shape functions and their derivatives
              call C3D8_N_dN(nnpel,numdim,g,h,r,N_C3D8,dN_C3D8)
!
!
!
!             Shape function mapping
              Nmat_C3D8(ipt,1:nnpel) = N_C3D8
!
!
!
!             Shape function derivative
              dNmat_C3D8(ipt,1:numdim,1:nnpel) = dN_C3D8
!
!
          end do
!
!
!
!
!         Compute inverse of the interpolation matrix
          call nolapinverse(Nmat_C3D8,invNmat_C3D8,8)
!
!
!         Shape function derivatives at the element center
          g = 0.
          h = 0.
          r = 0.
!
!         Calculate shape functions and their derivatives
          call C3D8_N_dN(nnpel,numdim,g,h,r,N_C3D8,dN_C3D8)
!
!         Assign the center value
          dNmatc_C3D8 = dN_C3D8
!
!
!
!
!
!         Allocate global arrays
          allocate(Nmat(numpt,nnpel))
          Nmat=0.
          allocate(dNmat(numpt,numdim,nnpel))
          dNmat=0.
          allocate(invNmat(nnpel,numpt))
          invNmat=0.
          allocate(dNmatc(numdim,nnpel))
          dNmatc=0.
!
!
!
!
!
!         Assign global arrays
          Nmat(:,:) = Nmat_C3D8(:,:)
          dNmat(:,:,:) = dNmat_C3D8(:,:,:)
          invNmat(:,:) = invNmat_C3D8(:,:)
          dNmatc(:,:) = dNmatc_C3D8(:,:)
!
!
!
!
!
!     Element type: C3D10
      case('C3D10')
!
!         Number of dimensions of the problem
          numdim = 3
!
!
!
!
!         IP integration weights
          allocate(ipweights(numpt))
          ipweights = 1./4./6.
!
!         Gauss points
          a = 0.138196601125010
          b = 0.585410196624968
!
!
          GPcoords_C3D10 = 0.
          GPcoords_C3D10(1,:) = (/ a, a, a /) 
          GPcoords_C3D10(2,:) = (/ b, a, a /)
          GPcoords_C3D10(3,:) = (/ a, b, a /)
          GPcoords_C3D10(4,:) = (/ a, a, b /)
!
!
!
!         Based on the selected option
!         Linear interpolation functions
          if (gndlinear==1) then
!
!             Number of nodes per element (linear element)
!             Assumed as a linear element
              nnpel = 4
!
!
!             Number of nodes per the sub element
              nnpel_sub = 0
              
!             Additional Gauss points for quadratic interpolation
              numpt_sub = 0
!
!
!             Use C3D4 (linear) element function
              write(*,*) 'Linear interpolation for will be used GNDs'
!
!             Compute shape functions and their derivatives at the integration points      
              do ipt = 1, numpt
!
!                 Gauss point coordinates
                  g = GPcoords_C3D10(ipt,1)
                  h = GPcoords_C3D10(ipt,2)
                  r = GPcoords_C3D10(ipt,3)
!
!                 Calculate shape functions and their derivatives
                  call C3D4_N_dN(nnpel,numdim,g,h,r,N_C3D4,dN_C3D4)
!
!                 Shape function mapping
                  Nmat_C3D4(ipt,1:nnpel) = N_C3D4
!
!                 Shape function derivative
                  dNmat_C3D4(ipt,1:numdim,1:nnpel) = dN_C3D4
!
              end do          
!
!
!             Compute inverse of the interpolation matrix
              call nolapinverse(Nmat_C3D4,invNmat_C3D4,4)
!
!
!             Shape function derivatives at the element center
              g = 1./4.
              h = 1./4.
              r = 1./4.
!
!             Calculate shape functions and their derivatives
              call C3D4_N_dN(nnpel,numdim,g,h,r,N_C3D4,dN_C3D4)    
!
!             Assign the center value
              dNmatc_C3D4 = dN_C3D4
!
!
!
!             Allocate global arrays
              allocate(Nmat(numpt,nnpel))
              Nmat=0.
              allocate(dNmat(numpt,numdim,nnpel))
              dNmat=0.
              allocate(invNmat(nnpel,numpt))
              invNmat=0.
              allocate(dNmatc(numdim,nnpel))
              dNmatc=0.
!
!             Assign global arrays
              Nmat(:,:) = Nmat_C3D4(:,:)
              dNmat(:,:,:) = dNmat_C3D4(:,:,:)
              invNmat(:,:) = invNmat_C3D4(:,:)
              dNmatc(:,:) = dNmatc_C3D4(:,:)
!
!
!
!
!         Quadratic interpolation functions
          else
!
!             Number of nodes per element
              nnpel = 10
!
!             Number of nodes per the sub element
              nnpel_sub = 4
!
!             Additional Gauss points for quadratic interpolation
              numpt_sub = nnpel-numpt          
!
!             Calculated using the linear sub-element
!             Global coordinates (in its quadratic form)
              GPcoords_sub_C3D10 = 0.
              do i=1,numdim
                  GPcoords_sub_C3D10(1,i) =
     + 0.5*(GPcoords_C3D10(1,i)+GPcoords_C3D10(2,i))
                  GPcoords_sub_C3D10(2,i) =
     + 0.5*(GPcoords_C3D10(2,i)+GPcoords_C3D10(3,i))
                  GPcoords_sub_C3D10(3,i) =
     + 0.5*(GPcoords_C3D10(3,i)+GPcoords_C3D10(1,i))
                  GPcoords_sub_C3D10(4,i) =
     + 0.5*(GPcoords_C3D10(4,i)+GPcoords_C3D10(1,i))
                  GPcoords_sub_C3D10(5,i) =
     + 0.5*(GPcoords_C3D10(2,i)+GPcoords_C3D10(4,i))
                  GPcoords_sub_C3D10(6,i) =
     + 0.5*(GPcoords_C3D10(3,i)+GPcoords_C3D10(4,i))
              end do
!
!
!             Local coordinates (in its linear form)
              GPcoords_sub_C3D4 = 0.
              GPcoords_sub_C3D4(1,:) = (/ 0.5, 0.,  0. /)
              GPcoords_sub_C3D4(2,:) = (/ 0.5, 0.5, 0. /)
              GPcoords_sub_C3D4(3,:) = (/  0., 0.5, 0. /)
              GPcoords_sub_C3D4(4,:) = (/  0., 0., 0.5 /)
              GPcoords_sub_C3D4(5,:) = (/ 0.5, 0., 0.5 /)
              GPcoords_sub_C3D4(6,:) = (/ 0., 0.5, 0.5 /)
!
!
!
!
!
!             Use C3D10 (quadratic) element function
              write(*,*) 'Quadratic interpolation for will be used GNDs'
!
!             Compute shape functions and their derivatives at the integration points
              do ipt = 1, numpt
!
!                 Gauss point coordinates
                  g = GPcoords_C3D10(ipt,1)
                  h = GPcoords_C3D10(ipt,2)
                  r = GPcoords_C3D10(ipt,3)
!
!                 Calculate shape functions and their derivatives
                  call C3D10_N_dN(nnpel,numdim,g,h,r,N_C3D10,dN_C3D10)
!
!                 Shape function mapping
                  Nmat_C3D10(ipt,1:nnpel) = N_C3D10
!
!                 Shape function derivative
                  dNmat_C3D10(ipt,1:numdim,1:nnpel) = dN_C3D10
!
!
              end do                 
!
!
!
!             Compute shape functions and their derivatives at the integration points
              do ipt = 1, numpt_sub
!
!                 Gauss point coordinates - global
                  g = GPcoords_sub_C3D10(ipt,1)
                  h = GPcoords_sub_C3D10(ipt,2)
                  r = GPcoords_sub_C3D10(ipt,3)
!
!                 Calculate shape functions and their derivatives
                  call C3D10_N_dN(nnpel,numdim,g,h,r,N_C3D10,dN_C3D10)
!
!                 Shape function mapping
                  Nmat_C3D10(numpt+ipt,1:nnpel) = N_C3D10
!
!                 
!                 Gauss point coordinates - local
                  g = GPcoords_sub_C3D4(ipt,1)
                  h = GPcoords_sub_C3D4(ipt,2)
                  r = GPcoords_sub_C3D4(ipt,3)
!
!                 Calculate shape functions and their derivatives
                  call C3D4_N_dN(nnpel_sub,numdim,g,h,r,N_C3D4,dN_C3D4)
!
!                 Shape function mapping
                  Nmat_sub_C3D4(ipt,1:nnpel_sub) = N_C3D4
!
!
              end do
!
!
!             Identity matrix
              IPmap_C3D10 = 0.
              do i=1,numpt
                  IPmap_C3D10(i,i)=1.
              end do
!
              IPmap_C3D10(numpt+1:numpt+numpt_sub,1:nnpel_sub) =
     + Nmat_sub_C3D4
!
!
!
!
!
!             Compute inverse of the interpolation matrix
              call nolapinverse(Nmat_C3D10,coeff_C3D10,10)
!
!
!             
              invNmat_C3D10 = matmul(coeff_C3D10,IPmap_C3D10)
!
!
!             Shape function derivatives at the element center
              g = 1./4.
              h = 1./4.
              r = 1./4.
!
!             Calculate shape functions and their derivatives
              call C3D10_N_dN(nnpel,numdim,g,h,r,N_C3D10,dN_C3D10)
!
!             Assign the center value
              dNmatc_C3D10 = dN_C3D10
!
!
!
!             Allocate global arrays
              allocate(Nmat(numpt+numpt_sub,nnpel))
              Nmat=0.
              allocate(dNmat(numpt,numdim,nnpel))
              dNmat=0.
              allocate(invNmat(nnpel,numpt))
              invNmat=0.
              allocate(dNmatc(numdim,nnpel))
              dNmatc=0.
!
!
!
!
!             Assign global arrays
              Nmat(:,:) = Nmat_C3D10(:,:)
              dNmat(:,:,:) = dNmat_C3D10(:,:,:)
              invNmat(:,:) = invNmat_C3D10(:,:)
              dNmatc(:,:) = dNmatc_C3D10(:,:)
!
!
!
!
!
!
!
!
          endif
!
!
!
!     Element type: C3D15   
      case('C3D15')
!
!         Number of dimensions of the problem
          numdim = 3
!
!
!
!         IP integration weights
          allocate(ipweights(numpt))
          ipweights=1./6./3.
!
!
!
!         Isoparametric coordinate (r-axis)
          a = 0.774596669241483
!
!
          GPcoords_C3D15 = 0.
          GPcoords_C3D15(2,:) = (/ 2./3., 1./6., -1. /)
          GPcoords_C3D15(1,:) = (/ 1./6., 1./6., -1. /)
          GPcoords_C3D15(3,:) = (/ 1./6., 2./3., -1. /)
          GPcoords_C3D15(5,:) = (/ 2./3., 1./6.,  0. /)
          GPcoords_C3D15(4,:) = (/ 1./6., 1./6.,  0. /)
          GPcoords_C3D15(6,:) = (/ 1./6., 2./3.,  0. /)
          GPcoords_C3D15(8,:) = (/ 2./3., 1./6.,  1. /)
          GPcoords_C3D15(7,:) = (/ 1./6., 1./6.,  1. /)
          GPcoords_C3D15(9,:) = (/ 1./6., 2./3.,  1. /)
!
          GPcoords_C3D15(:,3) = GPcoords_C3D15(:,3) * a
!
!
!
!         Based on the selected option
!         Linear interpolation functions
          if (gndlinear==1) then
!
!             Number of nodes per element (linear element)
!             Assumed as a linear element
              nnpel = 6
!
!
!             Number of nodes per the sub element
              nnpel_sub = 0
!
!             Additional Gauss points for quadratic interpolation
              numpt_sub = 0
!
!
!
!             Use C3D6 (linear) element function
              write(*,*) 'Linear interpolation will be used for GNDs'
!
!             Compute shape functions and their derivatives at the integration points         
              do ipt = 1, numpt
!
!                 Gauss point coordinates
                  g = GPcoords_C3D15(ipt,1)
                  h = GPcoords_C3D15(ipt,2)
                  r = GPcoords_C3D15(ipt,3)
!
!                 Calculate shape functions and their derivatives
                  call C3D6_N_dN(nnpel,numdim,g,h,r,N_C3D6,dN_C3D6)
!
!                 Shape function mapping
                  Nmat_C3D6(ipt,1:nnpel) = N_C3D6
!
!                 Shape function derivative
                  dNmat_C3D6(ipt,1:numdim,1:nnpel) = dN_C3D6
!
              end do
!
!
              NTN_C3D6 = matmul(transpose(Nmat_C3D6),Nmat_C3D6)
!
!
!             Compute inverse of the interpolation matrix
              call nolapinverse(NTN_C3D6,coeff_C3D6,6)
!
!
              invNmat_C3D6 = matmul(coeff_C3D6,transpose(Nmat_C3D6))
!
!
!             Shape function derivatives at the element center
              g = 1./3.
              h = 1./3.
              r = 0.
!
!             Calculate shape functions and their derivatives
              call C3D6_N_dN(nnpel,numdim,g,h,r,N_C3D6,dN_C3D6)
!
!             Assign the center value
              dNmatc_C3D6 = dN_C3D6
!
!
!
!
!             Allocate global arrays
              allocate(Nmat(numpt,nnpel))
              Nmat=0.
              allocate(dNmat(numpt,numdim,nnpel))
              dNmat=0.
              allocate(invNmat(nnpel,numpt))
              invNmat=0.
              allocate(dNmatc(numdim,nnpel))
              dNmatc=0.
!
!
!             Assign global arrays
              Nmat(:,:) = Nmat_C3D6(:,:)
              dNmat(:,:,:) = dNmat_C3D6(:,:,:)
              invNmat(:,:) = invNmat_C3D6(:,:)
              dNmatc(:,:) = dNmatc_C3D6(:,:)
!
!
!
!
!
!
!
!         Quadratic interpolation functions
          else
!
!             Number of nodes per element      
              nnpel = 15
!
!             Number of nodes of the sub-element
              nnpel_sub = 6
!
!             Additional Gauss points for quadratic interpolation
              numpt_sub = nnpel-numpt
!             Calculated using the linear sub-element
!             Global coordinates (in its quadratic form)
              GPcoords_sub_C3D15 = 0.
              do i=1,numdim
                  GPcoords_sub_C3D15(1,i) =
     + 0.5*(GPcoords_C3D15(1,i)+GPcoords_C3D15(2,i))
                  GPcoords_sub_C3D15(2,i) =
     + 0.5*(GPcoords_C3D15(2,i)+GPcoords_C3D15(3,i))
                  GPcoords_sub_C3D15(3,i) =
     + 0.5*(GPcoords_C3D15(3,i)+GPcoords_C3D15(1,i))
                  GPcoords_sub_C3D15(4,i) =
     + 0.5*(GPcoords_C3D15(7,i)+GPcoords_C3D15(8,i))
                  GPcoords_sub_C3D15(5,i) =
     + 0.5*(GPcoords_C3D15(8,i)+GPcoords_C3D15(9,i))
                  GPcoords_sub_C3D15(6,i) =
     + 0.5*(GPcoords_C3D15(7,i)+GPcoords_C3D15(9,i))
              end do
!
!
!             Local coordinates (in its linear form)
              GPcoords_sub_C3D6 = 0.
              GPcoords_sub_C3D6(1,:) = (/ 0.5,  0., -1. /)
              GPcoords_sub_C3D6(2,:) = (/ 0.5, 0.5, -1. /)
              GPcoords_sub_C3D6(3,:) = (/  0., 0.5, -1. /)
              GPcoords_sub_C3D6(4,:) = (/ 0.5,  0.,  1. /)
              GPcoords_sub_C3D6(5,:) = (/ 0.5, 0.5,  1. /)
              GPcoords_sub_C3D6(6,:) = (/  0., 0.5,  1. /)
!
!
!
!
!
!             Use C3D15 (quadratic) element function
              write(*,*) 'Quadratic interpolation will be used for GNDs'
!
!             Compute shape functions and their derivatives at the integration points        
              do ipt = 1, numpt
!
!                 Gauss point coordinates
                  g = GPcoords_C3D15(ipt,1)
                  h = GPcoords_C3D15(ipt,2)
                  r = GPcoords_C3D15(ipt,3)
!
!                 Calculate shape functions and their derivatives
                  call C3D15_N_dN(nnpel,numdim,g,h,r,N_C3D15,dN_C3D15)
!
!                 Shape function mapping
                  Nmat_C3D15(ipt,1:nnpel) = N_C3D15
!
!                 Shape function derivative
                  dNmat_C3D15(ipt,1:numdim,1:nnpel) = dN_C3D15
!
!
              end do                 
!
!
!
!      
!             Compute shape functions and their derivatives at the integration points
              do ipt = 1, numpt_sub
!
!                 Gauss point coordinates - global
                  g = GPcoords_sub_C3D15(ipt,1)
                  h = GPcoords_sub_C3D15(ipt,2)
                  r = GPcoords_sub_C3D15(ipt,3)
!
!                 Calculate shape functions and their derivatives
                  call C3D15_N_dN(nnpel,numdim,g,h,r,N_C3D15,dN_C3D15)
!
!                 Shape function mapping
                  Nmat_C3D15(numpt+ipt,1:nnpel) = N_C3D15
!
!
!
!
!                 Gauss point coordinates - local
                  g = GPcoords_sub_C3D6(ipt,1)
                  h = GPcoords_sub_C3D6(ipt,2)
                  r = GPcoords_sub_C3D6(ipt,3)
!
!                 Calculate shape functions and their derivatives
                  call C3D6_N_dN(nnpel_sub,numdim,g,h,r,N_C3D6,dN_C3D6)
!
!                 Shape function mapping
                  Nmat_sub_C3D6(ipt,1:nnpel_sub) = N_C3D6
!
!
!
              end do
!
!
!             Identity matrix
              IPmap_C3D15=0.
              do i=1,numpt
                  IPmap_C3D15(i,i)=1.
              end do
              IPmap_C3D15(numpt+1:numpt+numpt_sub,1:nnpel_sub) =
     + Nmat_sub_C3D6
!
!
!
!
!             Compute inverse of the interpolation matrix
              call nolapinverse(Nmat_C3D15,coeff_C3D15,15)
!
!
              invNmat_C3D15 = matmul(coeff_C3D15,IPmap_C3D15)
!
!
!             Shape function derivatives at the element center
              g = 1./3.
              h = 1./3.
              r = 0.
!
!             Calculate shape functions and their derivatives
              call C3D15_N_dN(nnpel,numdim,g,h,r,N_C3D15,dN_C3D15)      
!
!             Assign the center value
              dNmatc_C3D15 = dN_C3D15
!
!
!
!             Allocate global arrays
              allocate(Nmat(numpt+numpt_sub,nnpel))
              Nmat=0.
              allocate(dNmat(numpt,numdim,nnpel))
              dNmat=0.
              allocate(invNmat(nnpel,numpt))
              invNmat=0.
              allocate(dNmatc(numdim,nnpel))
              dNmatc=0.
!
!
!             Assign global arrays
              Nmat(:,:) = Nmat_C3D15(:,:)
              dNmat(:,:,:) = dNmat_C3D15(:,:,:)
              invNmat(:,:) = invNmat_C3D15(:,:)
              dNmatc(:,:) = dNmatc_C3D15(:,:)
!
!
!
!
!
!
          endif
!
!
!
!
!
!
!     Element type: C3D20
      case('C3D20')
!
!         Number of dimensions of the problem
          numdim = 3
!
!
!
!         IP integration weights
          allocate(ipweights(numpt))
          ipweights=0.
!
          ipweights(1) = 0.171467764060357
          ipweights(2) = 0.274348422496571
          ipweights(3) = 0.171467764060357
          ipweights(4) = 0.274348422496571
          ipweights(5) = 0.438957475994513
          ipweights(6) = 0.274348422496571
          ipweights(7) = 0.171467764060357
          ipweights(8) = 0.274348422496571
          ipweights(9) = 0.171467764060357
          ipweights(10) = 0.274348422496571
          ipweights(11) = 0.438957475994513
          ipweights(12) = 0.274348422496571
          ipweights(13) = 0.438957475994513
          ipweights(14) = 0.702331961591221
          ipweights(15) = 0.438957475994513
          ipweights(16) = 0.274348422496571
          ipweights(17) = 0.438957475994513
          ipweights(18) = 0.274348422496571
          ipweights(19) = 0.171467764060357
          ipweights(20) = 0.274348422496571
          ipweights(21) = 0.171467764060357
          ipweights(22) = 0.274348422496571
          ipweights(23) = 0.438957475994513
          ipweights(24) = 0.274348422496571
          ipweights(25) = 0.171467764060357
          ipweights(26) = 0.274348422496571
          ipweights(27) = 0.171467764060357
!
!         Number of nodes per the sub element
          nnpel_sub = 0
!
!         Additional Gauss points for quadratic interpolation
          numpt_sub = 0
!
!
          GPcoords_C3D20 = 0.
          GPcoords_C3D20(1,:)=  (/ -1., -1., -1. /)
          GPcoords_C3D20(2,:)=  (/  0., -1., -1. /)
          GPcoords_C3D20(3,:)=  (/  1., -1., -1. /)
          GPcoords_C3D20(4,:)=  (/ -1.,  0., -1. /)
          GPcoords_C3D20(5,:)=  (/  0.,  0., -1. /)
          GPcoords_C3D20(6,:)=  (/  1.,  0., -1. /)
          GPcoords_C3D20(7,:)=  (/ -1.,  1., -1. /)
          GPcoords_C3D20(8,:)=  (/  0.,  1., -1. /)
          GPcoords_C3D20(9,:)=  (/  1.,  1., -1. /)
          GPcoords_C3D20(10,:)= (/ -1., -1.,  0. /)
          GPcoords_C3D20(11,:)= (/  0., -1.,  0. /)
          GPcoords_C3D20(12,:)= (/  1., -1.,  0. /)
          GPcoords_C3D20(13,:)= (/ -1.,  0.,  0. /)
          GPcoords_C3D20(14,:)= (/  0.,  0.,  0. /)
          GPcoords_C3D20(15,:)= (/  1.,  0.,  0. /)
          GPcoords_C3D20(16,:)= (/ -1.,  1.,  0. /)
          GPcoords_C3D20(17,:)= (/  0.,  1.,  0. /)
          GPcoords_C3D20(18,:)= (/  1.,  1.,  0. /)
          GPcoords_C3D20(19,:)= (/ -1., -1.,  1. /)
          GPcoords_C3D20(20,:)= (/  0., -1.,  1. /)
          GPcoords_C3D20(21,:)= (/  1., -1.,  1. /)
          GPcoords_C3D20(22,:)= (/ -1.,  0.,  1. /)
          GPcoords_C3D20(23,:)= (/  0.,  0.,  1. /)
          GPcoords_C3D20(24,:)= (/  1.,  0.,  1. /)
          GPcoords_C3D20(25,:)= (/ -1.,  1.,  1. /)
          GPcoords_C3D20(26,:)= (/  0.,  1.,  1. /)
          GPcoords_C3D20(27,:)= (/  1.,  1.,  1. /)
          GPcoords_C3D20 = GPcoords_C3D20 * sqrt(0.6)
!
!
!         Based on the selected option
!         Linear interpolation functions
          if (gndlinear==1) then
!
!
              write(*,*) 'Linear interpolation for will be used GNDs'
!
!             Number of nodes per element (linear element)
!             Assumed as a linear element
              nnpel = 8
!
!
!
!
!             Use C3D8 (linear) element function
!
!             Compute shape functions and their derivatives at the integration points
              do ipt = 1, numpt
!
!                 Gauss point coordinates
                  g = GPcoords_C3D20(ipt,1)
                  h = GPcoords_C3D20(ipt,2)
                  r = GPcoords_C3D20(ipt,3)
!
!                 Calculate shape functions and their derivatives
                  call C3D8_N_dN(nnpel,numdim,g,h,r,N_C3D8,dN_C3D8)
!
!                 Shape function mapping
                  Nmat_C3D20lin(ipt,1:nnpel) = N_C3D8
!
!                 Shape function derivative
                  dNmat_C3D20lin(ipt,1:numdim,1:nnpel) = dN_C3D8
!
              end do
!
!
!             Make Nmat a square matrix
              NTN_C3D20lin =
     + matmul(transpose(Nmat_C3D20lin),Nmat_C3D20lin)
!
!
!             Compute inverse of the interpolation matrix
              call nolapinverse(NTN_C3D20lin,coeff_C3D20lin,8)
!
!             Overall mapping for mapping from IPs to nodes
              invNmat_C3D20lin =
     + matmul(coeff_C3D20lin,transpose(Nmat_C3D20lin))
!
!
!             Shape function derivatives at the element center
              g = 0.
              h = 0.
              r = 0.
!
!             Calculate shape functions and their derivatives
              call C3D8_N_dN(nnpel,numdim,g,h,r,N_C3D8,dN_C3D8)
!
!             Assign the center value
              dNmatc_C3D20lin = dN_C3D8
!
!
!
!             Allocate global arrays
              allocate(Nmat(numpt,nnpel))
              Nmat=0.
              allocate(dNmat(numpt,numdim,nnpel))
              dNmat=0.
              allocate(invNmat(nnpel,numpt))
              invNmat=0.
              allocate(dNmatc(numdim,nnpel))
              dNmatc=0.
!
!             Assign global arrays
              Nmat(:,:) = Nmat_C3D20lin(:,:)
              dNmat(:,:,:) = dNmat_C3D20lin(:,:,:)
              invNmat(:,:) = invNmat_C3D20lin(:,:)
              dNmatc(:,:) = dNmatc_C3D20lin(:,:)
!
!
!
!
          else
!
!
!
!
!             Number of nodes per element
              nnpel = 20
!
              write(*,*) 'Quadratic interpolation will be used for GNDs'
!
!
!
!             Compute shape functions and their derivatives at the integration points
              do ipt = 1, numpt
!
!                 Gauss point coordinates
                  g = GPcoords_C3D20(ipt,1)
                  h = GPcoords_C3D20(ipt,2)
                  r = GPcoords_C3D20(ipt,3)
!
!                 Calculate shape functions and their derivatives
                  call C3D20_N_dN(nnpel,numdim,g,h,r,N_C3D20,dN_C3D20)
!
!                 Shape function mapping
                  Nmat_C3D20(ipt,1:nnpel) = N_C3D20
!
!                 Shape function derivative
                  dNmat_C3D20(ipt,1:numdim,1:nnpel) = dN_C3D20
!
              end do
!
!
!             Make Nmat a square matrix
              NTN_C3D20 = matmul(transpose(Nmat_C3D20),Nmat_C3D20)
!
!             Compute inverse of the interpolation matrix
              call nolapinverse(NTN_C3D20,coeff_C3D20,20)
!
!             Overall mapping to map from the IPs to the nodes
              invNmat_C3D20 = matmul(coeff_C3D20,transpose(Nmat_C3D20))
!
!             Shape function derivatives at the element center
              g = 0.
              h = 0.
              r = 0.
!
!             Calculate shape functions and their derivatives
              call C3D20_N_dN(nnpel,numdim,g,h,r,N_C3D20,dN_C3D20)
!
!             Assign the center value
              dNmatc_C3D20 = dN_C3D20
!
!
!             Allocate global arrays
              allocate(Nmat(numpt,nnpel))
              Nmat=0.
              allocate(dNmat(numpt,numdim,nnpel))
              dNmat=0.
              allocate(invNmat(nnpel,numpt))
              invNmat=0.
              allocate(dNmatc(numdim,nnpel))
              dNmatc=0.
!
!             Assign global arrays
              Nmat(:,:) = Nmat_C3D20(:,:)
              dNmat(:,:,:) = dNmat_C3D20(:,:,:)
              invNmat(:,:) = invNmat_C3D20(:,:)
              dNmatc(:,:) = dNmatc_C3D20(:,:)
!
!
          end if
!
!
!
      case default
!
          call error(2)
!
!
      end select
!
!
      write(*,*) 'Dimension of the analysis: ', numdim
      write(*,*) 'Number of integration points per element: ', numpt
!      write(*,*) 'Number of nodes per element: ', nnpel
!
!
!
      return
      end subroutine feprop
!
!
!
!     Find the number of nodes for displacement analysis (not for GND analysis)
      subroutine findelementtype(numtens,numpt,eltyp)
      implicit none
      integer, intent(in) :: numtens
      integer, intent(in) :: numpt
      character(len=:), allocatable, intent(out) :: eltyp
      
!
!
!     Determine the element type
!     Based on the dimension of the problem (numdim=ntens)
!
!     Dimension of the problem (2D-plane stress)
      if (numtens==3) then
!
          select case(numpt)
!
!         2D - CPS3 / CPS6R / CPS4R
          case(1)
!
              eltyp = 'CPS3'
!
!         2D - CPS4 / CPS8R
          case(4)
!
              eltyp = 'CPS4'
!
!         2D - CPS6
          case(3)
!
              eltyp = 'CPS6'
!
!         2D - 8 node quadratic quadrilateral
          case(9)
!
              eltyp = 'CPS8'
!
          end select    
!
!     Dimension of the problem (2D - Plane strain)
      elseif (numtens==4) then
!
          select case(numpt)
!
!         2D - CPS3 / CPS6R / CPS4R
          case(1)
!
              eltyp = 'CPE3'
!
!         2D - CPS4 / CPS8R
          case(4)
!
              eltyp = 'CPE4'
!
!         2D - CPS6
          case(3)
!
              eltyp = 'CPE6'
!
!         2D - 8 node quadratic quadrilateral
          case(9)
!
              eltyp = 'CPE8'
!
          end select
!
!
!     Dimension of the problem (3D)
      elseif (numtens==6) then
!
          select case(numpt)
!
!         3D - C3D4 / C3D8R / C3D10R
          case(1)
!
              eltyp = 'C3D4'
!
!         3D - C3D6 / C3D15R
          case(2)
!
              eltyp = 'C3D6' 
!       
!         3D - C3D8 / C3D20R
          case(8)
!
              eltyp = 'C3D8'
!
!         3D - C3D10
          case(4)
!
              eltyp = 'C3D10'
!
!         3D - C3D15
          case(9)
!
              eltyp = 'C3D15'
!
!         3D - C3D20
          case(27)
!
              eltyp = 'C3D20'
!
          end select
!
      end if
!
      write(*,*) 'Element type identified as: ', eltyp    
!
      return
      end subroutine findelementtype
!
!
!     Contains the element-by-element initializations
!     Done only once at the first run
      subroutine initialize_gradientoperators
       use globalvariables, only : numel, 
     + numdim, numpt, nnpel, gradip2ip, ipweights,
     + ipcoords, ipdomain, Nmat, invNmat, dNmat, dNmatc
      use userinputs, only : gndlinear
      use utilities, only: nolapinverse, inv2x2, inv3x3
      implicit none
!     Gauss point coordinates in sample reference
      real(8) :: xIP(numpt,numdim)
!     Node point coordinates
      real(8) :: xnode(nnpel,numdim)
!     Jacobian transpose
      real(8) :: JT(numdim,numdim)
!     Jacobian inverse transpose
      real(8) :: invJT(numdim,numdim)
!     Determinant of the inversion
      real(8) :: det
!     Gradient operator
      real(8) :: grad(numdim,nnpel)
!     Overall mapping
      real(8) :: grad_invN(numdim,numpt)
!     Element number and ip number
      integer :: iel, ipt
!
!
!
!
!     Loop through each element
      do iel = 1, numel
!
!
!
!         Vectorize element ipcoordinates
          xIP = 0.
          do ipt = 1, numpt
!
              xIP(ipt,1:numdim) = ipcoords(iel,ipt,1:numdim)
!
          end do
!
!
!
!
!
          xnode = matmul(invNmat,xIP)
!
!
!
!
!
!         For each integration point
          do ipt = 1, numpt
!
!
!             Calculate jacobian (transpose) for isoparametric space to real reference
              JT = matmul(dNmat(ipt,1:numdim,1:nnpel),xnode)
!
!
!
!
!             Invert the Jacobian
              if (numdim==2) then
!
                  call inv2x2(JT,invJT,det)
!
              elseif (numdim==3) then
!
                  call inv3x3(JT,invJT,det)
!
              endif
!
!
!             IP domain size
              ipdomain(iel,ipt) = det * ipweights(ipt)
!
!
!             
!             Gradient operator
              grad = matmul(invJT,dNmat(ipt,1:numdim,1:nnpel))
!
!
!
!             Overall mapping
              grad_invN = matmul(grad,invNmat)
!
!
!
!
!             Store
              gradip2ip(iel,ipt,1:numdim,1:numpt) = grad_invN
!
!             
!
!
!
          end do
!
!          write(*,*) 'vol: ', sum(ipdomain(iel,:))
!          write(*,*) '******************'
!
!         For the center of the element
!
!         Calculate jacobian (transpose) for isoparametric space to real reference
          JT = matmul(dNmatc,xnode)
!
!         Invert the Jacobian
          if (numdim==2) then
!
              call inv2x2(JT,invJT,det)
!
          elseif (numdim==3) then
!
              call inv3x3(JT,invJT,det)
!
          endif
!
!
!         Gradient operator
          grad = matmul(invJT,dNmatc)
!
!
!         Overall mapping
          grad_invN = matmul(grad, invNmat)
!
!         Store
          gradip2ip(iel,numpt+1,1:numdim,1:numpt) = grad_invN
!
!
!
!
!
      end do    
!
!
      return
      end subroutine initialize_gradientoperators
!
!
!
!     Finite element shape functions and derivatives
!     For Element type: CPE3 or CPS3
      subroutine CP3_N_dN(nnpel,numdim,g,h,N,dN)
      implicit none
!     inputs
!     Number of nodes per element
      integer, intent(in) :: nnpel
!     Dimensions of the analysis
      integer, intent(in) :: numdim
!     Parametric coordinates
      real(8), intent(in) :: g
      real(8), intent(in) :: h
!
!     outputs
!     Shape functions
      real(8), dimension(nnpel), intent(out) :: N
!     Shape function derivatives
      real(8), dimension(numdim,nnpel), intent(out) :: dN
!
!
!     Shape functions for CP3 element
      N(1) = 1. - g - h
!
      N(2) = g
!
      N(3) = h
!
!
!
!
!
!     Shape function derivatives wrt natural coords for CP3 element
!     dN_dg
      dN(1,1) = -1.
!
      dN(1,2) = 1.
!
      dN(1,3) = 0.
!
!     dN_dh      
      dN(2,1) = -1.
!
      dN(2,2) = 0.
!
      dN(2,3) = 1.
!
!
!
!
      return
      end subroutine CP3_N_dN
!
!
!
!
!
!
!
!
!
!     Finite element shape functions and derivatives
!     For Element type: CPE4 or CPS4 or CPE8R or CPS8R
      subroutine CP4_N_dN(nnpel,numdim,g,h,N,dN)
      implicit none
!     inputs
!     Number of nodes per element
      integer, intent(in) :: nnpel
!     Dimensions of the analysis
      integer, intent(in) :: numdim
!     Parametric coordinates
      real(8), intent(in) :: g
      real(8), intent(in) :: h
!
!     outputs
!     Shape functions
      real(8), dimension(nnpel), intent(out) :: N
!     Shape function derivatives
      real(8), dimension(numdim,nnpel), intent(out) :: dN
!
!
!     Shape functions for CP4 element
      N(1) = (1. - g) * (1. - h) / 4.
!
      N(2) = (1. + g) * (1. - h) / 4.
!
      N(3) = (1. + g) * (1. + h) / 4.
!
      N(4) = (1. - g) * (1. + h) / 4.
!
!
!
!
!     Shape function derivatives wrt natural coords for CP4 element
!
!     dN_dg
      dN(1,1) = -1. * (1. - h) / 4.
!
      dN(1,2) = 1. * (1. - h) / 4.
!
      dN(1,3) = 1. * (1. + h) / 4.
!
      dN(1,4) = -1. * (1. + h) / 4.
!
!     dN_dh
      dN(2,1) = (1. - g) * -1. / 4.
!
      dN(2,2) = (1. + g) * -1. / 4.
!
      dN(2,3) = (1. + g) * 1. / 4.
 
      dN(2,4) = (1. - g) * 1. / 4.
!
!
!
!
!
      return
      end subroutine CP4_N_dN
!
!
!
!
!
!
!
!     Finite element shape functions and derivatives
!     For Element type: CPE6 or CPS6
      subroutine CP6_N_dN(nnpel,numdim,g,h,N,dN)
      implicit none
!     inputs
!     Number of nodes per element
      integer, intent(in) :: nnpel
!     Dimensions of the analysis
      integer, intent(in) :: numdim
!     Parametric coordinates
      real(8), intent(in) :: g
      real(8), intent(in) :: h
!
!     outputs
!     Shape functions
      real(8), dimension(nnpel), intent(out) :: N
!     Shape function derivatives
      real(8), dimension(numdim,nnpel), intent(out) :: dN
!
!
!     Shape functions for CP6 element
      N(1) = 2. * (0.5 - g - h) * (1. - g - h)
!
      N(2) = 2. * g * (g - 0.5)
!
      N(3) = 2. * h * (h - 0.5)
!
      N(4) = 4. * g * (1. - g - h)
!
      N(5) = 4. * g * h
!
      N(6) = 4. * h * (1. - g - h)
!
!
!
!
!     Shape function derivatives wrt natural coords for C3D4 element
!     dN_dg
      dN(1,1) = 2. * (-1.) * (1. - g - h) + 2. * (0.5 - g - h) * (-1.)
!
      dN(1,2) = 2. * (1.) * (g - 0.5) + 2. * g * (1.)
!
      dN(1,3) = 0.
!
      dN(1,4) = 4. * (1.) * (1. - g - h) + 4. * g * (-1.)
!
      dN(1,5) = 4. * (1.) * h
!
      dN(1,6) = 4. * h * (-1.)
!
!     dN_dh       
      dN(2,1) = 2. * (-1.) * (1. - g - h) + 2. * (0.5 - g - h) * (-1.)
!
      dN(2,2) = 0.
!
      dN(2,3) = 2. * (1.) * (h - 0.5)  + 2. * h * (1.)
!
      dN(2,4) = 4. * g * (-1.)
!
      dN(2,5) = 4. * g * (1.)
!
      dN(2,6) = 4. * (1.) * (1. - g - h) + 4. * h * (-1.)
!
!
!
!
      return
      end subroutine CP6_N_dN
!
!
!
!
!
!     Finite element shape functions and derivatives
!     For Element type: CPE8 or CPS8
      subroutine CP8_N_dN(nnpel,numdim,g,h,N,dN)
      implicit none
!     inputs
!     Number of nodes per element
      integer, intent(in) :: nnpel
!     Dimensions of the analysis
      integer, intent(in) :: numdim
!     Parametric coordinates
      real(8), intent(in) :: g
      real(8), intent(in) :: h
!
!     outputs
!     Shape functions
      real(8), dimension(nnpel), intent(out) :: N
!     Shape function derivatives
      real(8), dimension(numdim,nnpel), intent(out) :: dN
!
!
!     Shape functions for CP8 element
      N(1) = -1./4. * (1. - g) * (1. - h) * (1. + g + h)
!
      N(2) = -1./4. * (1. + g) * (1. - h) * (1. - g + h)
!
      N(3) = -1./4. * (1. + g) * (1. + h) * (1. - g - h)
!
      N(4) = -1./4. * (1. - g) * (1. + h) * (1. + g - h)
!
      N(5) = 1./2. * (1. - g) * (1. + g) * (1. - h)
!
      N(6) = 1./2. * (1. - h) * (1. + h) * (1. + g)
!
      N(7) = 1./2. * (1. - g) * (1. + g) * (1. + h)
!
      N(8) = 1./2. * (1. - h) * (1. + h) * (1. - g)
!
!
!
!
!
!     Shape function derivatives wrt natural coords for CP8 element
!     dN_dg
      dN(1,1) = (-1./4.) * (-1.) * (1. - h) * (1. + g + h) +
     + (-1./4.) * (1. - g) * (1. - h) * (1.)
!
      dN(1,2) = (-1./4.) * (1.) *  (1. - h) * (1. - g + h) +
     + (-1./4.) * (1. + g) * (1. - h) * (-1.)
!
      dN(1,3) = (-1./4.) * (1.) *  (1. + h) * (1. - g - h) +
     + (-1./4.) * (1. + g) * (1. + h) * (-1.)
!
      dN(1,4) = (-1./4.) * (-1.) * (1. + h) * (1. + g - h) +
     + (-1./4.) * (1. - g) * (1. + h) * (1.)
!
      dN(1,5) = 1./2. * (-1.) * (1. + g) * (1. - h) +
     + 1./2. * (1. - g) * (1.) * (1. - h)
!
      dN(1,6) = 1./2. * (1. - h) * (1. + h) * (1.)
!
      dN(1,7) = 1./2. * (-1.) * (1. + g) * (1. + h) +
     + 1./2. * (1. - g) * (1.) * (1. + h)
!
      dN(1,8) = 1./2. * (1. - h) * (1. + h) * (-1.)
!
!
!
!     dN_dh
      dN(2,1) = (-1./4.) * (1. - g) * (-1.) * (1. + g + h)  +
     + (-1./4.) * (1. - g) * (1. - h) * (1.)
!
      dN(2,2) = (-1./4.) * (1. + g) * (-1.) * (1. - g + h) +
     + (-1./4.) * (1. + g) * (1. - h) * (1.)
!
      dN(2,3) = (-1./4.) * (1. + g) *  (1.) * (1. - g - h) + 
     + (-1./4.) * (1. + g) * (1. + h) * (-1.)
 
      dN(2,4) = (-1./4.) * (1. - g) *  (1.) * (1. + g - h) +
     + (-1./4.) * (1. - g) * (1. + h) * (-1.)
!
      dN(2,5) = 1./2. * (1. - g) * (1. + g) * (-1.)
!
      dN(2,6) = 1./2. * (-1.) * (1. + h) * (1. + g) +
     + 1./2. * (1. - h) * (1.) * (1. + g)
!
      dN(2,7) = 1./2. * (1. - g) * (1. + g) * (1.)
!
      dN(2,8) = 1./2. * (-1.) * (1. + h) * (1. - g) +
     + 1./2. * (1. - h) * (1.) * (1. - g)
!
!
!
!
!          
!     
      return
      end subroutine CP8_N_dN
!
!
!
!
!
!     Finite element shape functions and derivatives
!     For Element type: C3D4
      subroutine C3D4_N_dN(nnpel,numdim,g,h,r,N,dN)
      implicit none
!     inputs
!     Number of nodes per element
      integer, intent(in) :: nnpel
!     Dimensions of the analysis
      integer, intent(in) :: numdim
!     Parametric coordinates
      real(8), intent(in) :: g
      real(8), intent(in) :: h
      real(8), intent(in) :: r
!
!     outputs
!     Shape functions
      real(8), dimension(nnpel), intent(out) :: N
!     Shape function derivatives
      real(8), dimension(numdim,nnpel), intent(out) :: dN
!
!
!     Shape functions for CP3 element
      N(1) = 1. - h - r - g
!
      N(2) = g
!
      N(3) = h
!
      N(4) = r
!
!
!
!
!
!     Shape function derivatives wrt natural coords for CP3 element
!     dN_dg
      dN(1,1) = -1.
!
      dN(1,2) = 1.
!
      dN(1,3) = 0.
!
      dN(1,4) = 0.
!
!
!     dN_dh      
      dN(2,1) = -1.
!
      dN(2,2) = 0.
!
      dN(2,3) = 1.
!
      dN(2,4) = 0.
!
!
!     dN_dr      
      dN(3,1) = -1.
!
      dN(3,2) = 0.
!
      dN(3,3) = 0.
!
      dN(3,4) = 1.
!
!
!
!
!
!
      return
      end subroutine C3D4_N_dN
!
!
!
!
!
!     Finite element shape functions and derivatives
!     For Element type: C3D6
      subroutine C3D6_N_dN(nnpel,numdim,g,h,r,N,dN)
      implicit none
!     inputs
!     Number of nodes per element
      integer, intent(in) :: nnpel
!     Dimensions of the analysis
      integer, intent(in) :: numdim
!     Parametric coordinates
      real(8), intent(in) :: g
      real(8), intent(in) :: h
      real(8), intent(in) :: r
!
!     outputs
!     Shape functions
      real(8), dimension(nnpel), intent(out) :: N
!     Shape function derivatives
      real(8), dimension(numdim,nnpel), intent(out) :: dN
!
!
!     Shape functions for C3D6 element
      N(1) = 1./2.*(1.-g-h)*(1.-r)
!
      N(2) = 1./2.*g*(1.-r)
!
      N(3) = 1./2.*h*(1.-r)
!
      N(4) = 1./2.*(1.-g-h)*(1.+r)
!
      N(5) = 1./2.*g*(1.+r)
!
      N(6) = 1./2.*h*(1.+r)
!
!
!
!
!
!
!     Shape function derivatives wrt natural coords for CP3 element
!     dN_dg
      dN(1,1) = 1./2.*(-1.)*(1.-r)
!
      dN(1,2) = 1./2.*(1.)*(1.-r)
!
      dN(1,3) = 0.
!
      dN(1,4) = 1./2.*(-1.)*(1.+r)
!
      dN(1,5) = 1./2.*(1.)*(1.+r)
!
      dN(1,6) = 0.
!
!
!
!     dN_dh      
      dN(2,1) =  1./2.*(-1.)*(1.-r)
!
      dN(2,2) = 0.
!
      dN(2,3) = 1./2.*(1.)*(1.-r)
!
      dN(2,4) = 1./2.*(-1.)*(1.+r)
!
      dN(2,5) = 0.
!
      dN(2,6) = 1./2.*(1.)*(1.+r)
!
!
!
!     dN_dr      
      dN(3,1) = 1./2.*(1.-g-h)*(-1.)
!
      dN(3,2) = 1./2.*g*(-1.)
!
      dN(3,3) = 1./2.*h*(-1.)
!
      dN(3,4) = 1./2.*(1.-g-h)*(1.)
!
      dN(3,5) = 1./2.*g*(1.)
!
      dN(3,6) = 1./2.*h*(1.)
!
!
!
!
!
!
      return
      end subroutine C3D6_N_dN      
!
!
!
!
!
!
!
!     Finite element shape functions and derivatives
!     For Element type: C3D8
      subroutine C3D8_N_dN(nnpel,numdim,g,h,r,N,dN)
      implicit none
!     inputs
!     Number of nodes per element
      integer, intent(in) :: nnpel
!     Dimensions of the analysis
      integer, intent(in) :: numdim
!     Parametric coordinates
      real(8), intent(in) :: g
      real(8), intent(in) :: h
      real(8), intent(in) :: r
!
!     outputs
!     Shape functions
      real(8), dimension(nnpel), intent(out) :: N
!     Shape function derivatives
      real(8), dimension(numdim,nnpel), intent(out) :: dN
!
!
!     Shape functions for C3D8 element
      N(1) = 1./8.*(1.-g)*(1.-h)*(1.-r)
!
      N(2) = 1./8.*(1.+g)*(1.-h)*(1.-r)
!
      N(3) = 1./8.*(1.+g)*(1.+h)*(1.-r)
!
      N(4) = 1./8.*(1.-g)*(1.+h)*(1.-r)
!
      N(5) = 1./8.*(1.-g)*(1.-h)*(1.+r)
!
      N(6) = 1./8.*(1.+g)*(1.-h)*(1.+r)
!
      N(7) = 1./8.*(1.+g)*(1.+h)*(1.+r)
!
      N(8) = 1./8.*(1.-g)*(1.+h)*(1.+r)
!
!
!
!
!
!
!     Shape function derivatives wrt natural coords for C3D8 element
!     dN_dg
      dN(1,1) = 1./8.*(-1.)*(1.-h)*(1.-r)
!
      dN(1,2) = 1./8.*(1.)*(1.-h)*(1.-r)
!
      dN(1,3) = 1./8.*(1.)*(1.+h)*(1.-r)
!
      dN(1,4) = 1./8.*(-1.)*(1.+h)*(1.-r)
!
      dN(1,5) = 1./8.*(-1.)*(1.-h)*(1.+r)
!
      dN(1,6) = 1./8.*(1.)*(1.-h)*(1.+r)
!
      dN(1,7) = 1./8.*(1.)*(1.+h)*(1.+r)
!
      dN(1,8) = 1./8.*(-1.)*(1.+h)*(1.+r)
!
!
!
!
!     dN_dh
      dN(2,1) = 1./8.*(1.-g)*(-1.)*(1.-r)
!
      dN(2,2) = 1./8.*(1.+g)*(-1.)*(1.-r)
!
      dN(2,3) = 1./8.*(1.+g)*(1.)*(1.-r)
!
      dN(2,4) = 1./8.*(1.-g)*(1.)*(1.-r)
!
      dN(2,5) = 1./8.*(1.-g)*(-1.)*(1.+r)
!
      dN(2,6) = 1./8.*(1.+g)*(-1.)*(1.+r)
!
      dN(2,7) = 1./8.*(1.+g)*(1.)*(1.+r)
!
      dN(2,8) = 1./8.*(1.-g)*(1.)*(1.+r)
!
!
!     dN_dr
      dN(3,1) = 1./8.*(1.-g)*(1.-h)*(-1.)
!
      dN(3,2) = 1./8.*(1.+g)*(1.-h)*(-1.)
!
      dN(3,3) = 1./8.*(1.+g)*(1.+h)*(-1.)
!
      dN(3,4) = 1./8.*(1.-g)*(1.+h)*(-1.)
!
      dN(3,5) = 1./8.*(1.-g)*(1.-h)*(1.)
!
      dN(3,6) = 1./8.*(1.+g)*(1.-h)*(1.)
!
      dN(3,7) = 1./8.*(1.+g)*(1.+h)*(1.)
!
      dN(3,8) = 1./8.*(1.-g)*(1.+h)*(1.)
!
!
!
!
!
!
      return
      end subroutine C3D8_N_dN
!
!
!
!
!
!
!
!
!
!
!
!
!
!     Finite element shape functions and derivatives
!     For Element type: C3D10
      subroutine C3D10_N_dN(nnpel,numdim,g,h,r,N,dN)
      implicit none
!     inputs
!     Number of nodes per element
      integer, intent(in) :: nnpel
!     Dimensions of the analysis
      integer, intent(in) :: numdim
!     Parametric coordinates
      real(8), intent(in) :: g
      real(8), intent(in) :: h
      real(8), intent(in) :: r
!
!     outputs
!     Shape functions
      real(8), dimension(nnpel), intent(out) :: N
!     Shape function derivatives
      real(8), dimension(numdim,nnpel), intent(out) :: dN
!
!
!     Shape functions for C3D10 element
      N(1) = (2.*(1.-g-h-r)-1.)*(1.-g-h-r)
!
      N(2) = (2.*g - 1.)*g
!
      N(3) = (2.*h - 1.)*h
!
      N(4) = (2.*r - 1.)*r
!
      N(5) = 4.*(1.-g-h-r)*g
!
      N(6) = 4.*g*h
!
      N(7) = 4.*(1.-g-h-r)*h
!
      N(8) = 4.*(1.-g-h-r)*r
!
      N(9) = 4.*g*r
!
      N(10) = 4.*h*r
!
!
!
!
!
!     Shape function derivatives wrt natural coords for C3D10 element
!     dN_dg
      dN(1,1) = (2.*(-1.))*(1.-g-h-r) + (2.*(1.-g-h-r)-1.)*(-1.)
!
      dN(1,2) = (2.*(1.))*g + (2.*g - 1.)*(1.)
!
      dN(1,3) = 0.
!
      dN(1,4) = 0.
!
      dN(1,5) = 4.*(-1.)*g + 4.*(1.-g-h-r)*(1.)
!
      dN(1,6) = 4.*h
!
      dN(1,7) = 4.*(-1.)*h
!
      dN(1,8) = 4.*(-1.)*r
!
      dN(1,9) = 4.*(1.)*r
!
      dN(1,10) = 0.
!
!
!
!
!     dN_dh
      dN(2,1) = (2.*(-1.))*(1.-g-h-r) + (2.*(1.-g-h-r)-1.)*(-1.)
!
      dN(2,2) = 0.
!
      dN(2,3) = (2.*(1.))*h + (2.*h - 1.)*(1.)
!
      dN(2,4) = 0.
!
      dN(2,5) = 4.*(-1.)*g
!
      dN(2,6) = 4.*g*(1.)
!
      dN(2,7) = 4.*(-1.)*h + 4.*(1.-g-h-r)*(1.)
!
      dN(2,8) = 4.*(-1.)*r
!
      dN(2,9) = 0.
!
      dN(2,10) = 4.*(1.)*r
!
!
!
!     dN_dr
      dN(3,1) = (2.*(-1.))*(1.-g-h-r) + (2.*(1.-g-h-r)-1.)*(-1.)
!
      dN(3,2) = 0.
!
      dN(3,3) = 0.
!
      dN(3,4) = (2.*(1.))*r + (2.*r - 1.)*(1.)
!
      dN(3,5) = 4.*(-1.)*g
!
      dN(3,6) = 0.
!
      dN(3,7) = 4.*(-1.)*h
!
      dN(3,8) = 4.*(-1.)*r + 4.*(1.-g-h-r)*(1.)
!
      dN(3,9) = 4.*g*(1.)
!
      dN(3,10) = 4.*h*(1.)
!
!
!
!
!
!
!
      return
      end subroutine C3D10_N_dN
!
!
!
!
!
!     Finite element shape functions and derivatives
!     For Element type: C3D15
      subroutine C3D15_N_dN(nnpel,numdim,g,h,r,N,dN)
      implicit none
!     inputs
!     Number of nodes per element
      integer, intent(in) :: nnpel
!     Dimensions of the analysis
      integer, intent(in) :: numdim
!     Parametric coordinates
      real(8), intent(in) :: g
      real(8), intent(in) :: h
      real(8), intent(in) :: r
!
!     outputs
!     Shape functions
      real(8), dimension(nnpel), intent(out) :: N
!     Shape function derivatives
      real(8), dimension(numdim,nnpel), intent(out) :: dN
!
!
!     Shape functions for C3D15 element
      N(1) = 1./2.*((1.-g-h)*(2.*(1.-g-h)-1.)*(1.-r)-(1.-g-h)*(1.-r**2))
!
      N(2) = 1./2.*(g*(2.*g-1.)*(1.-r)-g*(1.-r**2))
!
      N(3) = 1./2.*(h*(2.*h-1.)*(1.-r)-h*(1.-r**2))
!
      N(4) = 1./2.*((1.-g-h)*(2.*(1.-g-h)-1.)*(1.+r)-(1.-g-h)*(1.-r**2))
!
      N(5) = 1./2.*(g*(2.*g-1.)*(1.+r)-g*(1.-r**2))
!
      N(6) = 1./2.*(h*(2.*h-1.)*(1.+r)-h*(1.-r**2))
!
      N(7) = 2.*(1.-g-h)*g*(1.-r)
!
      N(8) = 2.*g*h*(1.-r)
!
      N(9) = 2.*h*(1.-g-h)*(1.-r)
!
      N(10) = 2.*(1.-g-h)*g*(1.+r)
!
      N(11) = 2.*g*h*(1.+r)
!
      N(12) = 2.*h*(1.-g-h)*(1.+r)
!
      N(13) = (1.-g-h)*(1.-r**2)
!
      N(14) = g*(1.-r**2)
!
      N(15) = h*(1.-r**2)
!
!
!
!
!
!
!     Shape function derivatives wrt natural coords for C3D15 element
!     dN_dg
      dN(1,1) = -((r - 1.)*(4.*g + 4.*h + r - 2.))/2.
!
      dN(1,2) = ((r - 1.)*(r - 4.*g + 2.))/2.
!
      dN(1,3) = 0.
!
      dN(1,4) = ((r + 1.)*(4.*g + 4.*h - r - 2.))/2.
!
      dN(1,5) = ((r + 1.)*(4.*g + r - 2.))/2.
!
      dN(1,6) = 0.
!
      dN(1,7) = 2.*(r - 1.)*(2.*g + h - 1.)
!
      dN(1,8) = -2.*h*(r - 1.)
!
      dN(1,9) = 2.*h*(r - 1.)
!
      dN(1,10) = -2.*(r + 1.)*(2.*g + h - 1.)
!
      dN(1,11) = 2.*h*(r + 1.)
!
      dN(1,12) = -2.*h*(r + 1.)
!
      dN(1,13) = r**2 - 1.
!
      dN(1,14) = 1. - r**2
!
      dN(1,15) = 0.
!
!
!
!     dN_dh
      dN(2,1) = -((r - 1.)*(4.*g + 4.*h + r - 2.))/2.
!
      dN(2,2) = 0.
!
      dN(2,3) = ((r - 1.)*(r - 4.*h + 2.))/2.
!
      dN(2,4) = ((r + 1.)*(4.*g + 4.*h - r - 2.))/2.
!
      dN(2,5) = 0.
!
      dN(2,6) = ((r + 1.)*(4.*h + r - 2.))/2.
!
      dN(2,7) = 2.*g*(r - 1.)
!
      dN(2,8) = -2.*g*(r - 1.)
!
      dN(2,9) = 2.*(r - 1.)*(g + 2.*h - 1.)
!
      dN(2,10) = -2.*g*(r + 1.)
!
      dN(2,11) = 2.*g*(r + 1.)
!
      dN(2,12) = -2.*(r + 1.)*(g + 2.*h - 1.)
!
      dN(2,13) = r**2 - 1.
!
      dN(2,14) = 0.
!
      dN(2,15) = 1. - r**2
!
!
!
!     dN_dr
      dN(3,1) = -((g + h - 1.)*(2.*g + 2.*h + 2.*r - 1.))/2.
!
      dN(3,2) = (g*(2.*r - 2.*g + 1.))/2.
!
      dN(3,3) = (h*(2.*r - 2.*h + 1.))/2.
!
      dN(3,4) = ((g + h - 1.)*(2.*g + 2.*h - 2.*r - 1.))/2.
!
      dN(3,5) = (g*(2.*g + 2.*r - 1.))/2.
!
      dN(3,6) = (h*(2.*h + 2.*r - 1.))/2.
!
      dN(3,7) = g*(2.*g + 2.*h - 2.)
!
      dN(3,8) = -2.*g*h
!
      dN(3,9) = 2.*h*(g + h - 1.)
!
      dN(3,10) = -g*(2.*g + 2.*h - 2.)
!
      dN(3,11) = 2.*g*h
!
      dN(3,12) = -2.*h*(g + h - 1.)
!
      dN(3,13) = 2.*r*(g + h - 1.)
!
      dN(3,14) = -2.*g*r
!
      dN(3,15) = -2.*h*r
!
!
!
!
!
!
!
!
      return
      end subroutine C3D15_N_dN
!
!
!
!     Finite element shape functions and derivatives
!     For Element type: C3D20
      subroutine C3D20_N_dN(nnpel,numdim,g,h,r,N,dN)
      implicit none
!     inputs
!     Number of nodes per element
      integer, intent(in) :: nnpel
!     Dimensions of the analysis
      integer, intent(in) :: numdim
!     Parametric coordinates
      real(8), intent(in) :: g
      real(8), intent(in) :: h
      real(8), intent(in) :: r
!
!     outputs
!     Shape functions
      real(8), dimension(nnpel), intent(out) :: N
!     Shape function derivatives
      real(8), dimension(numdim,nnpel), intent(out) :: dN
!
!
!     Shape functions for C3D20 element
      N(1) = (g/8. - 1./8.)*(h - 1.)*(r - 1.)*(g + h + r + 2.)
!
      N(2) = -(g/8. + 1./8.)*(h - 1.)*(r - 1.)*(h - g + r + 2.)
!
      N(3) = -(g/8. + 1./8.)*(h + 1.)*(r - 1.)*(g + h - r - 2.)
!
      N(4) = -(g/8. - 1./8.)*(h + 1.)*(r - 1.)*(g - h + r + 2.)
!
      N(5) = -(g/8. - 1./8.)*(h - 1.)*(r + 1.)*(g + h - r + 2.)
!
      N(6) = -(g/8. + 1./8.)*(h - 1.)*(r + 1.)*(g - h + r - 2.)
!
      N(7) = (g/8. + 1./8.)*(h + 1.)*(r + 1.)*(g + h + r - 2.)
!
      N(8) = (g/8. - 1./8.)*(h + 1.)*(r + 1.)*(g - h - r + 2.)
!
      N(9) = -(g/4. - 1./4.)*(g + 1.)*(h - 1.)*(r - 1.)
!
      N(10) = (h/4. - 1./4.)*(g + 1.)*(h + 1.)*(r - 1.)
!
      N(11) = (g/4. - 1./4.)*(g + 1.)*(h + 1.)*(r - 1.)
!
      N(12) = -(h/4. - 1./4.)*(g - 1.)*(h + 1.)*(r - 1.)
!
      N(13) = (g/4. - 1./4.)*(g + 1.)*(h - 1.)*(r + 1.)
!
      N(14) = -(h/4. - 1./4.)*(g + 1.)*(h + 1.)*(r + 1.)
!
      N(15) = -(g/4. - 1./4.)*(g + 1.)*(h + 1.)*(r + 1.)
!
      N(16) = (h/4. - 1./4.)*(g - 1.)*(h + 1.)*(r + 1.)
!
      N(17) = -(r/4. - 1./4.)*(g - 1.)*(h - 1.)*(r + 1.)
!
      N(18) = (r/4. - 1./4.)*(g + 1.)*(h - 1.)*(r + 1.)
!
      N(19) = -(r/4. - 1./4.)*(g + 1.)*(h + 1.)*(r + 1.)
!
      N(20) = (r/4. - 1./4.)*(g - 1.)*(h + 1.)*(r + 1.)
!
!
!
!
!
!     Shape function derivatives wrt natural coords for C3D20 element
!     dN_dg
      dN(1,1) = (g/8. - 1./8.)*(h - 1.)*(r - 1.) +
     + (h - 1.)*(r - 1.)*(g + h + r + 2.)/8.
!
      dN(1,2) = (g/8. + 1./8.)*(h - 1.)*(r - 1.) -
     + (h - 1.)*(r - 1.)*(h - g + r + 2.)/8.
!
      dN(1,3) = - (g/8. + 1./8.)*(h + 1.)*(r - 1.) -
     + (h + 1.)*(r - 1.)*(g + h - r - 2.)/8.
!
      dN(1,4) = - (g/8. - 1./8.)*(h + 1.)*(r - 1.) -
     + (h + 1.)*(r - 1.)*(g - h + r + 2.)/8.
!
      dN(1,5) = - (g/8. - 1./8.)*(h - 1.)*(r + 1.) -
     + (h - 1.)*(r + 1.)*(g + h - r + 2.)/8.
!
      dN(1,6) = - (g/8. + 1./8.)*(h - 1.)*(r + 1.) -
     + (h - 1.)*(r + 1.)*(g - h + r - 2.)/8.
!
      dN(1,7) = (g/8. + 1./8.)*(h + 1.)*(r + 1.) +
     + (h + 1.)*(r + 1.)*(g + h + r - 2.)/8.
!
      dN(1,8) = (g/8. - 1./8.)*(h + 1.)*(r + 1.) +
     + (h + 1.)*(r + 1.)*(g - h - r + 2.)/8.
!
      dN(1,9) = - (g/4. - 1./4.)*(h - 1.)*(r - 1.) -
     + (g + 1.)*(h - 1.)*(r - 1.)/4.
!
      dN(1,10) = (h/4. - 1./4.)*(h + 1.)*(r - 1.)

      dN(1,11) = (g/4. - 1./4.)*(h + 1.)*(r - 1.) +
     + (g + 1.)*(h + 1.)*(r - 1.)/4.
!
      dN(1,12) = -(h/4. - 1./4.)*(h + 1.)*(r - 1.)
!
      dN(1,13) = (g/4. - 1./4.)*(h - 1.)*(r + 1.) +
     + (g + 1.)*(h - 1.)*(r + 1.)/4.
!
      dN(1,14) = -(h/4. - 1./4.)*(h + 1.)*(r + 1.)
!
      dN(1,15) = - (g/4. - 1./4.)*(h + 1.)*(r + 1.) -
     + (g + 1.)*(h + 1.)*(r + 1.)/4.
!
      dN(1,16) = (h/4. - 1./4.)*(h + 1.)*(r + 1.)
!
      dN(1,17) = -(r/4. - 1./4.)*(h - 1.)*(r + 1.)
!
      dN(1,18) = (r/4. - 1./4.)*(h - 1.)*(r + 1.)
!
      dN(1,19) = -(r/4. - 1./4.)*(h + 1.)*(r + 1.)
!
      dN(1,20) = (r/4. - 1./4.)*(h + 1.)*(r + 1.)
!
!
!
!     dN_dh
      dN(2,1) = (g/8. - 1./8.)*(h - 1.)*(r - 1.) +
     + (g/8. - 1./8.)*(r - 1.)*(g + h + r + 2.)
!
      dN(2,2) = - (g/8. + 1./8.)*(h - 1.)*(r - 1.) -
     + (g/8. + 1./8.)*(r - 1.)*(h - g + r + 2.)
!
      dN(2,3) = - (g/8. + 1./8.)*(h + 1.)*(r - 1.) -
     + (g/8. + 1./8.)*(r - 1.)*(g + h - r - 2.)
!
      dN(2,4) = (g/8. - 1./8.)*(h + 1.)*(r - 1.) -
     + (g/8. - 1./8.)*(r - 1.)*(g - h + r + 2.)
!
      dN(2,5) = - (g/8. - 1./8.)*(h - 1.)*(r + 1.) -
     + (g/8. - 1./8.)*(r + 1.)*(g + h - r + 2.)
!
      dN(2,6) = (g/8. + 1./8.)*(h - 1.)*(r + 1.) -
     + (g/8. + 1./8.)*(r + 1.)*(g - h + r - 2.)
!
      dN(2,7) = (g/8. + 1./8.)*(h + 1.)*(r + 1.) +
     + (g/8. + 1./8.)*(r + 1.)*(g + h + r - 2.)
!
      dN(2,8) = (g/8. - 1./8.)*(r + 1.)*(g - h - r + 2.) -
     + (g/8. - 1./8.)*(h + 1.)*(r + 1.)
!
      dN(2,9) = -(g/4. - 1./4.)*(g + 1.)*(r - 1.)
!
      dN(2,10) = (h/4. - 1./4.)*(g + 1.)*(r - 1.) +
     + (g + 1.)*(h + 1.)*(r - 1.)/4.
!
      dN(2,11) = (g/4. - 1./4.)*(g + 1.)*(r - 1.)
!
      dN(2,12) = - (h/4. - 1./4.)*(g - 1.)*(r - 1.) -
     + (g - 1.)*(h + 1.)*(r - 1.)/4.
!
      dN(2,13) = (g/4. - 1./4.)*(g + 1.)*(r + 1.)
!
      dN(2,14) = - (h/4. - 1./4.)*(g + 1.)*(r + 1.) -
     + (g + 1.)*(h + 1.)*(r + 1.)/4.
!
      dN(2,15) = -(g/4. - 1./4.)*(g + 1.)*(r + 1.)

      dN(2,16) = (h/4. - 1./4.)*(g - 1.)*(r + 1.) +
     + (g - 1.)*(h + 1.)*(r + 1.)/4.
!
      dN(2,17) = -(r/4. - 1./4.)*(g - 1.)*(r + 1.)
!
      dN(2,18) = (r/4. - 1./4.)*(g + 1.)*(r + 1.)
!
      dN(2,19) = -(r/4. - 1./4.)*(g + 1.)*(r + 1.)
!
      dN(2,20) = (r/4. - 1./4.)*(g - 1.)*(r + 1.)
!
!
!
!     dN_dr
      dN(3,1) = (g/8. - 1./8.)*(h - 1.)*(r - 1.) +
     + (g/8. - 1./8.)*(h - 1.)*(g + h + r + 2.)
!
      dN(3,2) = - (g/8. + 1./8.)*(h - 1.)*(r - 1.) -
     + (g/8. + 1./8.)*(h - 1.)*(h - g + r + 2.)
!
      dN(3,3) = (g/8. + 1./8.)*(h + 1.)*(r - 1.) -
     + (g/8. + 1./8.)*(h + 1.)*(g + h - r - 2.)
!
      dN(3,4) = - (g/8. - 1./8.)*(h + 1.)*(r - 1.) -
     + (g/8. - 1./8.)*(h + 1.)*(g - h + r + 2.)
!
      dN(3,5) = (g/8. - 1./8.)*(h - 1.)*(r + 1.) -
     + (g/8. - 1./8.)*(h - 1.)*(g + h - r + 2.)
!
      dN(3,6) = - (g/8. + 1./8.)*(h - 1.)*(r + 1.) -
     + (g/8. + 1./8.)*(h - 1.)*(g - h + r - 2.)
!
      dN(3,7) = (g/8. + 1./8.)*(h + 1.)*(r + 1.) +
     + (g/8. + 1./8.)*(h + 1.)*(g + h + r - 2.)
!
      dN(3,8) = (g/8. - 1./8.)*(h + 1.)*(g - h - r + 2.) -
     + (g/8. - 1./8.)*(h + 1.)*(r + 1.)
!
      dN(3,9) = -(g/4. - 1./4.)*(g + 1.)*(h - 1.)
!
      dN(3,10) = (h/4. - 1./4.)*(g + 1.)*(h + 1.)
!
      dN(3,11) = (g/4. - 1./4.)*(g + 1.)*(h + 1.)
!
      dN(3,12) = -(h/4. - 1./4.)*(g - 1.)*(h + 1.)
!
      dN(3,13) = (g/4. - 1./4.)*(g + 1.)*(h - 1.)
!
      dN(3,14) = -(h/4. - 1./4.)*(g + 1.)*(h + 1.)
!
      dN(3,15) = -(g/4. - 1./4.)*(g + 1.)*(h + 1.)
!
      dN(3,16) = (h/4. - 1./4.)*(g - 1.)*(h + 1.)
!
      dN(3,17) = - (r/4. - 1./4.)*(g - 1.)*(h - 1.) -
     + (g - 1.)*(h - 1.)*(r + 1.)/4.
!
      dN(3,18) = (r/4. - 1./4.)*(g + 1.)*(h - 1.) +
     + (g + 1.)*(h - 1.)*(r + 1.)/4.
!
      dN(3,19) = - (r/4. - 1./4.)*(g + 1.)*(h + 1.) -
     + (g + 1.)*(h + 1.)*(r + 1.)/4.
!
      dN(3,20) = (r/4. - 1./4.)*(g - 1.)*(h + 1.) +
     + (g - 1.)*(h + 1.)*(r + 1.)/4.
!
!
!
!
!
!
!
!
      return
      end subroutine C3D20_N_dN
!
!
!
!
!
!
!
!
!
!
      end module meshprop