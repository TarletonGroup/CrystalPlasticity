c 3D UEL for ABAQUS Code (Compatible with 8 node Brick Elements)
c By: Daniel Spring
c October 3rd, 2012
c
c References
c 1) Seong Hyeok Song "Fracture of Asphalt Concrete: A Cohesive
c Zone Modeling Approach Considering Viscoelastic Effects." PhD Thesis,
c Department of Civil and Environmental Engineering, UIUC, 2006.
c
c 2) Kyoungsoo Park. "Concrete Fracture Mechanics and Size Effect Using
c a Specialized Cohesive Zone Model" MS Thesis, Department of Civil and
c Environmental Engineering, UIUC, 2005.
c
c 3) "Cohesive fracture model for functionally graded fiber reinforced
c concrete" K. Park, G.H. Paulino, J. Roesler. Cement and Concrete
c Research. Vol 40, No. 6, pp. 956-965, 2010.
c
c 4) "A growing library of three-dimensional cohesive elements for 
c use in ABAQUS" D. W. Spring, G. H. Paulino. Engineering Fracture 
c Mechanics. Volume 126, August 2014, Pages 190-216
c
c Bilinear force-separation law included
c Nicolo Grilli
c University of Oxford 
c 7 Luglio 2020
c
c The damage affects the stiffness tensor
c so that damaged stiffness is (1-D)K
c Introducing the off-diagonal Jacobian terms Dnt
c damage in shear included
c using the effetive displacement variable delta
c in Ortiz, Pandolfi, 1999
c
c =====================================================================
      SUBROUTINE UEL (RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
     1 PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,
     2 DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG,
     3 PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS,
     4 NJPROP, PERIOD)
c
      INCLUDE 'ABA_PARAM.INC'
c
      include 'mycommon.f'
c     
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
c
      DIMENSION ds1(4),ds2(4),dn(4),Trac(MCRD,NRHS),
     1 Trac_Jacob(MCRD,MCRD),R(MCRD,MCRD),coord_l(MCRD,NNODE),
     2 GP_coord(2),sf(4),B(MCRD,NDOFEL),co_de_m(3,4),
     3 B_t(NDOFEL,MCRD), Transformation_M(NDOFEL,NDOFEL),
     4 Transformation_M_T(NDOFEL,NDOFEL),temp1(MCRD,NDOFEL)
c
      DIMENSION stiff_l(NDOFEL,NDOFEL),temp2(NDOFEL,NDOFEL),
     1 stiff_g(NDOFEL,NDOFEL),residual_l(NDOFEL,NRHS),
     2 residual_g(NDOFEL,NRHS),aJacob_M(2,3),delu_loc_gp(mcrd),
     3 co_de(mcrd,nnode),damage(4)
c
      DOUBLE PRECISION delta0, deltaF, stiffK, stiffG, sigma0

c Normal and shear distance
      DOUBLE PRECISION opn, opt, deltaortiz

c Map between cohesive elements and neighbouring bulk elements
      real*8 :: cohelebulkmap(nElemUel,3)

c Twin volume fractions in the neighbouring elements
      real*8 :: tvfneigh1(2), tvfneigh2(2)
      integer :: tempelem, tempneighelem1, tempneighelem2

      include 'CoheleBulkMap.f'

      tempelem = JELEM - nElements ! number of the cohesive element

c Find indices of the neighbouring bulk elements
c and corresponding twin volume fraction
      tempneighelem1 = cohelebulkmap(tempelem,2)
      tempneighelem2 = cohelebulkmap(tempelem,3)

      tvfneigh1(1:2) = 0.0
      tvfneigh2(1:2) = 0.0

      do i = 1,nintpts
        tvfneigh1(1) = tvfneigh1(1) + kTwinVolFrac(tempneighelem1,i,1)
        tvfneigh1(2) = tvfneigh1(2) + kTwinVolFrac(tempneighelem1,i,2)
      end do
      tvfneigh1(1) = tvfneigh1(1) / nintpts
      tvfneigh1(2) = tvfneigh1(2) / nintpts

      do i = 1,nintpts
        tvfneigh2(1) = tvfneigh2(1) + kTwinVolFrac(tempneighelem2,i,1)
        tvfneigh2(2) = tvfneigh2(2) + kTwinVolFrac(tempneighelem2,i,2)
      end do
      tvfneigh2(1) = tvfneigh2(1) / nintpts
      tvfneigh2(2) = tvfneigh2(2) / nintpts

c
c Define Inputs========================================================
c 

      stiffK = PROPS(1) ! normal stiffness
      stiffG = PROPS(2) ! shear stiffness
      sigma0 = PROPS(3) ! max stress before failure

      do j = 1,nnode
        if (coords(1,j) .GT. 25.0) then
          sigma0 = sigma0 * (1.0 + ((coords(1,j)-25.0)/5.0))
          EXIT
        end if
        if (coords(1,j) .LT. 5.0) then
          sigma0 = sigma0 * (1.0 + ((5.0-coords(1,j))/5.0))
          EXIT
        end if
      end do

      sigma0 = sigma0 - 400.0*max(abs(tvfneigh1(1)-tvfneigh2(1)),abs(tvfneigh1(2)-tvfneigh2(2)))

      ! output sigma0 in the bulk elements as the
      ! minimum on the four interface elements connected to
      ! that element
      call MutexLock( 12 )      ! lock Mutex #12
      do i = 1,nintpts   
        if (sigma0 < kSigma0(tempneighelem1,i)) then
          kSigma0(tempneighelem1,i) = sigma0
        end if
        if (sigma0 < kSigma0(tempneighelem2,i)) then
          kSigma0(tempneighelem2,i) = sigma0
        end if
      end do   
      call MutexUnlock( 12 )   ! unlock Mutex #12

      delta0 = sigma0/stiffK ! max displacement before failure begins
      deltaF = PROPS(4) ! displacement for complete failure

! 4 Gauss points at the surface of the cohesive element
      GP_n=4d0

c
c Initialize Matrices and Vectors======================================
c 
      ! shear and normal openings at the four nodes
      ! of the midsurface
      call k_vector_zero(ds1,4) ! ds1 has dimension 4
      call k_vector_zero(ds2,4) ! ds2 has dimension 4
      call k_vector_zero(dn,4) ! dn has dimension 4
      ! traction vector and derivatives
      ! with respect to opening displacement vector
      call k_matrix_zero(Trac,mcrd,nrhs) ! nrhs = 1 (apart from Riks analysis)
      call k_matrix_zero(Trac_Jacob,mcrd,mcrd) ! Trac_Jacob(3,3)
      ! rotation matrix that transforms global coordinates
      ! into local coordinates of the midsurface frame of reference
      call k_matrix_zero(R,mcrd,mcrd) ! R(3,3)
      ! coordinates of the nodes in the local frame of reference
      call k_matrix_zero(coord_l,mcrd,nnode) ! coord_l(3,8)
      call k_vector_zero(GP_coord,2) ! surface Gauss point coordinates (2D)
      call k_vector_zero(sf,4) ! sf(4) = shape functions calculated at gauss points
      ! transformation matrices contain one R matrix for each node
      call k_matrix_zero(Transformation_M,ndofel,ndofel) ! Transformation_M(24,24)
      call k_matrix_zero(Transformation_M_T,ndofel,ndofel) ! Transformation_M_T(24,24)
! B(3,24) matrix connects global nodal displacement with local separation
      call k_matrix_zero(B,mcrd,ndofel)
      call k_matrix_zero(B_t,ndofel,mcrd)
      call k_matrix_zero(temp1,mcrd,ndofel) ! used in B_t * Trac_Jacob * B
      ! stiffness matrix expressed in the local frame of reference
      call k_matrix_zero(stiff_l,ndofel,ndofel) ! stiff_l(24,24)
      call k_matrix_zero(temp2,ndofel,ndofel) ! used in T' * K * T
      ! stiffness matrix expressed in the global frame of reference
      call k_matrix_zero(stiff_g,ndofel,ndofel) ! stiff_g(24,24)
      ! residual in the local and global coordinates
      call k_matrix_zero(residual_l,ndofel,nrhs) ! residual_l(24,1)
      call k_matrix_zero(residual_g,ndofel,nrhs) ! residual_g(24,1)
      ! rows contain vectors on the midsurface
      call k_matrix_zero(aJacob_M,2,3)
      ! right hand side of the overall system of equations
      ! output for Abaqus
      call k_matrix_zero(rhs,ndofel,nrhs) ! rhs(24,1)
! amatrx(24,24) is the stiffness matrix of cohesive element
! given by B^T R^T D_{loc} R B (Eq. 4 in Spring 2014)
! D_{loc} is the local tangent stiffness of the cohesive element
! R is the rotation matrix of nodal displacement
! local coordinates are determined by connecting the midside points of the cohesive elements
! R transforms the global coordinates into the local coordinates
      call k_matrix_zero(amatrx,ndofel,ndofel)
! co_de(3,8) are the deformed coordinates of the nodes
      call k_matrix_zero(co_de,mcrd,nnode) 
      a_Jacob=0.d0
	  
! The following will create an infinite loop,
! giving sufficient time to attach the debugger to 
! the running process once the program is running
       debug = 0
       DO WHILE(debug .eq. 1 .and. KINC == 1)
       !paused for debugging
       END DO
c
c Do local computations================================================
c Determine deformed coordinates
      do i = 1,mcrd
         do j = 1,nnode
            co_de(i,j)=coords(i,j)+U(3.0*(j-1.0)+i)
         end do
      end do

c
c Do Calculations at Gauss Points======================================
c
c Begin cycling over Gauss points
c
      do i = 1,GP_n
c
      call k_matrix_zero(aJacob_M,2,3)
c
      gpt = i ! Gauss point index

c Damage at Gauss point is stored into state variables
      damage(i) = SVARS(i)

c
c Define rotation matrix R that transforms global coordinates
c into local coordinates referred to the midplane
c
      call k_local_coordinates(gpt,co_de,R,coord_l,Transformation_M,
     & Transformation_M_T,a_Jacob,aJacob_M,coords,u,ndofel,nnode,
     & mcrd,SVARS,KINC)
c
c Compute shear and normal local opening displacements==================
c note that the order of the nodes is 1,2,3,4 on the bottom face
c corresponding to 5,6,7,8 on the top face
c Coordinates are in the midsurface reference frame
c 
      do j = 1,4
         ds1(j)=coord_l(1,j+4)-coord_l(1,j)
         ds2(j)=coord_l(2,j+4)-coord_l(2,j)
         dn(j) =coord_l(3,j+4)-coord_l(3,j)
      end do

c
c Determine the values of the shape function at Gauss Point i
c
         call k_shape_fun(i,sf)
c
         call k_vector_zero(delu_loc_gp,mcrd)
c
c Determine shear and normal opening displacements at Gauss points
c Coordinates are in the midsurface reference frame
c         
         do j = 1,4
            delu_loc_gp(1)=delu_loc_gp(1)+ds1(j)*sf(j)
            delu_loc_gp(2)=delu_loc_gp(2)+ds2(j)*sf(j)
            delu_loc_gp(3)=delu_loc_gp(3)+dn(j)*sf(j)
         end do
         
c
c Shear distance calculation
c
         opn=delu_loc_gp(3)
         opt=sqrt(delu_loc_gp(1)**2.0+delu_loc_gp(2)**2.0)

c
c Determine Traction vector and tangent modulus matrix
c Bilinear cohesive law
c Trac is a pressure
c Trac_Jacob is the derivative of Trac with respect to
c crack opening vector in the local frame of reference
c with coordinates in the midsurface reference frame

         call k_cohesive_law(Trac,Trac_Jacob,stiffK,stiffG,sigma0,
     & delta0,deltaF,delu_loc_gp,mcrd,nrhs,SVARS,damage(i),KINC,deltaortiz)

! assign new damage to state variables
! and to global variable
         SVARS(i) = damage(i)

      ! output effective opening in the bulk elements as the
      ! maximum on the four interface elements connected to
      ! that element
      call MutexLock( 15 )      ! lock Mutex #15
      do j = 1,nintpts
        if (deltaortiz > kDeltaEff(tempneighelem1,j)) then
          kDeltaEff(tempneighelem1,j) = deltaortiz
        end if
        if (deltaortiz > kDeltaEff(tempneighelem2,j)) then
          kDeltaEff(tempneighelem2,j) = deltaortiz
        end if
      end do
      call MutexUnlock( 15 )   ! unlock Mutex #15

c
c Determine B matrix and its transpose
c B connects the displacement at the nodes with
c the crack opening along the midsurface
c 
         call k_Bmatrix(sf,B,mcrd,ndofel)
c
         call k_matrix_transpose(B,B_t,mcrd,ndofel)
c
c Compute the stiffness matrix
c Local Stiffness = B_t * Trac_Jacob * B
c Calculated with coordinates in the midsurface reference frame
c
         call k_matrix_multiply(Trac_Jacob,B,temp1,mcrd,mcrd,
     & ndofel)
         call k_matrix_multiply(B_t,temp1,stiff_l,ndofel,
     & mcrd,ndofel)
c
c Compute Global stiffness matrix
c Global_K = T' * K * T
c Transformation_M rotates coordinates of each node
c from global to local system of the midsurface
c
         call k_matrix_multiply(Transformation_M_T,stiff_l,
     & temp2,ndofel,ndofel,ndofel)
         call k_matrix_multiply(temp2,Transformation_M,stiff_g,
     & ndofel,ndofel,ndofel)
c
c Multiply Jacobian with the Global stiffness and add contribution
c from each Gauss Point
c Equation 4 in Spring 2014
c a_Jacob is the scalar value of the area 
c of the midsurface (for 1 Gauss point)
c
         call k_matrix_plus_scalar(amatrx,stiff_g,a_Jacob,
     & ndofel,ndofel)
c
c Compute the global residual vector
c Local_residual = B_t * Trac
c Global_residual = T' * Local_residual
c Equation 5 in Spring 2014
c 
         call k_matrix_multiply(B_t,Trac,residual_l,ndofel,
     & mcrd,nrhs)
         call k_matrix_multiply(Transformation_M_T,residual_l,
     & residual_g,ndofel,ndofel,nrhs)
c
c Multiply the Global residual by the Jacobian and add the 
c contribution from each point
c   
         call k_matrix_plus_scalar(rhs,residual_g,a_Jacob,
     & ndofel,nrhs)

c
c End cycling over Gauss points
c
      end do
c
      return
      end
c======================================================================
c=============================SUBROUTINES==============================
c======================================================================
c
c Determine the strain-displacement (B) matrix
c B connects the displacement at the nodes with
c the crack opening along the midsurface
c
      subroutine k_Bmatrix(sf,B,mcrd,ndofel)
      INCLUDE 'ABA_PARAM.INC'
      dimension sf(4),B(mcrd,ndofel)
      B(1,1) =  sf(1)
      B(1,4) =  sf(2)
      B(1,7) =  sf(3)
      B(1,10)=  sf(4)
      B(1,13)= -sf(1)
      B(1,16)= -sf(2)
      B(1,19)= -sf(3)
      B(1,22)= -sf(4)
      B(2,2) =  sf(1)
      B(2,5) =  sf(2)
      B(2,8) =  sf(3)
      B(2,11)=  sf(4)
      B(2,14)= -sf(1)
      B(2,17)= -sf(2)
      B(2,20)= -sf(3)
      B(2,23)= -sf(4)
      B(3,3) =  sf(1)
      B(3,6) =  sf(2)
      B(3,9) =  sf(3)
      B(3,12)=  sf(4)
      B(3,15)= -sf(1)
      B(3,18)= -sf(2)
      B(3,21)= -sf(3)
      B(3,24)= -sf(4)
c
      return
      end
c======================================================================
c Bilinear cohesive law
c

      subroutine k_cohesive_law(T,T_d,stiffK,stiffG,sigma0,
     & delta0,deltaF,delu,mcrd,nrhs,SVARS,damage,KINC,deltaortiz)

      INCLUDE 'ABA_PARAM.INC'

      dimension T(mcrd,nrhs),T_d(mcrd,mcrd),delu(mcrd),SVARS(4)
       DOUBLE PRECISION stiffK, stiffG, sigma0, delta0, 
     & deltaF, popn, popt, Tn, Tt, damage, deltaD, maxS,
     & Dnn, Dnt, Dtt, Dtn

c weight factor in front of shear displacement
c determines how much shear contributes to damage
c corresponds to beta^2 in Ortiz, Pandolfi, 1999
      DOUBLE PRECISION smaxratio

c effective displacement delta which determines the damage
c see Ortiz, Pandolfi 1999
      DOUBLE PRECISION deltaortiz

c Derivatives of the damage with respect to popn and popt
c when damage is increase, otherwise zero
      DOUBLE PRECISION ddamagedpopn, ddamagedpopt
      DOUBLE PRECISION ddamagetemp

c
c Define normal and shear displacement
c
      popn=delu(3)
      popt=sqrt(delu(1)**2.0+delu(2)**2.0) ! always positive
      
! Damage accumulation
! equation holds for every case
! damage also in shear
! bilinear relationship determined by effective length deltaortiz

      smaxratio = 0.1
      if (popn .GT. 0.0) then
        deltaortiz = sqrt(smaxratio*popt*popt + popn*popn) ! effective displacement > 0
      else
	deltaortiz = popt*sqrt(smaxratio) ! shear/compression damage > 0
      end if

      ! deltaD is the value of deltaortiz
      ! above which the mechanical behaviour becomes softening
      ! calculated with the damage at the previous time step
      ! always > 0
      deltaD = delta0*deltaF/(deltaF-damage*(deltaF-delta0))

      if (damage .LT. 1.0) then
        if (deltaortiz .GT. deltaF) then ! fully damaged
          damage = 1.0
        else
          if (deltaortiz .GT. deltaD) then ! damage increase
            damage = deltaF*(deltaortiz-delta0)/(deltaortiz*(deltaF-delta0))
          end if ! otherwise constant
        end if
      else ! limit damage to 1.0
        damage = 1.0
      end if

c Determine load case and
c calculate the normal cohesive traction Tn
c and shear traction Tt
c
      if (damage .LE. 0.0) then ! undamaged

        ! uploading (popn > 0)
        ! compression (popn < 0)
        Tn = stiffK*popn
        Tt = stiffG*popt

        Dnn = stiffK
        Dtt = stiffG

        Dtn = 0.0 ! derivative of Tt with respect to popn
	Dnt = 0.0 ! derivative of Tn with respect to popt

      elseif (damage .GE. 1.0) then ! completely damaged

        if (popn .GT. 0.0) then ! crack open
          Tn = 0.0
          Dnn = 0.0
        else ! contact
          Tn = stiffK*popn
          Dnn = stiffK
        end if

	! shear traction always 0 when completely damaged
	! because of no friction
        Tt = 0.0
        Dtt = 0.0

        Dtn = 0.0 ! derivative of Tt with respect to popn
	Dnt = 0.0 ! derivative of Tn with respect to popt

      else ! partially damaged

        if (deltaortiz .GT. deltaD) then ! damage increase
          ddamagetemp = deltaF*delta0
          ddamagetemp = ddamagetemp/(deltaortiz*deltaortiz*deltaortiz)
          ddamagetemp = ddamagetemp/(deltaF-delta0)
          ddamagedpopn = ddamagetemp*popn
          ddamagedpopt = smaxratio*ddamagetemp*popt
        else ! damage is constant -> no derivatives of damage
          ddamagedpopn = 0.0
          ddamagedpopt = 0.0
        end if

        Tt = (1.0-damage)*stiffG*popt
        Dtt = (1.0-damage)*stiffG - stiffG*popt*ddamagedpopt
        Dtn = -stiffG*popt*ddamagedpopn

        if (popn .GT. 0.0) then
          Tn = (1.0-damage)*stiffK*popn
          Dnn = (1.0-damage)*stiffK - stiffK*popn*ddamagedpopn
          Dnt = -stiffK*popn*ddamagedpopt
        else ! contact
          Tn = stiffK*popn
          Dnn = stiffK
          Dnt = 0.0
        end if

      end if

      ! assign components of the traction vector
      T(3,1) = Tn

      if (popt .EQ. 0.0) then ! avoid division by zero
        T(1,1) = 0.0
        T(2,1) = 0.0
      else
        T(1,1) = Tt*delu(1)/popt
        T(2,1) = Tt*delu(2)/popt
      end if

      ! assign components of the derivatives of
      ! the traction vector

      T_d(3,3) = Dnn

      if (popt .EQ. 0.0) then
         T_d(1,1) = Dtt
         T_d(1,2) = 0.0d00
         T_d(1,3) = 0.0d00
         T_d(2,1) = 0.0d00
         T_d(2,2) = Dtt
         T_d(2,3) = 0.0d00
         T_d(3,1) = 0.0d00
         T_d(3,2) = 0.0d00
         
      else
         T_d(1,1)=Dtt*(delu(1)/popt)**2.0+Tt*((delu(2)**2.0)/(popt
     &    **3.0))
         T_d(1,2)=Dtt*delu(1)*delu(2)/(popt**2.0)
     &    -Tt*delu(1)*delu(2)/(popt**3.0)
	 T_d(1,3)=Dtn*delu(1)/popt
c
         T_d(2,1)=T_d(1,2)
         T_d(2,2)=Dtt*(delu(2)/popt)**2.0+Tt*((delu(1)**2.0)/(popt
     &    **3.0))
         T_d(2,3)=Dtn*delu(2)/popt

         T_d(3,1) = Dnt*delu(1)/popt
         T_d(3,2) = Dnt*delu(2)/popt
         
      end if

c
      return
      end

c=====================================================================
c determine the rotation matrix R that transforms
c global coordinates into local coordinates of the midsurface
c determine the Jacobian: a_Jacob, aJacob_M
c
      subroutine k_local_coordinates(gpt,co_de,R,coord_l,Transformation_M,
     & Transformation_M_T,a_Jacob,aJacob_M,coords,u,ndofel,nnode,
     & mcrd, SVARS, KINC)
      INCLUDE 'ABA_PARAM.INC'
      dimension R(mcrd,mcrd),coord_l(mcrd,nnode),aJacob_M(2,3),
     & Transformation_M(ndofel,ndofel),coords(mcrd,nnode),
     & Transformation_M_T(ndofel,ndofel),u(ndofel),
     & co_de(mcrd,nnode), co_de_m(3,4),SFD(2,4),SVARS(4)

! coordinate components of the 4 midpoints of the cohesive element
       DOUBLE PRECISION x1, x2, x3, x4, y1, y2, y3, y4, y5, y6, z1, z2,
     & z3, z4
c
! co_de_m(3,4) are the coordinates of the 4 midpoints of the cohesive element
      call k_matrix_zero(co_de_m,3,4)
c
      do i = 1,3
         co_de_m(i,1)=(co_de(i,1)+co_de(i,5))*0.5
         co_de_m(i,2)=(co_de(i,2)+co_de(i,6))*0.5
         co_de_m(i,3)=(co_de(i,3)+co_de(i,7))*0.5
         co_de_m(i,4)=(co_de(i,4)+co_de(i,8))*0.5
      end do
      
      !if (KINC == 1000) then
      !  write(*,*) 'co_de_m(1,2)'
      !  write(*,*) co_de_m(1,2)
      !  write(*,*) 'co_de_m'
      !  write(*,*) co_de_m
      !end if
c
      x1=co_de_m(1,1)
      x2=co_de_m(1,2)
      x3=co_de_m(1,3)
      x4=co_de_m(1,4)
c
      y1=co_de_m(2,1)
      y2=co_de_m(2,2)
      y3=co_de_m(2,3)
      y4=co_de_m(2,4)
c
      z1=co_de_m(3,1)
      z2=co_de_m(3,2)
      z3=co_de_m(3,3)
      z4=co_de_m(3,4)
c
c Coordinates of the Gauss points in the reference element
c with edges in plus/minus 1
c
      if (gpt .eq. 1) then
         c_r=-sqrt(1.0d0/3.0d0)
         c_s=-sqrt(1.0d0/3.0d0)
      elseif (gpt .eq. 2) then
         c_r= sqrt(1.0d0/3.0d0)
         c_s=-sqrt(1.0d0/3.0d0)
      elseif (gpt .eq. 3) then
         c_r= sqrt(1.0d0/3.0d0)
         c_s= sqrt(1.0d0/3.0d0)
      elseif (gpt .eq. 4) then
         c_r=-sqrt(1.0d0/3.0d0)
         c_s= sqrt(1.0d0/3.0d0)
      end if
c
c Shape function derivatives
c with respect to coordinates 1 and 2
c in the reference element = dphi/dxi
c
      SFD(1,1) =-0.25*(1-c_s)
      SFD(1,2) = 0.25*(1-c_s)
      SFD(1,3) = 0.25*(1+c_s)
      SFD(1,4) =-0.25*(1+c_s)
      SFD(2,1) =-0.25*(1-c_r)
      SFD(2,2) =-0.25*(1+c_r)
      SFD(2,3) = 0.25*(1+c_r)
      SFD(2,4) = 0.25*(1-c_r)
c
c Jacobian matrix: coordinates of the midpoints are multiplied
c by the derivatives dphi/dxi, since co_de_m * phi
c are the coordinates on the mid surface using shape function interpolation
c then aJacob_M = dx/dxi
c the two rows of aJacob_M are 3D vectors indicating the two edges
c of the mid surface
c
      do i = 1,2
         do j = 1,3
            do k =1,4
               aJacob_M(i,j) = aJacob_M(i,j) + SFD(i,k)*co_de_m(j,k)
            end do
         end do
      end do
c
c Vector product to find the normal of the mid surface
c
      dum1 = aJacob_M(1,2)*aJacob_M(2,3) - aJacob_M(1,3)*aJacob_M(2,2)
      dum2 = aJacob_M(1,3)*aJacob_M(2,1) - aJacob_M(1,1)*aJacob_M(2,3)
      dum3 = aJacob_M(1,1)*aJacob_M(2,2) - aJacob_M(1,2)*aJacob_M(2,1)
c
      a_Jacob = sqrt(dum1**2 + dum2**2 + dum3**2)
c
c Third row of the rotation matrix is the surface normal
c see Eq. 7 in Spring 2014
c
      R(3,1) = dum1/a_Jacob
      R(3,2) = dum2/a_Jacob
      R(3,3) = dum3/a_Jacob
c
c Definition of the first row of rotation matrix R
c as the first vector of aJacob_M, which is the vector on the mid surface
c determined by the coordinate c_s in the reference element
c
      aLen=sqrt(aJacob_M(1,1)**2.0d00+aJacob_M(1,2)**2.0d00+aJacob_M(1,3)**2.0d00)
      R(1,1)=aJacob_M(1,1)/aLen
      R(1,2)=aJacob_M(1,2)/aLen
      R(1,3)=aJacob_M(1,3)/aLen
c
c Definition of the second row of rotation matrix R
c as the vector product between midsurface normal and
c first vector of aJacob_M
c
      aLen2a=R(3,2)*R(1,3)-R(3,3)*R(1,2)
      aLen2b=R(3,3)*R(1,1)-R(3,1)*R(1,3)
      aLen2c=R(3,1)*R(1,2)-R(3,2)*R(1,1)
c
      aLen2 = sqrt(aLen2a**2.0d00 + aLen2b**2.0d00 + aLen2c**2.0d00)
c
      R(2,1) = aLen2a/aLen2
      R(2,2) = aLen2b/aLen2
      R(2,3) = aLen2c/aLen2
c
c Components of the Jacobian matrix aJacob_M
c are recalculated
c a_J11 = aJacob_M(1,1)
c a_J12 = aJacob_M(1,3)
c a_J21 = aJacob_M(2,1)
c a_J22 = aJacob_M(2,3)
c
      a_J11 = (-0.25*(1-c_s))*x1 + (0.25*(1-c_s))*x2 + (0.25*(1+c_s))*x3 +
     & (-0.25*(1+c_s))*x4
      a_J12 = (-0.25*(1-c_s))*z1 + (0.25*(1-c_s))*z2 + (0.25*(1+c_s))*z3 +
     & (-0.25*(1+c_s))*z4
      a_J21 = (-0.25*(1-c_r))*x1 + (-0.25*(1+c_r))*x2 + (0.25*(1+c_r))*x3 +
     & (0.25*(1-c_r))*x4
      a_J22 = (-0.25*(1-c_r))*z1 + (-0.25*(1+c_r))*z2 + (0.25*(1+c_r))*z3 +
     & (0.25*(1-c_r))*z4
c
      b_J11 = (-0.25*(1-c_s))*x1 + (0.25*(1-c_s))*x2 + (0.25*(1+c_s))*x3 +
     & (-0.25*(1+c_s))*x4
      b_J12 = (-0.25*(1-c_s))*y1 + (0.25*(1-c_s))*y2 + (0.25*(1+c_s))*y3 +
     & (-0.25*(1+c_s))*y4
      b_J21 = (-0.25*(1-c_r))*x1 + (-0.25*(1+c_r))*x2 + (0.25*(1+c_r))*x3 +
     & (0.25*(1-c_r))*x4
      b_J22 = (-0.25*(1-c_r))*y1 + (-0.25*(1+c_r))*y2 + (0.25*(1+c_r))*y3 +
     & (0.25*(1-c_r))*y4
c
      c_J11 = (-0.25*(1-c_s))*y1 + (0.25*(1-c_s))*y2 + (0.25*(1+c_s))*y3 +
     & (-0.25*(1+c_s))*y4
      c_J12 = (-0.25*(1-c_s))*z1 + (0.25*(1-c_s))*z2 + (0.25*(1+c_s))*z3 +
     & (-0.25*(1+c_s))*z4
      c_J21 = (-0.25*(1-c_r))*y1 + (-0.25*(1+c_r))*y2 + (0.25*(1+c_r))*y3 +
     & (0.25*(1-c_r))*y4
      c_J22 = (-0.25*(1-c_r))*z1 + (-0.25*(1+c_r))*z2 + (0.25*(1+c_r))*z3 +
     & (0.25*(1-c_r))*z4
c
c The following are the components of the normal vector
c to the mid surface, calculated using vector product of
c aJacob_M(1,1:3) and aJacob_M(2,1:3)
c
      a_Jacob1 = (a_J11*a_J22 - a_J12*a_J21)
      a_Jacob2 = (b_J11*b_J22 - b_J12*b_J21)
      a_Jacob3 = (c_J11*c_J22 - c_J12*c_J21)
c
c a_Jacob is the scalar value of the area 
c of the midsurface (for 1 Gauss point)
c note that this area change from Gauss point
c to Gauss point
c
      a_Jacob = sqrt(a_Jacob1**2.0d00 + a_Jacob2**2.0d00 + a_Jacob3**2.0d00)
c
c======================================================================
      num=nnode
c
c Transformation_M(24,24) contains 8 times the rotation matrix
c on the diagonal 3x3 regions of the matrix
c It can rotate the coordinates of each node independently
c
      do i = 1,num
         dum=3.0*(i-1.0)
         Transformation_M(dum+1,dum+1)=R(1,1)
         Transformation_M(dum+1,dum+2)=R(1,2) 
         Transformation_M(dum+1,dum+3)=R(1,3)
         Transformation_M(dum+2,dum+1)=R(2,1)
         Transformation_M(dum+2,dum+2)=R(2,2)
         Transformation_M(dum+2,dum+3)=R(2,3)
         Transformation_M(dum+3,dum+1)=R(3,1)
         Transformation_M(dum+3,dum+2)=R(3,2)
         Transformation_M(dum+3,dum+3)=R(3,3)
      end do
c
      call k_matrix_transpose(Transformation_M,Transformation_M_T,
     $ ndofel,ndofel)
c
c coord_l(3,8) are the coordinates of the cohesive element nodes
c expressed in the local coordinate system
c which is the coordinate system of the midsurface
c
      do i = 1, nnode
         coord_l(1,i)=(R(1,1)*co_de(1,i)+R(1,2)*co_de(2,i)
     & +R(1,3)*co_de(3,i))
         coord_l(2,i)=(R(2,1)*co_de(1,i)+R(2,2)*co_de(2,i)
     & +R(2,3)*co_de(3,i))
         coord_l(3,i)=(R(3,1)*co_de(1,i)+R(3,2)*co_de(2,i)
     & +R(3,3)*co_de(3,i))
      end do
c
      !if (KINC == 2) then
      !    write(*,*) 'R(2,1)'
      !    write(*,*) R(2,1)
      !    write(*,*) 'R'
      !    write(*,*) R
      !end if
c    
      return
      end
c======================================================================
c Calculate shape function at gauss point i
c of the midsurface element
c
      subroutine k_shape_fun(i,sf)
      INCLUDE 'ABA_PARAM.INC'
      dimension sf(4), GP_coord(2)
c
      if (i .eq. 1) then
         GP_coord(1)=-sqrt(1.0d0/3.0d0)
         GP_coord(2)=-sqrt(1.0d0/3.0d0)
      elseif (i .eq. 2) then
         GP_coord(1)= sqrt(1.0d0/3.0d0)
         GP_coord(2)=-sqrt(1.0d0/3.0d0)
      elseif (i .eq. 3) then
         GP_coord(1)= sqrt(1.0d0/3.0d0)
         GP_coord(2)= sqrt(1.0d0/3.0d0)
      elseif (i .eq. 4) then
         GP_coord(1)=-sqrt(1.0d0/3.0d0)
         GP_coord(2)= sqrt(1.0d0/3.0d0)
      end if
c
      sf(1)=(1-GP_coord(1))*(1-GP_coord(2))*0.25
      sf(2)=(1+GP_coord(1))*(1-GP_coord(2))*0.25
      sf(3)=(1+GP_coord(1))*(1+GP_coord(2))*0.25
      sf(4)=(1-GP_coord(1))*(1+GP_coord(2))*0.25
c
      return
      end
c======================================================================
c Multiply matrices A and B to give C
c
      subroutine k_matrix_multiply(A,B,C,l,n,m)
      INCLUDE 'ABA_PARAM.INC'
      dimension A(l,n),B(n,m),C(l,m)
c
      call k_matrix_zero(C,l,m)
c
      do i = 1, l
         do j = 1, m
            do k = 1, n
               C(i,j)=C(i,j)+A(i,k)*B(k,j)
            end do
         end do
      end do
c
      return
      end
c======================================================================
c Sum matrices with prefactor
c
      subroutine k_matrix_plus_scalar(A,B,c,n,m)
      INCLUDE 'ABA_PARAM.INC'
      dimension A(n,m),B(n,m)
c
      do i = 1, n
         do j = 1, m
            A(i,j)=A(i,j)+c*B(i,j)
         end do
      end do
c
      return
      end
c======================================================================
c Transpose matrix
      subroutine k_matrix_transpose(A,B,n,m)
      INCLUDE 'ABA_PARAM.INC'
      dimension A(n,m),B(m,n)
c
      do i = 1,n
         do j = 1,m
            B(j,i)=A(i,j)
         end do
      end do
c
      return
      end
c======================================================================
c Initialize all the components of a matrix to zero
      subroutine k_matrix_zero(A,n,m)
      INCLUDE 'ABA_PARAM.INC'
      dimension A(n,m)
c
      do i = 1, n
         do j = 1, m
            A(i,j)=0.d0
         end do
      end do
c
      return
      end
c======================================================================
c Initialize all the components of a vector to zero
      subroutine k_vector_zero(A,n)
      INCLUDE 'ABA_PARAM.INC'
      dimension A(n)
c
      do i = 1, n
         A(i)=0.d0
      end do
c
      return
      end
c======================================================================
c Calculate Macaulay brackets
      subroutine k_Mac(pM,a,b)
      INCLUDE 'ABA_PARAM.INC'
c
      if ((a-b) .GE. 0.0) then
         pM=a-b
      elseif ((a-b) .LT. 0.0) then
         pM=0.d0
      end if
c
      return
      end
c======================================================================
c=================================END==================================
c======================================================================









      
