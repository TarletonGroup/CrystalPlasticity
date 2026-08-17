C ***Subroutine module for modelling dislocation source controlled CP***
C Chris Hardie 2026    
      module sources
      implicit none
      contains
      
      !      
!     *******************************************************
!     *     This routine calculates source strengths        *
!     *     based on statistical distribution               *
!     *******************************************************      
      subroutine sourcestress(noel,npt, nslip)!, celent)
       use globalvariables, only: statev_source_stress, numpt, ipdomain, 
     +    statev_tauc, sourceparam_all, materialid
       use meshprop, only: initialize_gradientoperators
!     INPUTS
!     Element number
      integer, intent(in) :: noel
!     Integration point
      integer, intent(in) :: npt
!     Number of slip systems
      integer, intent(in) :: nslip
!     Source parameters
      real(8) :: refstress, refvolume, weibullmodulus
!     Random numbers
      real(8) :: r(nslip), ip_vol
!     Iteration integer
      integer :: is
!     Debug tools
      integer :: debug, debugwait, t
      
      debug=0
      do while (debug==1)
          debugwait = 1
      end do 
      
      refstress = sourceparam_all(materialid(noel,npt), 1)
      refvolume = sourceparam_all(materialid(noel,npt), 2)
      weibullmodulus = sourceparam_all(materialid(noel,npt), 3)
          
          CALL RANDOM_NUMBER(r) ! generate nslip random numbers
         DO is = 1, nslip 
!             Calculate the gradient mapping and effective volume for each integration point
             ip_vol=ipdomain(noel,npt)
             statev_source_stress(noel,npt,is)=refstress*(-(refvolume/ip_vol)*log(1.0-r(is)))**(1.0/weibullmodulus ) 
         END DO
         
      end subroutine sourcestress
!      
!     *******************************************************
!     *     This routine determines whether a source has    *
!     *     activated at this gauss point                   *
!     *******************************************************      
      subroutine sourceactivation(kinc, sigma, nslip, Schmidvec, tauceff_t, 
     +                        matid,noel,npt,dt,nors_t, COORDS, pnewdt)
      
      use globalvariables, only : statev_source_stress, 
     +    statev_source_list, statev_source_activation, numdim, 
     +    ipcoords, ipnors, statev_max_overstress, statev_source_list_2,
     +    featureid, statev_source_overstress
      
      use userinputs, only : tsource
      
!     INPUTS
      
!     Increment Number
      integer, intent(in) :: kinc
!     Material ID (grain number)
      integer, intent(in) :: matid
!     Number of slip systems
      integer, intent(in) :: nslip
!     Element number
      integer, intent(in) :: noel
!     Integration point number
      integer, intent(in) :: npt
!     Time increment
      real(8), intent(in) :: dt
!     Stress Components
      real(8), intent(in) :: sigma(6)
!     Array of Schmid Tensor components for each slip system
      real(8), intent(in)  :: Schmidvec(nslip,6)
!     CRSS
      real(8), intent(in) :: tauceff_t(nslip)
!     Rotated slip system normals
      real(8), intent(in) :: nors_t(nslip,numdim)
!     Integration point coordinates
      real(8), intent(in) :: COORDS(numdim)
      
!     OUTPUTS
      real(8), intent(out) :: pnewdt

!     Local Variables
!     Vector of source strengths for each slip system
      real(8) :: source_stress(nslip)
!     Maximum source overstress
      real(8) :: source_overstress
!     Shear stresses on each slip system
      real(8) :: tau(nslip), abstau(nslip), signtau(nslip)
!     Temporary variable for overstress
      real(8) :: temp
!     Slip system number for highest stressed source
      integer :: active_SS, sourcesign, grainid
      
      integer :: is, debug, debugwait
      
!     Initialise source list
      statev_source_list_2(noel, npt, :) = 0
      
!     Record up-to-date data
      
      ipcoords(noel,npt,1:numdim) = COORDS(1:numdim) ! ip coordinates
      ipnors(noel,npt,1:nslip,1:numdim) = nors_t(:,:)
      ! slip plane normals
      
      grainid=featureid(noel,npt)
           
      source_stress=statev_source_stress(noel,npt,:)
      source_overstress=0.
      active_SS=0
!         CALCULATE RESOLVED SHEAR STRESS ON SLIP SYSTEMS
!         rss and its sign
          do is = 1, nslip
              tau(is) = dot_product(Schmidvec(is,:),sigma)
              signtau(is) = sign(1.0,tau(is))
              abstau(is) = abs(tau(is))
              temp=abstau(is)-(source_stress(is)) !***Fraction of source stress for activation***
              if (temp .GT. source_overstress) then
                  source_overstress=temp
                  active_SS=is
                  sourcesign=signtau(is)
              end if
          
          end do  
  
      debug=0
      do while (debug==1)
          debugwait = 1
      end do
      
  
      if (active_SS .GT. 0) then
          IF (((source_overstress .LT. 0.01*source_stress(active_SS)) 
     +     .OR. (source_overstress .LT. 1.0)) ! Threshold 1 MPa
     +     .OR. (dt .LE. tsource)) THEN
            
          IF (source_overstress .GT. statev_max_overstress(grainid))THEN
                  statev_max_overstress(grainid)=source_overstress
          END IF
                              
              !record activation event
            !call MutexLock(2) ! lock Mutex 2
              !statev_source_activation(grainid)=statev_source_activation(grainid)+1
              !write(*,*) 'Source activation no.: ', statev_source_activation(grainid)
              statev_source_list_2(noel, npt, :) = (/kinc, grainid,statev_source_activation(grainid), active_SS, sourcesign/)
              statev_source_overstress(noel, npt) = source_overstress
            !call MutexUnlock(2)   ! unlock Mutex #2
            
            ELSE !IF (source_overstress .GT. source_stress(active_SS))THEN
                pnewdt=0.5 ! Half time increment
            END IF
      end if
           
      end subroutine sourceactivation
      
!     *******************************************************
!     *     This routine searches for local activated       *
!     *     dislocation sources and assigns relevant        *
!     *     state variables                                 *
!     *     (i.e. mobile dislocation density)               *
!     *******************************************************      
      subroutine sourceseek(kinc, matid, noel, npt, nslip, taus, active, Schmidvec, sigma_t)
      
      
      use globalvariables, only : statev_source_count, 
     +     statev_source_register, ipcoords, ipnors, numdim, ipdomain, 
     +     featureid, statev_source_stress, statev_active_flag,
     +     statev_source_newregister, numpt, numel
      
      use userinputs, only : tsource 
!     INPUTS
      
!     Increment number
      integer, intent(in) :: kinc
      
!     Material (phase) ID
      integer, intent(in) :: matid
      
!     Element number
      integer, intent(in) :: noel

!     Integration point number
      integer, intent(in) :: npt
      
!     Number of slip systems
      integer, intent(in) :: nslip
      
!     Schmid tensors in Voigt notation
      real(8), intent(in) :: Schmidvec(nslip, 6)
      
!     Stress at the start of the increment (previous stress)
      real(8), intent(in) :: sigma_t(6)
      
!     lattice friction/source stress
      real(8), intent(inout) :: taus(nslip)
      
!     OUTPUTS
!     Active flag per slip system
      real(8), intent(out) :: active(nslip)
!      
!!     Sign of shear stress for sources
!      integer, intent(out) :: source_sign(nslip)
      
!     LOCAL VARIABLES

!     Coordinates of Source
      real(8) :: sourceCoords(numdim)  
      
!     Coordinates of local position
      real(8) :: pointCoords(numdim)   
      
!     Vector between local point and source
      real(8) :: relativeCoords(numdim) 
      
!     Slip plane normal for source
      real(8) :: sourceNormal(numdim)
      
!     Absolute distance between point and source
      real(8) :: SP_distance
      
!     Source plane capture radius
      real(8) :: SP_bandwidth
      
!     Shear stress
      real(8) :: tau_t(nslip)
!     Absolute & sign of rss at the end of the previous increment
      real(8) :: abstau_t(nslip), signtau_t(nslip)
            
      real(8) :: test
      
!     source stress to implement as CRSS
      real(8) :: source_crss
      
!     Source activation increment
      integer :: source_inc
      
!     Element number for Source
      integer :: source_noel
      
!     Gauss point number for Source
      integer :: source_npt
      
!     Slip system number for Source
      integer :: source_slip
      
!     Sign of RSS of source
      integer :: source_sign
      
      integer :: J, is, debug, debugwait
            
      test=ipdomain(noel,npt)
      SP_bandwidth=2.0*ipdomain(noel,npt)**(1.0/3.0)
      do is = 1, nslip
          tau_t(is) = dot_product(Schmidvec(is,:),sigma_t)
          signtau_t(is) = sign(1.0,tau_t(is))
          abstau_t(is) = abs(tau_t(is))
      end do
      
      J=0
      ! Look for local activated sources elsewhere in grain
      DO WHILE (J .GE. 0)
          
             J=J+1 
             if ((statev_source_newregister(featureid(noel,npt),J,1)
     +        .NE. 0) .AND. (J .LE. numpt*numel)) then
              
              source_inc = statev_source_newregister(featureid(noel,npt),J,1)
              source_noel = statev_source_newregister(featureid(noel,npt),J,2)
              source_npt = statev_source_newregister(featureid(noel,npt),J,3)
              source_slip = statev_source_newregister(featureid(noel,npt),J,4)
              source_sign = statev_source_newregister(featureid(noel,npt),J,5)
                            
              sourceCoords(1:numdim) = ipcoords(source_noel,source_npt,1:numdim)
              pointCoords = ipcoords(noel,npt,1:numdim)
              
              source_crss=statev_source_stress(source_noel,source_npt,source_slip)
              
              relativeCoords=pointCoords-sourceCoords
              
              sourceNormal(1:numdim)=ipnors(source_noel,source_npt,source_slip,:)
              
              SP_distance=abs(dot_product(relativeCoords, sourceNormal))
              
              IF(SP_distance .LT. SP_bandwidth) THEN !If this integration point is close enough to an activated slip plane
                  IF (signtau_t(source_slip) == 
     +                sign(1,source_sign)) THEN
                      
                      if ((source_crss .LT. taus(source_slip)) .OR.        ! If the source is weaker than anything previously
     +       (statev_active_flag(noel, npt, source_slip) .EQ. 0.0)) then   ! Or nothing has facilitated slip on this system & location before
                      taus(source_slip)=source_crss ! take the lowest source stress
                      active(source_slip)=1.0 
                      statev_active_flag(noel, npt, source_slip)=1.0!statev_active_flag(noel, npt, source_slip)+1.0
                      end if
                     
                  END IF                  
              END IF
          else
          J=-1
          end if
          
      END DO
      
      
      end subroutine sourceseek
       
      
      end module sources