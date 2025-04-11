!
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
!     Subroutine for initialization 
      use initializations, only : initialize_variables,
     + initialize_once 
      use userinputs, only : gndmodel, backstressmodel
      use globalvariables, only: numel, numpt, 
     + init_once, statev_gmatinv, statev_gmatinv_t,
     + statev_gammasum_t, statev_gammasum,
     + statev_gammadot_t, statev_gammadot,
     + statev_Fp_t, statev_Fp, statev_sigma_t, statev_sigma,
     + statev_jacobi_t, statev_jacobi, statev_Fth_t, statev_Fth,
     + statev_tauc_t, statev_tauc, statev_maxx_t, statev_maxx,
     + statev_Eec_t, statev_Eec, statev_gnd_t, statev_gnd,
     + statev_ssd_t, statev_ssd, statev_forest_t, statev_forest,
     + statev_substructure_t, statev_substructure, statev_evmp_t,
     + statev_evmp, statev_totgammasum_t, statev_totgammasum,
     + statev_ssdtot_t, statev_ssdtot, statev_loop_t, statev_loop,
     + statev_Lambda_t, statev_Lambda, statev_sigma_t2,
     + statev_tausolute_t, statev_tausolute, time_old, dt_t,
     + statev_backstress, statev_backstress_t,
     + statev_plasdiss, statev_plasdiss_t, 
     + ip_count, calculategradient, ip_init, grad_init
!
      use straingradients, only: gndmodel1, gndmodel2, 
     + gndmodel3, gndmodel4
      use backstress, only: backstressmodel2
      use meshprop, only: initialize_gradientoperators
      use useroutputs, only: write_statev_legend   
!
      implicit none
!
!
!
      INTEGER,                        INTENT(IN) ::
     + LOP,
     + LRESTART,
     + KSTEP,
     + KINC
      REAL(8), DIMENSION(1),          INTENT(IN) ::
     + DTIME
      REAL(8), DIMENSION(2),          INTENT(IN) ::
     + TIME
!
!
      integer :: i
!
!
!
!     at the start of the analysis (only ONCE!)
!     if there are initializations/calculations that are needed once,
!     and that are independepent of element properties, you may use here.
      if (LOP==0) then
!
!         0. Initialize variables (instead of using data statements)
          call initialize_variables
!
!
      end if
!
!     Initialization of gradient operators
!     Check if the element has more than one integration point
!     And the gradient initialization was done before or not
      if ((calculategradient==1).and.(grad_init==0)) then
!
!         Check if IP coordinates were collected
          if ((sum(ip_init)==numel*numpt).and.
     + (numel>0).and.(numpt>0)) then
!
!
!             Calculate the gradient mapping for each element once.
              call initialize_gradientoperators
!
              grad_init=1
              write(*,*) '7. Gradient operators initialized!'
!
          endif
!
      end if
!
!
!
!
!
      if (LOP==1) then
!
!
!         in case of force BC
          if ((KINC==1).and.(KSTEP==1)) then
!
!             
!             
!             check if the one-time initialization is done or not
              if ((init_once==0).and.(numel>0).and.(numpt>0)) then
!            
!
!
!                 check the number of elements
                  i = sum(ip_count(1,:))
                  if (i>numpt) numpt = i
!
!                 check the number of elements
                  i = sum(ip_count(:,1))
                  if (i>numel) numel = i
!
!                 
!                 write element number
                  write(*,*) 'total elements in the mesh: ', numel
!
!                 write integration point number
                  write(*,*) 'integration points per element: ', numpt
!
!
                  write(*,*) '--------------------------------'
!
!                 call the initializations
                  call initialize_once
!
!                 display upon completion
                  write(*,*) 'one-time initialization is complete!'
!
!
                  write(*,*) '--------------------------------'
!
!
!                 write the legend to a text file
                  call write_statev_legend
!
!                 message for initializatoin
                  write(*,*) '6. "STATEV_legend.txt" file is ready!'
!
!                 set the one-time initilaziation flag (at the very end)
                  init_once=1
!
              end if
!
!
!
!
!
          end if
!
!
      end if
!
!
!
!     at the end of each increment
      if (LOP==2) then
!
!
!         in case of displacement BC
          if ((KINC==1).and.(KSTEP==1)) then
!
!             
!             
!             check if the one-time initialization is done or not
              if ((init_once==0).and.(numel>0).and.(numpt>0)) then
!            
!
!
!                 check the number of elements
                  i = sum(ip_count(1,:))
                  if (i>numpt) numpt = i
!
!                 check the number of elements
                  i = sum(ip_count(:,1))
                  if (i>numel) numel = i
!
!                 
!                 write element number
                  write(*,*) 'total elements in the mesh: ', numel
!
!                 write integration point number
                  write(*,*) 'integration points per element: ', numpt
!
!
                  write(*,*) '--------------------------------'
!
!                 call the initializations
                  call initialize_once
!
!                 display upon completion
                  write(*,*) 'one-time initialization is complete!'
!
!
                  write(*,*) '--------------------------------'
!
!
!                 write the legend to a text file
                  call write_statev_legend
!
!                 message for initializatoin
                  write(*,*) '6. "STATEV_legend.txt" file is ready!'
!
!
!
!                 set the one-time initilaziation flag (at the very end)
                  init_once=1
!
!             
              end if
!
!
!
!
!
!
          end if

!
!
!
          write(*,*) '--------------------------------'
!
          write(*,*) 'At the end of increment: ', KINC
!
!   
!     
!
!
!         GND calculations
!         do this only if the initiliazation complete and 
!         element gradients are calculatable (calculategradient==1)
          if ((init_once==1).and.(grad_init==1).and.
     + (calculategradient==1).and.(KINC>1)) then
!
!             update and calculate GNDs (nonlocal calculations)
!             this is done at the end of calculations.
!             GNDs that belong to the PREVIOUS time step are used.
!             initially GNDs are assumed to have "0" values.
!
!             calculate GNDs (if initialization complete)
              if (gndmodel>0) then
!
!                 curl followed by L2 minimization - SVD -
!                 constrained to active slip systems
                  if (gndmodel==1) then
!
                      call gndmodel1
!
!                 Cermelli-Gurtin's incompatibility measure - L2 minimization - SVD
!                 constrained to active slip systems (same as model-1)
                  elseif (gndmodel==2) then
!
                      call gndmodel2
!
!                 rate form followed by direct projections
                  elseif (gndmodel==3) then
!
                      call gndmodel3(DTIME(1))
!
!                 slip gradients
                  elseif (gndmodel==4) then
!
                      call gndmodel4(DTIME(1))
!
!
                  endif
!
!
                  write(*,*) 'GNDs are computed!'
!
!
                  if (backstressmodel==2) then
!
                      call backstressmodel2
!
                      write(*,*) 'Backstresses are computed!'
!
                  end if
!
              end if
!
          end if
!
!         update the time
          time_old = time(2)
!         store former converged time increment
          dt_t = DTIME(1)
!
!
!         update state variables
!         assign the current results to the former values
          statev_gmatinv_t(:,:,:,:) = statev_gmatinv(:,:,:,:)
          statev_evmp_t(:,:) = statev_evmp(:,:)
          statev_gammasum_t(:,:,:) = statev_gammasum(:,:,:)
          statev_totgammasum_t(:,:) = statev_totgammasum(:,:)
          statev_gammadot_t(:,:,:) = statev_gammadot(:,:,:)
          statev_Fp_t(:,:,:,:) = statev_Fp(:,:,:,:)
          statev_Fth_t(:,:,:,:) = statev_Fth(:,:,:,:)
          statev_sigma_t2(:,:,:) = statev_sigma_t(:,:,:)
          statev_sigma_t(:,:,:) = statev_sigma(:,:,:)
          statev_jacobi_t(:,:,:,:) = statev_jacobi(:,:,:,:)
          statev_tauc_t(:,:,:) = statev_tauc(:,:,:)
          statev_maxx_t(:,:) = statev_maxx(:,:)
          statev_Eec_t(:,:,:) = statev_Eec(:,:,:)
          statev_gnd_t(:,:,:) = statev_gnd(:,:,:)
          statev_ssd_t(:,:,:) = statev_ssd(:,:,:)
          statev_ssdtot_t(:,:) = statev_ssdtot(:,:)
          statev_forest_t(:,:,:) = statev_forest(:,:,:)
          statev_loop_t(:,:,:) = statev_loop(:,:,:)
          statev_substructure_t(:,:) = statev_substructure(:,:)
          statev_tausolute_t(:,:) = statev_tausolute(:,:)
          statev_Lambda_t(:,:,:) = statev_Lambda(:,:,:)
          statev_backstress_t(:,:,:) = statev_backstress(:,:,:)
          statev_plasdiss_t(:,:) = statev_plasdiss(:,:)
!
          write(*,*) 'States are updated!'
!
          write(*,*) '--------------------------------'
!
!
      endif
!
!
!
      RETURN
      END
!
!