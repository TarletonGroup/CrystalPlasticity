!     Oct. 03rd, 2022
!     Eralp Demir
!
      module hardening
      implicit none
      contains
!     Since crss is computed considering the phase
!     Hardening shall evolve the proper state variables
!     
!
!     ********************************************
!     ** HARDENING updates the state variables  **
!     **  that determine mechanical hardening   **
!     ********************************************
      subroutine hardeningrules(iphase,nslip,
     + temperature,dt,G12,
     + burgerv,totgammasum,gammadot,pdot,
     + irradiationmodel,irradiationparam,
     + hardeningmodel,hardeningparam,
     + hintmat1,hintmat2, 
     + tauc_t,ssd_t,
     + loop_t, rhofor_t,rhosub_t,
     + tausolute,dtauc,
     + drhotot,drhofor,
     + drhosub, dssd, dloop)
      use globalvariables, only : KB
      use userinputs, only: maxnparam, maxnloop
      implicit none
!
!     INPUTS
!     crystal type
      integer,intent(in) :: iphase
!
!     number of slip systems
      integer,intent(in) :: nslip
!
!     current temperature
      real(8),intent(in) :: temperature
!
!     time increment
      real(8),intent(in) :: dt
!
!     shear modulus
      real(8),intent(in) :: G12
!
!     Burger's vector
      real(8), intent(in) :: burgerv(nslip)
!
!
!     cumulative crystallographic slip at the current tiemstep
!     needed for model with irradiation
      real(8),intent(in) :: totgammasum
!
!     plastic shear rate on slip systems
!     and absolute value
      real(8),intent(in) :: gammadot(nslip)
!
!     Von-mises equivalent plastic strain rate
      real(8),intent(in) :: pdot
!
!     irradiation effect
      integer,intent(in) :: irradiationmodel
!
!     irradiation model parameters
      real(8),intent(in) :: irradiationparam(maxnparam)
!
!     hardening model
      integer,intent(in) :: hardeningmodel
!
!     hardening parameters
      real(8),intent(in) :: hardeningparam(maxnparam)
!
!     Hardening interaction matrices
!     Latent hardening
      real(8), intent(in) :: hintmat1(nslip,nslip)
!     Hardening interaction matrix between dislocations
      real(8), intent(in) :: hintmat2(nslip,nslip)
!
!     crss
      real(8),intent(in) :: tauc_t(nslip)
!
!     ssd density
      real(8),intent(in) :: ssd_t(nslip)
!
!     loop density
      real(8),intent(in) :: loop_t(maxnloop) 
!
!     forest density
      real(8),intent(in) :: rhofor_t(nslip)
      
!     substructure density
      real(8),intent(in) :: rhosub_t
      
!
!
!     OUTPUTS
!     increase in tauc due to solute force
      real(8),intent(out) :: tausolute      
      
!     crss increment
      real(8),intent(out) :: dtauc(nslip)
!

!
!
!     total ssd density increment
      real(8),intent(out) :: drhotot
!
!     forest dislocation density increment
      real(8),intent(out) :: drhofor(nslip)
!
!     substructure dislocation density increment
      real(8),intent(out) :: drhosub
!
!     ssd density increment
      real(8),intent(out) :: dssd(nslip)
!
!     loop density increment
      real(8),intent(out) :: dloop(maxnloop)
!
!     Variables used within this subroutine
!     absolute value of the plastic shear rate on slip systems
      real(8), dimension(nslip) :: absgammadot
!     Parameters for irradiation hardening
      real(8) :: taus0, gammasat, gamma
!     Parameters for irradiation model = 2
      integer :: nloop
!     Parameters for Voce typ hardening model
      real(8) :: h0, ss, m, q, hb(nslip), Hab(nslip,nslip)
!     Parameter for linear hardening model
      real(8) :: k
!     Parameters for Kocks-Mecking hardening model
      real(8) :: k1, b, X, g, D, pdot0, k2, KBT
!     Additional parameters for substructure hardening
      real(8) f, ksub
!
      integer :: is, i, j
      integer :: il
!
!      
!     Set outputs zero initially
      dtauc = 0.
      dssd = 0.
      drhofor = 0.
      drhosub = 0.
      drhotot = 0.
      tausolute = 0.
      dloop = 0.
      
!     absolute value of slip rates
      absgammadot = dabs(gammadot)
!

!
!     Irradiation hardening      
!     =========================================================================      

       
!     update solute force
      if (irradiationmodel == 1) then
!
!         Prefactor in irradiation hardening
          taus0 = irradiationparam(1)
!
!         Saturation strain
          gammasat = irradiationparam(2)
!
!         Solute strength
          tausolute = taus0*dexp(-totgammasum/gammasat)
!
      end if
      
      
      
!     update loop density
      if (irradiationmodel == 2) then
!
!         Number of type of the loop defects
          nloop = int(irradiationparam(1))
!
!         Compute the evolution (annihiliation)
          do il = 1, nloop


              do is = 1, nslip

                  dloop(il) = dloop(il) - 0.5 * hintmat2(il,is) *
     + sqrt(loop_t(il)) / burgerv(is) * absgammadot(is) * dt

              end do


          end do


!
      end if      
      
      

!     =========================================================================



!     Other strainhardening models
!     No hardening
      if (hardeningmodel==0) then
          
!         Do not harden!

          
      elseif (hardeningmodel==1) then
          
!         read in the parameters
          h0 = hardeningparam(1)
          ss = hardeningparam(2)
          m = hardeningparam(3)
          q = hardeningparam(4)
          

          
!         Construct the latent hardening matrix
!         Note that this is a matrix different than latent hardening
          Hab = hintmat1
          
          hb=0.   
          do is = 1,nslip
              
              hb(is) = h0*(1.-tauc_t(is)/ss)**m*
     + absgammadot(is)*dt
              
              
          end do
          
          dtauc = matmul(Hab,hb)
          
          
!     linear hardening model        
      elseif (hardeningmodel==2) then
          
!         read in the parameters
          k = hardeningparam(1)
          
          
!         total density
          drhotot = k*pdot*dt
          
!         SSD evolution
          do is = 1, nslip
              dssd(is) = k*absgammadot(is)*dt
          end do
          
          
          
!     Kocks-Mecking hardening
      elseif (hardeningmodel==3) then
          
!         input parameters
          k1 = hardeningparam(1)
          X = hardeningparam(2)
          g = hardeningparam(3)
          D = hardeningparam(4)
          pdot0 = hardeningparam(5)

          
          KBT = KB*temperature
          


          
          do is = 1,nslip
              
!             Burgers vector
              b = burgerv(is)
              
              k2 = k1*X*b/g*(1.-KBT/D/b**3*dlog(pdot/pdot0))
              
              dssd(is) = (k1*dsqrt(rhofor_t(is))-k2*ssd_t(is))*
     + absgammadot(is)*dt
              
          end do
          
          drhotot = sum(dssd)
          
          
!     Kocks-Mecking hardening with substructure evolution
      elseif (hardeningmodel==4) then    
          
          
     
          
          
          
          
          
          

!         for alpha-uranium material parameters are built in
          if (iphase == 4) then



!             forest dislocations evolution
!             using constants from calibration of tensile bar 3
!             using twin-slip interaction model
              drhofor(1) = 43.2*max(dsqrt(rhofor_t(1))-
     + (0.17100+2.6093e-03*temperature)
     + *rhofor_t(1),0.0)*absgammadot(1)*dt
              drhofor(2) = 6320.0*max(dsqrt(rhofor_t(2))-
     + (0.25650+5.8708e-04*temperature)
     + *rhofor_t(2),0.0)*absgammadot(2)*dt
              drhofor(3) = 0.24*max(dsqrt(rhofor_t(3))-
     + (0.11718+1.9289e-04*temperature)
     + *rhofor_t(3),0.0)*absgammadot(3)*dt
              drhofor(4) = 0.24*max(dsqrt(rhofor_t(4))-
     + (0.11718+1.9289e-04*temperature)
     + *rhofor_t(4),0.0)*absgammadot(4)*dt
              drhofor(5) = 800.0*max(dsqrt(rhofor_t(5))-
     + (0.12+1.5e-05*temperature)
     + *rhofor_t(5),0.0)*absgammadot(5)*dt
              drhofor(6) = 800.0*max(dsqrt(rhofor_t(6))-
     + (0.12+1.5e-05*temperature)
     + *rhofor_t(6),0.0)*absgammadot(6)*dt
              drhofor(7) = 800.0*max(dsqrt(rhofor_t(7))-
     + (0.12+1.5e-05*temperature)
     + *rhofor_t(7),0.0)*absgammadot(7)*dt
              drhofor(8) = 800.0*max(dsqrt(rhofor_t(8))-
     + (0.12+1.5e-05*temperature)
     + *rhofor_t(8),0.0)*absgammadot(8)*dt

!             substructure dislocations evolution
              drhosub =  0.216*(17.545+0.26771*temperature)*
     + rhofor_t(1)*dsqrt(rhosub_t)*absgammadot(1)*dt
              
!         If any other material
          else
              
              
!             input parameters
              k1 = hardeningparam(1)
              X = hardeningparam(2)
              g = hardeningparam(3)
              D = hardeningparam(4)
              pdot0 = hardeningparam(5)
              q = hardeningparam(6)
              f = hardeningparam(7)
              ksub = hardeningparam(8)
          
          
              do is = 1,nslip
              
!                 Burgers vector
                  b = burgerv(is)
              
                  k2=k1*X*b/g*(1.-KBT/D/b**3*dlog(pdot/pdot0))
              
                  drhofor(is) = (k1*dsqrt(rhofor_t(is))-k2*rhofor_t(is))
     + *absgammadot(is)*dt    
              
                  drhosub = q*f*dsqrt(rhosub_t)* 
     + k2*rhofor_t(is)*absgammadot(is)*dt
              
              end do
                        
              
              
              
              
          end if
          
          
          
          


      end if

      return

      end subroutine hardeningrules
     
     

      

      
      
      end module hardening



