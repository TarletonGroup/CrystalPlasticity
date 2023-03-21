!     Oct. 03rd, 2022
!     Eralp Demir
!
      module crss
      implicit none
      contains
!
!     ********************************************************
!     ** crss calculates the crss of slip systems           **
!     ********************************************************
      subroutine slipresistance(iphase,nslip,gf,G12,
     + burgerv,sintmat1,sintmat2,X,
     + tauc0,rhotot,rhofor,
     + rhosub,tausolute,loop,
     + irradiationmodel,irradiationparam,
     + temperature,
     + tauc)
      use userinputs, only: useaveragestatevars, 
     + maxnloop, maxnparam
      implicit none

!     crystal type
      integer,intent(in) :: iphase

!     number of slip systems
      integer,intent(in) :: nslip

!     geometric factor
      real(8),intent(in) :: gf

!     shear modulus for taylor dislocation law
      real(8),intent(in) :: G12
            
!     burgers vectors
      real(8),intent(in) :: burgerv(nslip)

!     Strength interaction matrices
!     Strength interaction between dislocations
      real(8), intent(in) :: sintmat1(nslip,nslip)
!     Strength interaction dislocation loops related with irradiation
      real(8), intent(in) :: sintmat2(nslip,nslip)      
      
!     fraction of hard obstacles
      real(8),intent(in) :: X      
      
!     critical resolved shear stress of slip systems
      real(8),intent(in) :: tauc0(nslip)

!     total scalar dislocation density (immobile)
      real(8),intent(in) :: rhotot

!     forest dislocation density
      real(8),intent(in) :: rhofor(nslip)

!     substructure dislocation density
      real(8),intent(in) :: rhosub

!     increase in tauc due to solute force
      real(8),intent(in) :: tausolute
      
!     increase in tauc due to loop defects
      real(8),intent(in) :: loop(maxnloop)
      
!     activate irradiation effect
      integer,intent(in) :: irradiationmodel

!     irradiation model parameters
      real(8),intent(in) :: irradiationparam(maxnloop) 
      
!     current temperature
      real(8),intent(in) :: temperature

!     critical resolved shear stress of slip systems
      real(8),intent(out) :: tauc(nslip)

!     variables used in this subroutine
      integer :: is, il
      real(8) :: Ploop(nslip)
      
!     Reset the density of loops
      Ploop = 0.

!     If irradiation model-2 is set ON
!     Compute the strength contribution
      if (irradiationmodel==2) then
          
          do is = 1, nslip
              
              
              do il = 1, 3
          
                  
                  Ploop(is) = Ploop(is) + 
     + sintmat2(is,il)**2. * loop(il)
                  
                  
              end do
              
          end do
          
          
          
      end if
      
      
      
      
      
      

!     Check crystal type
!     Only for alpha-uranium there is an exception

          
         
!     Defined by data fitting
!     as a function of rhofor, rhosub, and temperature
      if (iphase == 4) then

!         alpha-uranium model with forest and substructure dislocations
!         R.J. McCabe, L. Capolungo, P.E. Marshall, C.M. Cady, C.N. Tome
!         Deformation of wrought uranium: experiments and modeling
!         Acta Materialia 58 (2010) 5447–5459
          tauc(1) = tauc0(1) + 19.066 * dsqrt(rhofor(1)) +
     +    1.8218 * dsqrt(rhosub) * 
     +    log(1. / burgerv(1) / dsqrt(rhosub))
          tauc(2) = tauc0(2) + 18.832 * dsqrt(rhofor(2)) +
     +    1.7995 * dsqrt(rhosub) * 
     +    log(1. / burgerv(2) / dsqrt(rhosub))
          tauc(3) = tauc0(3) + 54.052 * dsqrt(rhofor(3)) + 
     +    5.1650 * dsqrt(rhosub) * 
     +    log(1. / burgerv(3) / dsqrt(rhosub))
          tauc(4) = tauc0(4) + 54.052 * dsqrt(rhofor(4)) + 
     +    5.1650 * dsqrt(rhosub) * 
     +    log(1. / burgerv(4) / dsqrt(rhosub))
          tauc(5) = tauc0(5) + 123.357 * dsqrt(rhofor(5)) + 
     +    11.7875 * dsqrt(rhosub) * 
     +    log(1. / burgerv(5) / dsqrt(rhosub))
          tauc(6) = tauc0(6) + 123.357 * dsqrt(rhofor(6)) + 
     +    11.7875 * dsqrt(rhosub) * 
     +    log(1. / burgerv(6) / dsqrt(rhosub))
          tauc(7) = tauc0(7) + 123.357 * dsqrt(rhofor(7)) + 
     +    11.7875 * dsqrt(rhosub) * 
     +    log(1. / burgerv(7) / dsqrt(rhosub))
          tauc(8) = tauc0(8) + 123.357 * dsqrt(rhofor(8)) + 
     +    11.7875 * dsqrt(rhosub) * 
     +    log(1. / burgerv(8) / dsqrt(rhosub))

          ! zecevic 2016 temperature dependence
          tauc(1) = tauc(1) * dexp(-(temperature-293.0)/140.0)
          tauc(2) = tauc(2) * dexp(-(temperature-293.0)/140.0)
          tauc(3) = tauc(3) * dexp(-(temperature-293.0)/140.0)
          tauc(4) = tauc(4) * dexp(-(temperature-293.0)/140.0)
          tauc(5) = tauc(5) * dexp(-(temperature-293.0)/140.0)
          tauc(6) = tauc(6) * dexp(-(temperature-293.0)/140.0)
          tauc(7) = tauc(7) * dexp(-(temperature-293.0)/140.0)
          tauc(8) = tauc(8) * dexp(-(temperature-293.0)/140.0)

          ! daniel, lesage, 1971 minimum value as minima
          tauc(1) = max(tauc(1),4.0)
          tauc(2) = max(tauc(2),4.0)
          tauc(3) = max(tauc(3),4.0)
          tauc(4) = max(tauc(4),4.0)
          tauc(5) = max(tauc(5),4.0)
          tauc(6) = max(tauc(6),4.0)
          tauc(7) = max(tauc(7),4.0)
          tauc(8) = max(tauc(8),4.0)


          
!     For any other phase
!     User forest spacing to determine the strength
      else    
          
          
           if (useaveragestatevars == 0) then
              
              do is = 1, nslip


!                 Taylor strength - forest spacing
                  tauc(is) = tauc0(is) * (1.-X) + 
     + G12*burgerv(is)*dsqrt(X*tauc0(is)/G12/burgerv(is) +
     + gf**2.*rhofor(is) + Ploop(is))
                  
              end do


              
          elseif (useaveragestatevars == 1) then
              

              do is = 1, nslip
                  
!                 Taylor strength - total density
                  tauc(is) = tauc0(is)*(1.-X) +
     + G12*burgerv(is)*dsqrt(X*tauc0(is)/G12/burgerv(is) + 
     + gf**2.*rhotot + Ploop(is))
              
              end do
              
          endif    

          
          
          
          
          
          
      end if ! check crystal type

      
      
!     Effect of irradiation
!     True for all crystal types
           
      if (irradiationmodel == 1) then
          tauc = tauc + tausolute
      end if
      
      
      
      return

      end subroutine slipresistance

     
      end module crss