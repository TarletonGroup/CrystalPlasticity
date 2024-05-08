!     Nov. 11th, 2022
!     Eralp Demir
!
!     This is written to define the outputs for post-processing
!
      module useroutputs
      implicit none
!
!
!     GLOBAL VARIABLES USED IN POST-PROCESSING
!     -------------------------------------------------------------------------------
!
!     Flags for different outputs (22-30 custom outputs)
      integer, public :: statev_outputs(30)
!
!     Set the desired outputs by setting  0 / 1 
!     in the "defineoutputs" subroutine in the below!
!
!
!
!     Number of outputs - will be computed
      integer, public :: nstatv_outputs
!
!
!     -------------------------------------------------------------------------------
!
      contains
!
!     *******************************************************
!     *     This routine writes error statements in case    *
!     *     execution terminate or continue!                *
!     *******************************************************
      subroutine defineoutputs
!
      use userinputs, only : maxnslip, maxnloop
!
      implicit none
!
!     Variables used in the subroutine
      integer :: i, ind
!
!
!     List of variables
!     Set the flag to "1" if requested
!
!     1st State-variable output / number of outputs: 9
!     Crystal to Sample tranformation: statev_gmatinv
      statev_outputs(1) = 0
!
!     2nd State-variable output / number of outputs: 1
!     Equivalent Von-Mises plastic total strain: statev_evmp
      statev_outputs(2) = 0
!
!     3rd State-variable output / number of outputs: 1
!     Maximum ratio of rss to crss: statev_maxx
      statev_outputs(3) = 0
!
!     4th State-variable output / number of outputs: 6
!     Elastic strains in the crystal frame: statev_Eec
      statev_outputs(4) = 0
!
!     5th State-variable output / number of outputs: 9
!     Lattice curvature: statev_curvature
      statev_outputs(5) = 0
!
!     6th State-variable output / number of outputs: 1
!     Total statistically-stored dislocation density: statev_ssdtot
      statev_outputs(6) = 0
!
!     7th State-variable output / number of outputs: 1
!     Substructure dislocation density: statev_substructure
      statev_outputs(7) = 0
!
!     8th State-variable output / number of outputs: 1
!     Solute strength: statev_tausolute
      statev_outputs(8) = 0
!
!     9th State-variable output / number of outputs: 1
!     Cumulative slip: statev_totgammasum
      statev_outputs(9) = 0
!
!     10th State-variable output / number of outputs: maxnslip
!     Total slip per slip system: statev_gammasum
      statev_outputs(10) = 1
!
!     11th State-variable output / number of outputs: maxnslip
!     Slip rate: statev_gammadot
      statev_outputs(11) = 0
!
!     12nd State-variable output / number of outputs: maxnslip
!     Critical Resolved Shear Stress: statev_tauc
      statev_outputs(12) = 0
!
!     13rd State-variable output / number of outputs: maxnslip
!     Statistically-Stored Dislocation Density: statev_ssd
      statev_outputs(13) = 0
!
!     14th State-variable output / number of outputs: maxnslip*2
!     Geometrically Necessary Dislocation Density: statev_gnd
      statev_outputs(14) = 0
!
!     15th State-variable output / number of outputs: maxnslip
!     Foresty Dislocation Density: statev_forest
      statev_outputs(15) = 0
!
!     16th State-variable output / number of outputs: maxnloop
!     Defect Loop Density: statev_loop
      statev_outputs(16) = 0
!
!     17th State-variable output / number of outputs: maxnslip
!     Backstress: statev_backstress
      statev_outputs(17) = 0
!
!     18th State-variable output / number of outputs: 1
!     Total GND density
      statev_outputs(18) = 0
!
!     19th State-variable output / number of outputs: 1
!     Plastic dissipation power density
      statev_outputs(19) = 0
!     
!     20st State-variable output / number of outputs: 1
!     Fatemi Socie parameter
      statev_outputs(20) = 0
!
!     21-30 custom outputs
!     Need to be defined here!
      statev_outputs(21) = 0
      statev_outputs(22) = 0
      statev_outputs(23) = 0
      statev_outputs(24) = 0
      statev_outputs(25) = 0
      statev_outputs(26) = 0
      statev_outputs(27) = 0
      statev_outputs(28) = 0
      statev_outputs(29) = 0
      statev_outputs(30) = 0
!
!
!
!
!
!
!     
!     Count the user-defined outputs
      nstatv_outputs=0
!     Find the total number of state variable outputs
      if (statev_outputs(1)==1) then
          nstatv_outputs=nstatv_outputs+9
      endif
      if (statev_outputs(2)==1) then
          nstatv_outputs=nstatv_outputs+1
      endif        
      if (statev_outputs(3)==1) then
          nstatv_outputs=nstatv_outputs+1
      endif      
      if (statev_outputs(4)==1) then
          nstatv_outputs=nstatv_outputs+6
      endif
      if (statev_outputs(5)==1) then
          nstatv_outputs=nstatv_outputs+9
      endif
      if (statev_outputs(6)==1) then
          nstatv_outputs=nstatv_outputs+1
      endif
      if (statev_outputs(7)==1) then
          nstatv_outputs=nstatv_outputs+1
      endif            
      if (statev_outputs(8)==1) then
          nstatv_outputs=nstatv_outputs+1
      endif            
      if (statev_outputs(9)==1) then
          nstatv_outputs=nstatv_outputs+1
      endif    
      if (statev_outputs(10)==1) then
          nstatv_outputs=nstatv_outputs+maxnslip
      endif
      if (statev_outputs(11)==1) then
          nstatv_outputs=nstatv_outputs+maxnslip
      endif         
      if (statev_outputs(12)==1) then
          nstatv_outputs=nstatv_outputs+maxnslip
      endif         
      if (statev_outputs(13)==1) then
          nstatv_outputs=nstatv_outputs+maxnslip
      endif                
      if (statev_outputs(14)==1) then
          nstatv_outputs=nstatv_outputs+maxnslip*2
      endif       
      if (statev_outputs(15)==1) then
          nstatv_outputs=nstatv_outputs+maxnslip
      endif     
      if (statev_outputs(16)==1) then
          nstatv_outputs=nstatv_outputs+maxnloop
      endif     
      if (statev_outputs(17)==1) then
          nstatv_outputs=nstatv_outputs+maxnslip
      endif
      if (statev_outputs(18)==1) then
          nstatv_outputs=nstatv_outputs+1
      endif
      if (statev_outputs(19)==1) then
          nstatv_outputs=nstatv_outputs+1
      endif
      if (statev_outputs(20)==1) then
          nstatv_outputs=nstatv_outputs+1
      endif
!
!
!     Custom outputs need to be filled here!      
!
!
!
!
!
      end subroutine defineoutputs
!
!
!
!
!     *******************************************************
!     *     This routine writes error statements in case    *
!     *     execution terminate or continue!                *
!     *******************************************************
      subroutine checkoutputs(nstatv)
      use errors, only: error
      implicit none
!     Number of state variables
      integer, intent(in) :: nstatv
!
!     Check if defined correctly
      if (nstatv_outputs > nstatv) then
          write(*,*) 'There are ', 
     + nstatv_outputs, ' number of outputs!'
          write(*,*) 'But, ', 
     + nstatv, ' number of outputs were defined in DEPVAR!'
!         error message in .dat file
          call error(6)
      end if
!
      end subroutine checkoutputs
!
!
!
!
!
!
!
      subroutine write_statev_legend
      use userinputs, only : maxnslip, maxnloop
!
	implicit none
!
      integer count, i
      character*2 ij
!
!     Write the legend of the output variables (nstatv)
!     Outputs to extract: 
!     1: Transformation matrix from crystal to sample (x9)
!     2: Equivalent Von-Mises plastic strain (x1)
!     3: Maximum ratio of rss to crss (x1)
!     4: Elastic strain in crystal frame (x6)
!     5: Lattice curvature (x9)
!     6: SSD total (x1)
!     7: Substructure density (x1)
!     8: Solute strength (x1)
!     9: cumulative slip (x1)
!     10: total slip per slip system (x maxnslip)
!     11: slip rates per slip system (x maxnslip)
!     12: CRSS (x maxnslip)
!     13: SSD (x maxnslip)
!     14: GND (x maxnslip)
!     15: Forest (x maxnslip)
!     16: Loop (x maxnloop)
!     17: Backstress (x maxnslip)
!     18: GND density (x1)
!     19: Plastic dissipation power density (x1)
!     20: Fatemi Socie parameter (x1)
!     21-30: Custom outputs
!
!
!
!
!
!
      count = 0
      open(100,file='../STATEV_legend.txt',action='write',
     + status='replace')
!
!
!
!
!
!
!
!     State variable-1
      if (statev_outputs(1) == 1) then
!
          do i = 1, 9
!
              count = count + 1
!
              select case(i)
!
              case(1)
                  ij = 'xx'
              case(2)
                  ij = 'xy'
              case(3)
                  ij = 'xz'
              case(4)
                  ij = 'yx'
              case(5)
                  ij = 'yy'
              case(6)
                  ij = 'yz'
              case(7)
                  ij = 'zx'
              case(8)
                  ij = 'zy'
              case(9)
                  ij = 'zz'  
!
              end select
!
!
!
              write(100,'(A7,I3,A37,A2,A4)')
     + 'STATEV-', count, 
     + ':  Crystal to Sample transformation -',  ij,  ' [-]'
!
          end do
!
!
      endif
!
!
!
!
!     State variable-2
      if (statev_outputs(2) == 1) then
!
!
!
          count = count + 1
!
!
!
          write(100,'(A7,I3,A39,A4)')
     + 'STATEV-', count,
     + ':   Equivalent Von-Mises plastic strain', ' [-]'
!
!
!
!
      endif
! 
!
!     State variable-3
      if (statev_outputs(3) == 1) then
!
!
!
          count = count + 1
!
!
!
!
          write(100,'(A7,I3,A32,A4)')
     + 'STATEV-', count, 
     + ':   Maximum ratio of rss to crss', ' [-]'
!              
!
!
!
      endif
!
!
!
!     State variable-4
      if (statev_outputs(4) == 1) then
!
          do i = 1, 6
!
              count = count + 1
!
              select case(i)
!
              case(1)
                  ij = 'xx'
              case(2)
                  ij = 'yy'
              case(3)
                  ij = 'zz'
              case(4)
                  ij = 'xy'
              case(5)
                  ij = 'xz'
              case(6)
                  ij = 'yz'                  
!
              end select
!
!             
!
              write(100,'(A7,I3,A36,A2,A6)')
     + 'STATEV-', count, 
     + ':   Elastic strain in crystal frame-',  ij,  ' [MPa]'
!
          end do
!
!
      endif
!
!
!
!
!
!     State variable-5
      if (statev_outputs(5) == 1) then
!
          do i = 1, 9
!
              count = count + 1
!
              select case(i)
!
              case(1)
                  ij = 'xx'
              case(2)
                  ij = 'xy'
              case(3)
                  ij = 'xz'
              case(4)
                  ij = 'yx'
              case(5)
                  ij = 'yy'
              case(6)
                  ij = 'yz'
              case(7)
                  ij = 'zx'
              case(8)
                  ij = 'zy'
              case(9)
                  ij = 'zz'  
!
              end select
!
!
!
              write(100,'(A7,I3,A21,A2,A15)')
     + 'STATEV-', count, 
     + ':   Lattice curvature',  ij,  ' [1/micrometer]'
!
          end do
!
!
      endif
!
!
!
!
!
!
!
!
!
!     State variable-6   
      if (statev_outputs(6) == 1) then
!
!
!
          count = count + 1
!
!
!
!
          write(100,'(A7,I3,A21,A17)')
     + 'STATEV-', count, 
     + ':   total SSD density', ' [1/micrometer^2]'
!
!
!
!
      endif
!
!
!
!     State variable-7
      if (statev_outputs(7) == 1) then
!
!
!
          count = count + 1
!
!
!
!
          write(100,'(A7,I3,A24,A17)')
     + 'STATEV-', count, 
     + ':   substructure density', ' [1/micrometer^2]'
!
!
!
!
      endif 
!
!
!     State variable-8  
      if (statev_outputs(8) == 1) then
!
!
!
          count = count + 1
!
!
!
          write(100,'(A7,I3,A19,A6)')
     + 'STATEV-', count, 
     + ':   solute strength', ' [MPa]'
!
!
!
!
      endif
! 
!
!     State variable-9 
      if (statev_outputs(9) == 1) then
!
!
!
          count = count + 1
!
!
!
!
          write(100,'(A7,I3,A19,A4)')
     + 'STATEV-', count, 
     + ':   cumulative slip', ' [-]'
!
!
!
!
      endif
!
!
!
!
!     State variable-10
      if (statev_outputs(10) == 1) then
!
          do i = 1, maxnslip
!
              count = count + 1
!
!
!
              write(100,'(A7,I3,A30,I2,A4)')
     + 'STATEV-', count, 
     + ':   Total slip on slip system-',  i,  ' [-]'
!
          end do
!
!
      endif
!
!
!
!
!     State variable-11
      if (statev_outputs(11) == 1) then
!
          do i = 1, maxnslip
!
              count = count + 1
!
!
!
              write(100,'(A7,I3,A29,I2,A6)')
     + 'STATEV-', count, 
     + ':   Slip rate of slip system-',  i,  ' [1/s]'
!
          end do
!
!
      endif
!
!
!     State variable-12
      if (statev_outputs(12) == 1) then
!
          do i = 1, maxnslip
!
              count = count + 1
!
!
!
              write(100,'(A7,I3,A24,I2,A6)')
     + 'STATEV-', count, 
     + ':   CRSS on slip system-',  i,  ' [MPa]'
!
          end do
!
!
      endif
!
!
!
!
!     State variable-13
      if (statev_outputs(13) == 1) then
!
          do i = 1, maxnslip
!
              count = count + 1
!
!
!
              write(100,'(A7,I3,A31,I2,A17)')
     + 'STATEV-', count, 
     + ':   SSD density of slip system-',  i,  ' [1/micrometer^2]'
!
          end do
!
!
      endif
!
!
!
!
!     State variable-14
      if (statev_outputs(14) == 1) then
!
          do i = 1, maxnslip*2
!
              count = count + 1
!
!             Edge dislocation
              if (i <= maxnslip) then
!
                  write(100,'(A7,I3,A31,I2,A17)')
     + 'STATEV-', count, 
     + ':   Edge dislocation density-',  i,  ' [1/micrometer^2]'
!
              else
!
                  write(100,'(A7,I3,A31,I2,A17)')
     + 'STATEV-', count, 
     + ':   Screw dislocation density-', i-maxnslip, ' [1/micrometer^2]'
!
              end if
!
!
          end do
!
!
      endif
!
!
!     State variable-15
      if (statev_outputs(15) == 1) then
!
          do i = 1, maxnslip
!
              count = count + 1
!
!
!
              write(100,'(A7,I3,A46,I2,A17)')
     + 'STATEV-', count, 
     + ':   Forest dislocation density of slip system-',
     + i,  ' [1/micrometer^2]'
!
          end do
!
!
      endif
!
!
!     State variable-16
      if (statev_outputs(16) == 1) then
!
          do i = 1, maxnloop
!
              count = count + 1
!
!
!
              write(100,'(A7,I3,A46,I2,A17)')
     + 'STATEV-', count, 
     + ':   Loop dislocation density of defect system-',
     + i,  ' [1/micrometer^2]'
!
          end do
!
!
      endif
!
!
!     State variable-17
      if (statev_outputs(17) == 1) then
!
          do i = 1, maxnslip
!
              count = count + 1
!
!
!
              write(100,'(A7,I3,A30,I2,A6)')
     + 'STATEV-', count, 
     + ':   Backstress on slip system-',
     + i,  ' [MPa]'
!
          end do
!
!
      endif
!
!     State variable-18
      if (statev_outputs(18) == 1) then
!
!
!
          count = count + 1
!
!
!
!
          write(100,'(A7,I3,A21,A17)')
     + 'STATEV-', count, 
     + ':   Total GND density', ' [1/micrometer^2]'
!
!
!
!
      endif
!
!
!     State variable-19
      if (statev_outputs(19) == 1) then
!
!
!
          count = count + 1
!
!
!
!
          write(100,'(A7,I3,A37,A5)')
     + 'STATEV-', count, 
     + ':   Plastic dissipation power density', ' [nW]'
!
!
!
!
      endif
!
!     State variable-20
      if (statev_outputs(20) == 1) then
!
!
!
          count = count + 1
!
!
!
!
          write(100,'(A7,I3,A26,A4)')
     + 'STATEV-', count, 
     + ':   Fatemi Socie parameter', ' [-]'
!
!
!
!
      endif
!
      close(100)
!
!
!
!
!
!
	return
      end subroutine write_statev_legend
!
!
!
!
!
!
!
      subroutine assignoutputs(noel,npt,nstatv,statev)
      use globalvariables, only: statev_gmatinv,
     + statev_evmp, statev_maxx, statev_Eec,
     + statev_curvature, statev_backstress_t,
     + statev_ssdtot, statev_substructure,
     + statev_tausolute, statev_totgammasum,
     + statev_gammasum, statev_gammadot,
     + statev_tauceff, statev_ssd, statev_loop, statev_gnd_t,
     + statev_forest, statev_plasdiss
      use userinputs, only: maxnslip, maxnloop
      use utilities, only: matvec9
      implicit none
!     Element number
      integer, intent(in) :: noel
!     Integration point
      integer, intent(in) :: npt
!     Number of state variables
      integer, intent(in) :: nstatv
!     Values of state variables
      real(8), intent(inout) :: statev(nstatv)
!     Other variables
      integer :: i, j
      real(8) :: d6(6), d9(9), d1
      real(8) :: gnd(maxnslip*2)
!     Outputs to extract: 
!     1: Transformation matrix from crystal to sample (x9)
!     2: Equivalent Von-Mises plastic strain (x1)
!     3: Maxium ratio of rss to crss (x1)
!     4: Elastic strain in crystal frame (x6)
!     5: Lattice curvature (x9)
!     6: SSD total (x1)
!     7: Substructure density (x1)
!     8: Solute strength (x1)
!     9: Cumulative slip (x1)
!     10: Total slip per slip system (x maxnslip)
!     11: Slip rates per slip system (x maxnslip)
!     12: Effective CRSS (x maxnslip)
!     13: SSD (x maxnslip) 
!     14: GND (x maxnslip) 
!     15: Forest (x maxnslip)
!     16: Loop density (x maxnloop)
!     17: Backstress (x maxnslip)
!     18: Total GND density (x1)
!     19: Plastic dissipation (x1)
!     20: Fatemi Socie parameter (x1)      
!     21-30: Custom outputs
!
!
!     Reset the counter
      i=0
!
!     State variable-1
!     Transformation matrix from crystal to sample (x9)
      if (statev_outputs(1)==1) then
!
!
!
          call matvec9(statev_gmatinv(noel,npt,:,:),d9)
!
          do j = 1, 9
!
              i = i + 1
              statev(i) = d9(j)
!
          end do
!
      end if
!
!
!
!     State variable-2
!     Equivalent Von-Mises plastic strain (x1)
      if (statev_outputs(2)==1) then
!
          i = i + 1
          statev(i) = statev_evmp(noel,npt)
!
      end if
!
!
!     State variable-3
!     Maxium ratio of rss to crss (x1)
      if (statev_outputs(3)==1) then
!
          i = i + 1
          statev(i) = statev_maxx(noel,npt)
!
      end if
!
!
!     State variable-4
!     Elastic strain in crystal frame (x6)
      if (statev_outputs(4)==1) then
!
!
          do j = 1, 6
              i = i + 1
              statev(i) = statev_Eec(noel,npt,j)
          end do
!
!
      end if
!
!
!     State variable-5
!     Lattice curvature (x9)
      if (statev_outputs(5)==1) then
!
          d9 = statev_curvature(noel,npt,:)
!
          do j = 1, 9
!
              i = i + 1
              statev(i) = d9(j)
!
          end do
!
!
      end if
!
!
!     State variable-6
!     SSD total (x1)
      if (statev_outputs(6)==1) then
!
          i = i + 1
          statev(i) = statev_ssdtot(noel,npt)     
!
      end if
!
!
!
!     State variable-7
!     Substructure density (x1)
      if (statev_outputs(7)==1) then
!
          i = i + 1
          statev(i) = statev_substructure(noel,npt)     
!
      end if
!
!
!     State variable-8
!     Solute strength (x1)
      if (statev_outputs(8)==1) then
!
          i = i + 1
          statev(i) = statev_tausolute(noel,npt)     
!
      end if
!
!
!     State variable-9
!     Cumulative slip (x1)
      if (statev_outputs(9)==1) then
!
          i = i + 1
          statev(i) = statev_totgammasum(noel,npt)     
!
      end if
!
!
!     State variable-10
!     Total slip per slip system (x maxnslip)
      if (statev_outputs(10)==1) then
!
          do j = 1, maxnslip
              i = i + 1
              statev(i) = statev_gammasum(noel,npt,j)
          end do
!
      end if
!
!
!     State variable-11
!     Slip rates per slip system (x maxnslip)
      if (statev_outputs(11)==1) then
!
          do j = 1, maxnslip
              i = i + 1
              statev(i) = statev_gammadot(noel,npt,j)
          end do
!
      end if
!
!     State variable-12
!     CRSS (x maxnslip)
      if (statev_outputs(12)==1) then
!
          do j = 1, maxnslip
              i = i + 1
              statev(i) = statev_tauceff(noel,npt,j)
          end do
!      
      end if    
!
!
!     State variable-13
!     SSD (x maxnslip) 
      if (statev_outputs(13)==1) then
!
          do j = 1, maxnslip
              i = i + 1
              statev(i) = statev_ssd(noel,npt,j)
          end do
!  
      end if
!
!
!     State variable-14
!     GND (x maxnslip)       
      if (statev_outputs(14)==1) then
!
          do j = 1, maxnslip*2
              i = i + 1
              statev(i) = statev_gnd_t(noel,npt,j)
          end do
!
      end if
!
!
!     State variable-15
!     Forest density (x maxnslip)       
      if (statev_outputs(15)==1) then
!
          do j = 1, maxnslip
              i = i + 1
              statev(i) = statev_forest(noel,npt,j)
          end do
!
      end if
!
!     State variable-16
!     Loop density (x maxnloop)       
      if (statev_outputs(16)==1) then
!
          do j = 1, maxnloop
              i = i + 1
              statev(i) = statev_loop(noel,npt,j)
          end do
!
      end if
!
!
!     State variable-17
!     Backstress (x maxnslip)       
      if (statev_outputs(17)==1) then
!
          do j = 1, maxnslip
              i = i + 1
              statev(i) = statev_backstress_t(noel,npt,j)
          end do
!
      end if
!
!     State variable-18
!     GND total (x1)
      if (statev_outputs(18)==1) then
!
          i = i + 1
          gnd = statev_gnd_t(noel,npt,1:2*maxnslip)
          statev(i) = sqrt(sum(gnd*gnd))
!
      end if
!
!     State variable-19
!     Plastic dissipation power density (x1)
      if (statev_outputs(19)==1) then
!
          i = i + 1
          statev(i) = statev_plasdiss(noel,npt)
!
      end if
!
!     State variable-20
!     Fatemi Socie parameter (x1)
      if (statev_outputs(20)==1) then
!
          i = i + 1
          call  FatemiSocie_parameter(noel,npt,d1)
          statev(i) = d1
!
      end if
!
!
!     Custom state variables 21-30
!     Please add custom outputs here!
!
!
!
!
      return
      end subroutine assignoutputs
!
!
!!
!
!
!
!     Subroutine to calculate Fatemi Socie parameter
      subroutine FatemiSocie_parameter(noel,npt,FSp)
      use userinputs, only: maxnslip
      use globalvariables, only: statev_sigma, statev_gammasum,
     + statev_gmatinv, statev_tauceff, materialid, norc_0_all
      use utilities, only: vecmat6
      implicit none
      integer, intent(in) :: noel
      integer, intent(in) :: npt
      real(8), intent(out) :: FSp
!     Variables used in the subroutine
      real(8) :: sig6(6), sig3x3(3,3)
      real(8) :: gammasum(maxnslip)
      real(8) :: gmatinv(3,3)
      real(8) :: tauceff(maxnslip), tauceffmax
      integer :: matid
      real(8) :: norc(3), nors(3)
      real(8) :: shmax, sigmax, k
      integer :: ismax, is, i, j
!
!     Input parameter "k"
!     Reference: https://doi.org/10.1016/j.ijfatigue.2011.01.003
      k = 0.5
!
!
      FSp=0.
      sig6 = statev_sigma(noel,npt,:)
      gammasum = statev_gammasum(noel,npt,:)
      gmatinv = statev_gmatinv(noel,npt,:,:)
      tauceff = statev_tauceff(noel,npt,:)
      matid = materialid(noel,npt)
!
!
!     Find maximum slip
      shmax=0.
      do is=1,maxnslip
          if (abs(gammasum(is)).gt.shmax) then
              shmax = abs(gammasum(is))
              ismax = is
          end if
      end do
!
!     Maximum CRSS
      tauceffmax = tauceff(ismax)
!
!     Slip plane normal of maximum slip
      norc = norc_0_all(matid,ismax,:)
!
!     Transform to sample reference
      nors = matmul(gmatinv,norc)
!
!     Convert stress to 3x3 matrix
      call vecmat6(sig6,sig3x3)
!
!     Project stress to the slip plane of maximum slip
      sigmax = 0.
      do i = 1, 3
          do j = 1, 3
              sigmax = sigmax +
     + nors(i)*sig3x3(i,j)*nors(j)
          end do
      end do
!
!     Parameter calculation
      FSp = shmax/2.*(1.+k*sigmax/tauceffmax) 
!             
!
      return
      end subroutine FatemiSocie_parameter
!
!
!
      end module useroutputs