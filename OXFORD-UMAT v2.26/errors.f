!     Sept. 26th, 2022
!     Eralp Demir
!
!     This is written point the user to the sources of errors
!
      module errors
      implicit none
!
      contains
!
!     *******************************************************
!     *     This routine writes error statements in case    *
!     *     execution terminate or continue!                *
!     *******************************************************
      subroutine error(errorid)
      implicit none
!
!
!
      integer :: errorid
!
!
!     Message only when exiting the subroutine
!
      if(errorid.eq.1) then
          write(*,*) '*******************ERROR*******************'
          write(*,*) 'ERROR-001: Material-ID is out of range (1-9).
     +    Pls. check materials.f!'
          write(*,*) 'Exiting!'
          write(*,*) '*******************************************'
          call xit
      else if (errorid.eq.2) then
          write(*,*) '*******************ERROR*******************'
          write(*,*) 'ERROR-002: Element type is not defined! 
     +    Pls. check userinputs.f!'
          write(*,*) 'Exiting!'
          write(*,*) '*******************************************'
          call xit
      else if (errorid.eq.3) then
          write(*,*) '*******************ERROR*******************'
          write(*,*) 'ERROR-003: "cubicslip" integer parameter is not
     + properly defined. Pls. check userinputs.f!'
          write(*,*) 'Exiting!'
          write(*,*) '*******************************************'
          call xit
      else if (errorid.eq.4) then
          write(*,*) '*******************ERROR*******************'
          write(*,*) 'ERROR-004: phase id of the material is not
     + within the possible limits!. Pls. check PROPS variable!'
          write(*,*) 'Exiting!'
          write(*,*) '*******************************************'
          call xit
      else if (errorid.eq.5) then
          write(*,*) '*******************ERROR*******************'
          write(*,*) 'ERROR-005: "temperatureflag" is not
     + within the possible limits. Pls. check userinputs.f!'
          write(*,*) 'Exiting!'
          write(*,*) '*******************************************'
          call xit
      else if (errorid.eq.6) then
          write(*,*) '*******************ERROR*******************'
          write(*,*) 'ERROR-006: "NSTATV" is less than the
     + defined outputs.
     + Pls. check DEPVAR in ABAQUS!'
          write(*,*) 'Exiting!'
          write(*,*) '*******************************************'
      else if (errorid.eq.7) then
          write(*,*) '*******************ERROR*******************'
          write(*,*) 'ERROR-007: GND calculation is not possible
     + for the selected element type.
     + Pls. change the element type in ABAQUS!'
          write(*,*) 'Exiting!'
          write(*,*) '*******************************************'
          call xit
      else if (errorid.eq.8) then
          write(*,*) '*******************ERROR*******************'
          write(*,*) 'ERROR-008: Zero activation volume!
     + Pls. check the slip parameters and make sure that the
     + initial dislocation density has non-zero value
     + in usermaterials.f!'
          write(*,*) 'Exiting!'
          write(*,*) '*******************************************'
          call xit
      else if (errorid.eq.9) then
          write(*,*) '*******************ERROR*******************'
          write(*,*) 'ERROR-009: Could not read the element type
     + and the total number of elements from *.INP file.
     + Pls. check the file name in usermaterials.f!'
          write(*,*) 'Exiting!'
          write(*,*) '*******************************************'
      else if (errorid.eq.10) then
!          write(*,*) '*******************ERROR*******************'
!          write(*,*) 'ERROR-010: "Bmat" for L2 method
!     + is not invertible.
!     + GND model will give zero GND densities!.'
!          write(*,*) '*******************************************'
      else if (errorid.eq.11) then
!          write(*,*) '******************WARNING******************'
!          write(*,*) 'ERROR-006: WARNING DFGRD1  = NaN!'
!          write(*,*) 'Cut-back time by one-half!'
!          write(*,*) '*******************************************'
      else if (errorid.eq.12) then
!          write(*,*) '******************WARNING********************'
!          write(*,*) 'ERROR-012: Lp iteration did not converge'
!          write(*,*) 'Using forward-gradient Euler predictor'
!          write(*,*) '**********************************************'
      else if (errorid.eq.13) then
!          write(*,*) '******************WARNING********************'
!          write(*,*) 'ERROR-013: singular matrix in Lp iteration'
!          write(*,*) 'Will continue if alternate solution exists!'
!          write(*,*) '**********************************************'
      else if (errorid.eq.14) then
!          write(*,*) '******************WARNING********************'
!          write(*,*) 'ERROR-014: forward-gradient Euler predictor
!     + did not converge or it does not exist'
!          write(*,*) 'Cut-back time by one-half!'
!          write(*,*) '**********************************************'
      else if (errorid.eq.15) then
!          write(*,*) '******************WARNING********************'
!          write(*,*) 'ERROR-015: tmat array contains NaN elements'
!          write(*,*) 'Cut-back time by one-half!'
!          write(*,*) '**********************************************'
      else if (errorid.eq.16) then
!          write(*,*) '******************WARNING********************'
!          write(*,*) 'ERROR-016: NR scheme reached maximum number of
!     + iterations'
!          write(*,*) 'Cut-back time by one-half!'
!          write(*,*) '**********************************************'
      else if (errorid.eq.17) then
!          write(*,*) '******************WARNING********************'
!          write(*,*) 'ERROR-017: stress array contains NaN elements'
!          write(*,*) 'Cut-back time by one-half!'
!          write(*,*) '**********************************************'
      else if (errorid.eq.18) then
!          write(*,*) '******************WARNING********************'
!          write(*,*) 'ERROR-018: jacobian array contains NaN elements'
!          write(*,*) 'Cut-back time by one-half!'
!          write(*,*) '**********************************************'
      else if (errorid.eq.19) then
!          write(*,*) '******************WARNING********************'
!          write(*,*) 'ERROR-019: det(Fp) is close to zero'
!          write(*,*) 'Cut-back time by one-half!'
!          write(*,*) '**********************************************'
      else if (errorid.eq.20) then
!          write(*,*) '******************WARNING********************'
!          write(*,*) 'ERROR-020: det(Fe) is close to zero'
!          write(*,*) 'Cut-back time by one-half!'
!          write(*,*) '**********************************************'
      else if (errorid.eq.21) then
!          write(*,*) '******************WARNING********************'
!          write(*,*) 'ERROR-021: state variables become negative'
!          write(*,*) 'Cut-back time by one-half!'
!          write(*,*) '**********************************************'
      else if (errorid.eq.22) then
!          write(*,*) '******************WARNING********************'
!          write(*,*) 'ERROR-022: inversion in predictor did not work'
!          write(*,*) ''
!          write(*,*) '**********************************************'
      else
          write(*,*) '*******************ERROR*******************'
          write(*,*) 'ERROR-???: Unknown error!'
          write(*,*) 'Exiting!'
          write(*,*) '*******************************************'
          call xit
      end if
!
!
      return
      end subroutine error
!
!
!
      end module errors