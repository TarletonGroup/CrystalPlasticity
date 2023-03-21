
!     Global variables used for irradiation model-2
      module irradiation
      implicit none      
      
      
!     This is used in irradiationmodel=2 only
!     Types of defects : 3
!     Strength interaction coefficients for dislocation-defect interactions in HCP
      real(8), public ::  sintmat_hcp_irm2(30,3)
      data sintmat_hcp_irm2(1,:)  /1.25,1.25,1.25/
      data sintmat_hcp_irm2(2,:)  /1.25,1.25,1.25/
      data sintmat_hcp_irm2(3,:)  /1.25,1.25,1.25/
           
      data sintmat_hcp_irm2(4,:)  /1.25,1.875,1.875/
      data sintmat_hcp_irm2(5,:)  /1.875,1.25,1.875/
      data sintmat_hcp_irm2(6,:)  /1.875,1.875,1.25/
           
      data sintmat_hcp_irm2(7,:)  /1.25,1.875,1.875/
      data sintmat_hcp_irm2(8,:)  /1.875,1.25,1.875/
      data sintmat_hcp_irm2(9,:)  /1.875,1.875,1.25/
      data sintmat_hcp_irm2(10,:) /1.25,1.875,1.875/
      data sintmat_hcp_irm2(11,:) /1.875,1.25,1.875/
      data sintmat_hcp_irm2(12,:) /1.875,1.875,1.25/
           
      data sintmat_hcp_irm2(13,:) /1.875,1.25,1.875/
      data sintmat_hcp_irm2(14,:) /1.875,1.875,1.25/
      data sintmat_hcp_irm2(15,:) /1.875,1.25,1.875/
      data sintmat_hcp_irm2(16,:) /1.875,1.875,1.25/
      data sintmat_hcp_irm2(17,:) /1.25,1.875,1.875/
      data sintmat_hcp_irm2(18,:) /1.25,1.875,1.875/

      data sintmat_hcp_irm2(19,:) /1.875,1.25,1.875/
      data sintmat_hcp_irm2(20,:) /1.875,1.875,1.25/
      data sintmat_hcp_irm2(21,:) /1.875,1.25,1.875/
      data sintmat_hcp_irm2(22,:) /1.875,1.875,1.25/
      data sintmat_hcp_irm2(23,:) /1.25,1.875,1.875/
      data sintmat_hcp_irm2(24,:) /1.25,1.875,1.875/
           
      data sintmat_hcp_irm2(25,:) /1.875,1.875,1.875/
      data sintmat_hcp_irm2(26,:) /1.875,1.875,1.875/
      data sintmat_hcp_irm2(27,:) /1.875,1.875,1.875/
      data sintmat_hcp_irm2(28,:) /1.875,1.875,1.875/
      data sintmat_hcp_irm2(29,:) /1.875,1.875,1.875/
      data sintmat_hcp_irm2(30,:) /1.875,1.875,1.875/

           
           
           
!     This is used in irradiationmodel=2 only           
!     Annihilation Constants for defect-dislocation interactions in HCP
!     Softening (negative hardening) interaction coefficients
      real(8), public ::  hintmat_hcp_irm2(3,30)
      data hintmat_hcp_irm2(:,1)  /0.5,0.5,0.5/
      data hintmat_hcp_irm2(:,2)  /0.5,0.5,0.5/
      data hintmat_hcp_irm2(:,3)  /0.5,0.5,0.5/

      data hintmat_hcp_irm2(:,4)  /0.5,0.0,0.0/
      data hintmat_hcp_irm2(:,5)  /0.0,0.5,0.0/
      data hintmat_hcp_irm2(:,6)  /0.0,0.0,0.5/
          
      data hintmat_hcp_irm2(:,7)  /0.5,0.0,0.0/
      data hintmat_hcp_irm2(:,8)  /0.0,0.5,0.0/
      data hintmat_hcp_irm2(:,9)  /0.0,0.0,0.5/
      data hintmat_hcp_irm2(:,10) /0.5,0.0,0.0/
      data hintmat_hcp_irm2(:,11) /0.0,0.5,0.0/
      data hintmat_hcp_irm2(:,12) /0.0,0.0,0.5/

      data hintmat_hcp_irm2(:,13) /0.0,0.5,0.0/
      data hintmat_hcp_irm2(:,14) /0.0,0.0,0.5/
      data hintmat_hcp_irm2(:,15) /0.0,0.5,0.0/
      data hintmat_hcp_irm2(:,16) /0.0,0.0,0.5/
      data hintmat_hcp_irm2(:,17) /0.5,0.0,0.0/
      data hintmat_hcp_irm2(:,18) /0.5,0.0,0.0/

      data hintmat_hcp_irm2(:,19) /0.0,0.5,0.0/
      data hintmat_hcp_irm2(:,20) /0.0,0.0,0.5/
      data hintmat_hcp_irm2(:,21) /0.0,0.5,0.0/
      data hintmat_hcp_irm2(:,22) /0.0,0.0,0.5/
      data hintmat_hcp_irm2(:,23) /0.5,0.0,0.0/
      data hintmat_hcp_irm2(:,24) /0.5,0.0,0.0/
      
      data hintmat_hcp_irm2(:,25) /0.0,0.0,0.0/
      data hintmat_hcp_irm2(:,26) /0.0,0.0,0.0/
      data hintmat_hcp_irm2(:,27) /0.0,0.0,0.0/
      data hintmat_hcp_irm2(:,28) /0.0,0.0,0.0/
      data hintmat_hcp_irm2(:,29) /0.0,0.0,0.0/
      data hintmat_hcp_irm2(:,30) /0.0,0.0,0.0/
           
      
      
      
      end module irradiation