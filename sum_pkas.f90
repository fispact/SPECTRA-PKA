!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk,                                                     *
!******************************************************************************
!end rubric

  SUBROUTINE sum_pkas
  use globals
   IMPLICIT NONE
   
   REAL (KIND=DBL), ALLOCATABLE :: tdam_pka_temp(:),pka_temp(:)  


    !24/4/2018 changes to allow for user output grid for sums and globals
    ! also removes error check for mis-matched grids - now we can just proceed and flag
    IF(num_pka_recoil_points.NE.num_pka_recoil_points_master) THEN
     IF (.NOT.do_user_output_energy_grid) THEN
      PRINT *,'grid of recoil energies do not match while summing over mt numbers'
      PRINT *,'interpolation will be used to collapse the current grid onto the mt-sum one'
     END IF
    END IF
    ALLOCATE(pka_temp(MAX(num_pka_recoil_points_master,num_pka_recoil_points)))


    pka_temp(1:num_pka_recoil_points)=pka(1,:)
    CALL collapse_fluxes(pka_temp, &
                  num_pka_recoil_points,&
                  num_pka_recoil_points_master,&
                  pka_recoil_energies_master,&
                  pka_recoil_energies)   


! define tdam contributions
   IF(do_tdam) THEN
    ! 29/6/2016 collapse pka rates onto global tdam energy spectrum using tdam energy spectrum for
    ! current channel
    ALLOCATE(tdam_pka_temp(MAX(num_pka_recoil_points_master,num_pka_recoil_points)))

     tdam_pka_temp(:)=pka(1,:)
     CALL collapse_fluxes(tdam_pka_temp, &
                  num_pka_recoil_points,&
                  num_pka_recoil_points_master,&
                  tdam_energies_master,&
                  tdam_energies)
    
   END IF


   
!secondary particle section
   SELECT CASE(mtd)
   CASE(5,22:25,29:30,35:36,45,107:109,112:114,117,800:849)
    IF(INDEX(pka_element,"alpha").NE.0) THEN
     ! contains alpha recoils

      !26/6/2013
      !alpha_sum_flag - if false then don't add 800:849
      ! if mtd==107 then set to false
      ! assumption - 107 would come first??
      IF(mtd==107) alpha_sum_flag=.false.
      IF((mtd.LT.800).OR.(mtd.GT.849)) THEN
       pka_sums(1,:)=pka_sums(1,:)+pka_temp(:)
       IF(do_tdam) THEN
         mtd_disp_sums(1)=mtd_disp_sums(1)+displacements
         pka_sums_tdam(1,:)=pka_sums_tdam(1,:)+tdam_pka_temp(:)
       END IF
      ELSE IF(alpha_sum_flag) THEN
       pka_sums(1,:)=pka_sums(1,:)+pka_temp(:)
       IF(do_tdam) THEN
         mtd_disp_sums(1)=mtd_disp_sums(1)+displacements
         pka_sums_tdam(1,:)=pka_sums_tdam(1,:)+tdam_pka_temp(:)
       END IF
      END IF
      
     END IF    

   CASE DEFAULT
   
   END SELECT
   SELECT CASE(mtd)
   CASE(5,28,41:45,103,111:112,115:116,600:649)
    IF(INDEX(pka_element,"proton").NE.0) THEN
     ! contains proton recoils !26/6/2013

      !26/6/2013
      !proton_sum_flag - if false then don't add 600:649
      ! if mtd==103 then set to false
      ! assumption - 103 would come first??
      IF(mtd==103) proton_sum_flag=.false.
      IF((mtd.LT.600).OR.(mtd.GT.649)) THEN
       pka_sums(5,:)=pka_sums(5,:)+pka_temp(:)
       IF(do_tdam) THEN
         mtd_disp_sums(5)=mtd_disp_sums(5)+displacements
         pka_sums_tdam(5,:)=pka_sums_tdam(1,:)+tdam_pka_temp(:)
       END IF
      ELSE IF(proton_sum_flag) THEN
       pka_sums(5,:)=pka_sums(5,:)+pka_temp(:)
       IF(do_tdam) THEN
         mtd_disp_sums(5)=mtd_disp_sums(5)+displacements
         pka_sums_tdam(5,:)=pka_sums_tdam(5,:)+tdam_pka_temp(:)
       END IF
      END IF     
     
    END IF     
   CASE DEFAULT
    ! do nothing - not a summable reaction
   END SELECT

   
! heavy recoils section   
   SELECT CASE(mtd)
   CASE(107,800:849)
    IF(INDEX(pka_element,"recoil").NE.0) THEN
     ! 26/6/2013 - total (n,a) recoil of daughter

      !26/6/2013
      !alpha_sum_flag - if false then don't add 800:849
      ! if mtd==107 then set to false
      ! assumption - 107 would come first??
      IF(mtd==107) alpha_sum_flag=.false.
      IF((mtd.LT.800).OR.(mtd.GT.849)) THEN
       pka_sums(4,:)=pka_sums(4,:)+pka_temp(:)
       IF(do_tdam) THEN
         mtd_disp_sums(4)=mtd_disp_sums(4)+displacements
         pka_sums_tdam(4,:)=pka_sums_tdam(4,:)+tdam_pka_temp(:)
       END IF
      ELSE IF(alpha_sum_flag) THEN
       pka_sums(4,:)=pka_sums(4,:)+pka_temp(:)
       IF(do_tdam) THEN
         mtd_disp_sums(4)=mtd_disp_sums(4)+displacements
         pka_sums_tdam(4,:)=pka_sums_tdam(4,:)+tdam_pka_temp(:)
       END IF
      END IF     
    END IF
   CASE(51:91) ! inelastic scattering
    IF(INDEX(pka_element,"recoil").NE.0) THEN

      pka_sums(2,:)=pka_sums(2,:)+pka_temp(:)
       IF(do_tdam) THEN
         mtd_disp_sums(2)=mtd_disp_sums(2)+displacements
         pka_sums_tdam(2,:)=pka_sums_tdam(2,:)+tdam_pka_temp(:)
       END IF
      ! 21/6/2013 - also sum with elastic part
      pka_sums(3,:)=pka_sums(3,:)+pka_temp(:)
       IF(do_tdam) THEN
         mtd_disp_sums(3)=mtd_disp_sums(3)+displacements
         pka_sums_tdam(3,:)=pka_sums_tdam(3,:)+tdam_pka_temp(:)
       END IF
    END IF
   CASE(2) !total scattering

      ! 21/6/2013 - sum with inelastic part
      pka_sums(3,:)=pka_sums(3,:)+pka_temp(:)
       IF(do_tdam) THEN
         mtd_disp_sums(3)=mtd_disp_sums(3)+displacements
         pka_sums_tdam(3,:)=pka_sums_tdam(3,:)+tdam_pka_temp(:)
       END IF
   CASE(103,600:649)
    IF(INDEX(pka_element,"recoil").NE.0) THEN
     ! total (n,p) recoils !26/6/2013

      !26/6/2013
      !proton_sum_flag - if false then don't add 600:649
      ! if mtd==103 then set to false
      ! assumption - 103 would come first??
      IF(mtd==103) proton_sum_flag=.false.
      IF((mtd.LT.600).OR.(mtd.GT.649)) THEN
       pka_sums(6,:)=pka_sums(6,:)+pka_temp(:)
       IF(do_tdam) THEN
         mtd_disp_sums(6)=mtd_disp_sums(6)+displacements
         pka_sums_tdam(6,:)=pka_sums_tdam(6,:)+tdam_pka_temp(:)
       END IF
      ELSE IF(proton_sum_flag) THEN
       pka_sums(6,:)=pka_sums(6,:)+pka_temp(:)
       IF(do_tdam) THEN
         mtd_disp_sums(6)=mtd_disp_sums(6)+displacements
         pka_sums_tdam(6,:)=pka_sums_tdam(6,:)+tdam_pka_temp(:)
       END IF
      END IF
    END IF     
   CASE DEFAULT
    ! do nothing - not a summable reaction
   END SELECT
   
   !print *,'h2',mtd,do_tdam
   IF(do_tdam) THEN
    !PRINT *,allocated(tdam_pka_temp),size(tdam_pka_temp)
    DEALLOCATE(tdam_pka_temp,STAT=deallocerr)
    !PRINT *,allocated(tdam_pka_temp),size(tdam_pka_temp),deallocerr 
   END IF  
   DEALLOCATE(pka_temp)
  !print *,'h3',mtd,do_tdam
  !STOP
  END SUBROUTINE sum_pkas