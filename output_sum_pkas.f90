!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk,                                                     *
!******************************************************************************
!end rubric

  SUBROUTINE output_sum_pkas()
  use globals
  IMPLICIT NONE
  CHARACTER (LEN=20) :: array(10)=(/"total alpha         ", &
                                    "(n,inelastic) recoil", &
                                    "scatter recoil      ", &
                                    "total (z,a) recoil  ", &
                                    "total proton        ", &
                                    "total (z,p) recoil  ", &
                                    "alpha               ", &
                                    "alpha               ", &
                                    "alpha               ", &
                                    "alpha               "/)
 INTEGER :: numarray(10)=(/107,2,2,107,103,103,0,0,0,0/)  
 LOGICAL :: logarray(10)=(/.false.,.true.,.true.,.true.,.false., &
                           .true.,.false.,.false.,.false.,.false./)
   CHARACTER (LEN=100) :: fmt1='(4ES14.6,2x,ES14.6,3x,ES11.4,5x,ES11.4)'
   CHARACTER (LEN=100) :: fmt2='(1x,a31,a8,4x,a11,5x,a11,a19,a11)' 
   CHARACTER (LEN=100) :: fmt3='(1x,a31,a8,a11)'  
  REAL (KIND=DBL) :: energy_min   


       !23/4/2018 - new lower bin bound for first bin - 
       ! if using user grid then set to zero - otherwise half of lowest bin bound.
       IF(do_user_output_energy_grid) THEN
        energy_min=0_DBL
       ELSE
        energy_min=pka_recoil_energies_master(1)/2._DBL
       END IF

  
  WRITE(number_string2,'(I5)') parent_num(filenum)
  
  DO j=1,6
   IF(j==1) pka_element='alpha'
   IF(j==5) pka_element='proton'
   mtd=numarray(j)
   CALL define_daughter(logarray(j))
   IF(sum(pka_sums(j,1:num_pka_recoil_points_master)).NE.0) THEN
      WRITE(results_unit,*) '### index ',file_index,' ##### ',TRIM(ADJUSTL(pka_filename(filenum)))
      
      WRITE(number_string,'(I5)') daughter_num
      
      WRITE(results_unit,*) '#  '//TRIM(ADJUSTL(array(j)))//' matrix [ '//TRIM(ADJUSTL(daughter_ele)) &
       //TRIM(ADJUSTL(number_string))//' ] from [ '//TRIM(ADJUSTL(parent_ele(filenum))) &
       //TRIM(ADJUSTL(number_string2)) &
       //' ]'   
    WRITE(index_summary,*) file_index,' '//TRIM(ADJUSTL(array(j)))//' matrix [ '//TRIM(ADJUSTL(daughter_ele)) &
       //TRIM(ADJUSTL(number_string))//' ] from [ '//TRIM(ADJUSTL(parent_ele(filenum))) &
       //TRIM(ADJUSTL(number_string2)) &
       //' ]'//TRIM(ADJUSTL(pka_filename(filenum)))  
      
   ! WRITE(results_unit,*) '#DIFF SPEC AVG. PKA RECOIL DISTRIBUTIONS'
   WRITE(results_unit,*) '#PKA RECOIL DISTRIBUTIONS'
      
      
      IF (do_tdam) THEN
       WRITE(results_unit,fmt2) '#RECOIL energy (MeV low & high)','PKAs/s','norm_sum',&
        'tdam-pkas','disp_energy_(eV/s)','NRT_dpa/s'
       WRITE(results_unit,'(1x,a30)') '#(or T_dam energy low+high for' 
       WRITE(results_unit,'(1x,a30)') '#          tdam-pkas/disp/dpa)'
       i=1
       IF(pka_sums(j,i).NE.0) WRITE(results_unit,fmt1) &
               energy_min,pka_recoil_energies_master(i), &
                 pka_sums(j,i),pka_sums(j,i)/SUM(pka_sums(j,1:num_pka_recoil_points_master)), &
                 pka_sums_tdam(j,i), &
                   (pka_sums_tdam(j,i))*1e6_DBL*((pka_recoil_energies_master(i)+&
                   pka_recoil_energies_master(i)/2._DBL)/2._DBL), &
                   (pka_sums_tdam(j,i))*1e6_DBL*(0.8_DBL/(2._DBL*assumed_ed))*&
                    ((pka_recoil_energies_master(i)+&
                   pka_recoil_energies_master(i)/2._DBL)/2._DBL)       
       DO i=2,num_pka_recoil_points_master
        IF(pka_sums(j,i).NE.0) WRITE(results_unit,fmt1) &
                  pka_recoil_energies_master(i-1),pka_recoil_energies_master(i), &
                 pka_sums(j,i),pka_sums(j,i)/SUM(pka_sums(j,1:num_pka_recoil_points_master)), &
                 pka_sums_tdam(j,i), &
                   (pka_sums_tdam(j,i))*1e6_DBL*((pka_recoil_energies_master(i)+&
                   pka_recoil_energies_master(i-1))/2._DBL), &
                   (pka_sums_tdam(j,i))*1e6_DBL*(0.8_DBL/(2._DBL*assumed_ed))*&
                    ((pka_recoil_energies_master(i)+&
                   pka_recoil_energies_master(i-1))/2._DBL)
       END DO       
      
      ELSE
       WRITE(results_unit,fmt3) '#RECOIL energy (MeV low & high)','PKAs/s','norm_sum'
       i=1
       IF(pka_sums(j,i).NE.0) WRITE(results_unit,'(4ES11.4)') &
               energy_min,pka_recoil_energies_master(i), &
                 pka_sums(j,i),pka_sums(j,i)/SUM(pka_sums(j,1:num_pka_recoil_points_master))       
       DO i=2,num_pka_recoil_points_master
        IF(pka_sums(j,i).NE.0) WRITE(results_unit,'(4ES11.4)') &
                  pka_recoil_energies_master(i-1),pka_recoil_energies_master(i), &
                 pka_sums(j,i),pka_sums(j,i)/SUM(pka_sums(j,1:num_pka_recoil_points_master))
       END DO   
      END IF
  
      ! compute average pka energy
      pka_ave=0.5_DBL*(pka_recoil_energies_master(1)+energy_min)*pka_sums(j,1)
      DO i=2,num_pka_recoil_points_master  ! 24/4/2018 wrongly said i=1, before
       pka_ave=pka_ave+0.5_DBL*(pka_recoil_energies_master(i)+pka_recoil_energies_master(i-1))*pka_sums(j,i)
      END DO
      pka_ave=pka_ave*1000000._DBL
      
      WRITE(results_unit,'(a1,20x,A,ES11.4,A5)') '#','AVERAGE PKA ENERGY = ',pka_ave/SUM(pka_sums(j,:)),' (eV)'
   IF (do_tdam) then
   
    WRITE(results_unit,'(a1,20x,A,ES11.4)') '#',&
               'displacement energy eV/s = ',mtd_disp_sums(j)*1e6_DBL 
    WRITE(results_unit,'(a1,20x,A,ES11.4,A7,F4.1,a)') '#','equivalent NRT dpa/s = ',&
              0.8_DBL*mtd_disp_sums(j)/(2._DBL*assumed_ed*1e-6_DBL),&
            ' (E_d=',assumed_ed,' eV)'   
   

   END IF      
    WRITE(results_unit,*)
    WRITE(results_unit,*)
    file_index=file_index+1
   END IF !non zero check  

  
END DO !j




 

   
 
   




  
  END SUBROUTINE output_sum_pkas
