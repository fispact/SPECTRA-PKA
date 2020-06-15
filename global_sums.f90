!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk,                                                     *
!******************************************************************************
!end rubric

  SUBROUTINE add_to_globals(input_pka,recoil_points,tdam_vector,energy_vector)
   use globals
   IMPLICIT NONE
   INTEGER, INTENT (in) :: recoil_points
   REAL (KIND=DBL), INTENT(in) :: input_pka(recoil_points)
   REAL (KIND=DBL), INTENT(in) :: tdam_vector(recoil_points),energy_vector(recoil_points)
   REAL (KIND=DBL), ALLOCATABLE :: tdam_pka_temp(:),pka_temp(:)
   LOGICAL :: found
   INTEGER :: ii,jj !ii stores isotope position, jj the element position
   
   ii=1
   found=.false.
   DO WHILE ((ii.LE.number_global_recoils).AND.(.not.found))
    IF((daughter_num==global_daughter_nums(ii)).AND. &
    (TRIM(ADJUSTL(daughter_ele))==TRIM(ADJUSTL(global_daughter_eles(ii))))) THEN
      found=.true.
    ELSE
     ii=ii+1
    END IF
   END DO
   
   jj=1
   found=.false.
   DO WHILE ((jj.LE.number_global_recoil_elements).AND.(.not.found))
    IF((TRIM(ADJUSTL(daughter_ele))==TRIM(ADJUSTL(global_elements(jj))))) THEN
      found=.true.
    ELSE
     jj=jj+1
    END IF
   END DO
   
   IF(ii.GT.number_global_recoils) THEN
    number_global_recoils=ii
    IF(number_global_recoils.GT.max_global_recoils) THEN
     PRINT *,'max number of global nuclide recoils exceeded!'
     PRINT *,' increase max_global_recoils - STOP'
     WRITE(log_unit,*) 'max number of global nuclide recoils exceeded!'
     WRITE(log_unit,*) ' increase max_global_recoils - STOP'     
     STOP
    END IF
    global_daughter_eles(ii)=daughter_ele
    global_daughter_nums(ii)=daughter_num
    IF(do_outputs) WRITE(log_unit,*) 'number global nuclide recoils: ',number_global_recoils,TRIM(ADJUSTL(daughter_ele))
    
   END IF   
   
   
   IF(jj.GT.number_global_recoil_elements) THEN
    number_global_recoil_elements=jj
    IF(number_global_recoil_elements.GT.max_global_recoils) THEN
     WRITE(log_unit,*) 'max number of global recoil elements exceeded- STOP'
     STOP
    END IF
    global_elements(jj)=daughter_ele
    WRITE(log_unit,*) 'number global element recoils: ',number_global_recoils
   END IF
   

   IF(do_outputs) WRITE(log_unit,*) daughter_ele,daughter_num,'adding to globals'
    !24/4/2018 changes to allow for user output grid for sums and globals
    ! also removes error check for mis-matched grids - now we can just proceed and flag   
  IF(recoil_points.NE.global_num_pka_recoil_points_master) THEN
  ! 24/4/2018 allow for possibility of different grid and flag if do_user_output_energy_grid is false
   IF ((.NOT.do_user_output_energy_grid).OR.(do_user_output_energy_grid.AND.(user_grid_option==3))) THEN
      WRITE(log_unit,*) 'grid of recoil energies do not match global energy grid during summing'
      WRITE(log_unit,*) 'interpolation will be used to collapse the current grid onto the global one'
   END IF
  END IF
    
    ALLOCATE(pka_temp(MAX(global_num_pka_recoil_points_master,recoil_points)))
    pka_temp(1:recoil_points)=input_pka(:)*pka_ratios(filenum)
    IF(do_outputs) WRITE(log_unit,*) 'global sum allocate done',size(pka_recoil_energies), &
       size(global_pka_recoil_energies_master),recoil_points,global_num_pka_recoil_points_master, &
       size(pka_temp)
       
    CALL collapse_fluxes(pka_temp, &
                  recoil_points,&
                  global_num_pka_recoil_points_master,&
                  global_pka_recoil_energies_master,&
                  energy_vector) 
    IF(do_outputs) WRITE(log_unit,*) 'global sum flux collapse done'
   global_pka_sums(ii,:)=global_pka_sums(ii,:)+pka_temp(1:global_num_pka_recoil_points_master)
   global_pka_sums_element(jj,:)=global_pka_sums_element(jj,:)+pka_temp(1:global_num_pka_recoil_points_master)
   IF(do_tdam) THEN
    global_disp_sums(ii)=global_disp_sums(ii)+displacements*pka_ratios(filenum)
    global_disp_sums_element(jj)=global_disp_sums_element(jj)+displacements*pka_ratios(filenum)
    
    
    ! 29/6/2016 collapse pka rates onto global tdam energy spectrum using tdam energy spectrum for
    ! current channel
    !24/4/2018 - this was wrong before - the grids could be different now
    ALLOCATE(tdam_pka_temp(MAX(global_num_pka_recoil_points_master,recoil_points)))
    tdam_pka_temp(1:recoil_points)=input_pka(:)*pka_ratios(filenum)

    CALL collapse_fluxes(tdam_pka_temp, &
                  recoil_points,&
                  global_num_pka_recoil_points_master,&
                  global_tdam_energies_master,&
                  tdam_vector)


   global_pka_sums_tdam(ii,:)=global_pka_sums_tdam(ii,:)+&
        tdam_pka_temp(1:global_num_pka_recoil_points_master)
   global_pka_sums_element_tdam(jj,:)=global_pka_sums_element_tdam(jj,:)+&
        tdam_pka_temp(1:global_num_pka_recoil_points_master)
    
    
   END IF


   
   !20/5/2014 - add to total recoil spectra
   IF((do_exclude_light_from_total).AND.&
    ((TRIM(ADJUSTL(daughter_ele))=='He').OR.(TRIM(ADJUSTL(daughter_ele))=='H'))) THEN
    WRITE(log_unit,*) 'omitting gas particle from total pkas'
    !total_pka_sum(:)=total_pka_sum(:)+input_pka(:)*pka_ratios(filenum)
   ELSE IF((do_exclude_unknown_from_total).AND.&
    (TRIM(ADJUSTL(daughter_ele))=='unknown')) THEN
    WRITE(log_unit,*) 'omitting unknown daughter from total pkas'
   ELSE 
   
  !24/4/2018 - handling of user_grid_option 
  IF(do_user_output_energy_grid.AND.(user_grid_option==3)) THEN

    IF(do_tdam) DEALLOCATE(tdam_pka_temp)
    DEALLOCATE(pka_temp)
    ALLOCATE(pka_temp(MAX(totalglobal_num_pka_recoil_points_master,recoil_points)))
    pka_temp(1:recoil_points)=input_pka(:)*pka_ratios(filenum)
    IF(do_outputs) WRITE(log_unit,*) SUM(pka_temp(1:recoil_points)), &
      SUM(totalglobal_pka_recoil_energies_master),SUM(energy_vector)
    CALL collapse_fluxes(pka_temp, &
                  recoil_points,&
                  totalglobal_num_pka_recoil_points_master,&
                  totalglobal_pka_recoil_energies_master,&
                  energy_vector)  
IF(do_outputs) WRITE(log_unit,*) SUM(pka_temp(1:totalglobal_num_pka_recoil_points_master))                  
                  
    IF(do_tdam) THEN
     ! 24/4/2018 - assume tdam totalglobal matrix is same as totalglobal pka energy matrix
     ALLOCATE(tdam_pka_temp(MAX(totalglobal_num_pka_recoil_points_master,recoil_points)))
     tdam_pka_temp(1:recoil_points)=input_pka(:)*pka_ratios(filenum)

      CALL collapse_fluxes(tdam_pka_temp, &
                  recoil_points,&
                  totalglobal_num_pka_recoil_points_master,&
                  totalglobal_pka_recoil_energies_master,&
                  tdam_vector)    
    END IF
  END IF
                  
   
    ! 24/4/2018 - need to compute tdam_pka_temp, pka_temp, etc seperately
    total_pka_sum(:)=total_pka_sum(:)+pka_temp(1:totalglobal_num_pka_recoil_points_master)
    IF(do_tdam) THEN
        
      total_disp=total_disp+displacements*pka_ratios(filenum)
    ! 29/6/2016 collapse pka rates onto global tdam energy spectrum using tdam energy spectrum for
    ! current channel
     total_pka_sum_tdam=total_pka_sum_tdam+&
        tdam_pka_temp(1:totalglobal_num_pka_recoil_points_master)
    END IF !do_tdam
           

   END IF

  
  
  IF(do_tdam) DEALLOCATE(tdam_pka_temp)
  DEALLOCATE(pka_temp)
  
  
  END SUBROUTINE add_to_globals
 
 
 
   SUBROUTINE output_global_sums()
   use globals
   IMPLICIT NONE
   
   INTEGER, allocatable :: order(:),positions(:)
   LOGICAL :: found
   CHARACTER (LEN=100) :: fmt1='(ES14.6,2x,ES14.6,4(3x,ES11.4),3x,ES11.4,5x,ES11.4)'
   CHARACTER (LEN=100) :: fmt2='(1x,a31,a8,6x,a11,3x,a15,a11,a19,a11)' 
   CHARACTER (LEN=100) :: fmt3='(1x,a31,a8,6x,a11,3x,a15)'
   REAL (KIND=DBL) :: energy_min
    
    DO i=1,number_global_recoils
    IF(sum(global_pka_sums(i,1:global_num_pka_recoil_points_master)).NE.0) THEN
       WRITE(results_unit,*) '### index ',file_index,' ##### ( totals )'
       WRITE(number_string,'(I5)') global_daughter_nums(i)
       WRITE(results_unit,*) '#  recoil matrix of '//TRIM(ADJUSTL(global_daughter_eles(i))) &
        //TRIM(ADJUSTL(number_string))
    WRITE(index_summary,*) file_index,' recoil matrix of '//TRIM(ADJUSTL(global_daughter_eles(i))) &
        //TRIM(ADJUSTL(number_string))  
       
       
   ! WRITE(results_unit,*) '#DIFF SPEC AVG. PKA RECOIL DISTRIBUTIONS - per target atom'
   WRITE(results_unit,*) '#PKA RECOIL DISTRIBUTIONS - per target atom'       

       !23/4/2018 - new lower bin bound for first bin - 
       ! if using user grid then set to zero - otherwise half of lowest bin bound.
       IF(do_user_output_energy_grid) THEN
        energy_min=0_DBL
       ELSE
        energy_min=global_pka_recoil_energies_master(1)/2._DBL
       END IF
       
       IF (do_tdam) THEN
       
       WRITE(results_unit,fmt2) '#RECOIL energy (MeV low & high)','PKAs/s','norm_sum',&
        'cumulative_sum','tdam-pkas','disp_energy_(eV/s)','NRT_dpa/s'
       WRITE(results_unit,'(1x,a50)') '#(or T-dam energy low+high for tdam-pkas/disp/dpa)' 
       
       
       j=1
       IF((global_pka_sums(i,j).NE.0).or.(global_pka_sums_tdam(i,j).NE.0)) &
              WRITE(results_unit,fmt1) &
               energy_min,global_pka_recoil_energies_master(j), &
                 global_pka_sums(i,j), &
                   global_pka_sums(i,j)/SUM(global_pka_sums(i,1:global_num_pka_recoil_points_master)), &
                   SUM(global_pka_sums(i,1:j))/SUM(global_pka_sums(i,1:global_num_pka_recoil_points_master)), &
                   global_pka_sums_tdam(i,j), &
                   (global_pka_sums_tdam(i,j))*1e6_DBL*((global_pka_recoil_energies_master(j)+&
                   global_pka_recoil_energies_master(j)/2._DBL)/2._DBL), &
                   (global_pka_sums_tdam(i,j))*1e6_DBL*(0.8_DBL/(2._DBL*assumed_ed))*&
                    ((global_pka_recoil_energies_master(j)+&
                   global_pka_recoil_energies_master(j)/2._DBL)/2._DBL)
       DO j=2,global_num_pka_recoil_points_master
        IF((global_pka_sums(i,j).NE.0).or.(global_pka_sums_tdam(i,j).NE.0)) &
               WRITE(results_unit,fmt1) &
                   global_pka_recoil_energies_master(j-1),global_pka_recoil_energies_master(j), &
                 global_pka_sums(i,j), &
                   global_pka_sums(i,j)/SUM(global_pka_sums(i,1:global_num_pka_recoil_points_master)), &
                   SUM(global_pka_sums(i,1:j))/SUM(global_pka_sums(i,1:global_num_pka_recoil_points_master)), &
                   global_pka_sums_tdam(i,j), &
                   (global_pka_sums_tdam(i,j))*1e6_DBL*((global_pka_recoil_energies_master(j)+&
                   global_pka_recoil_energies_master(j-1))/2._DBL), &
                   (global_pka_sums_tdam(i,j))*1e6_DBL*(0.8_DBL/(2._DBL*assumed_ed))*&
                    ((global_pka_recoil_energies_master(j)+&
                   global_pka_recoil_energies_master(j-1))/2._DBL)
       END DO  
       
       ELSE
       
       WRITE(results_unit,fmt3) '#RECOIL energy (MeV low & high)','PKAs/s','norm_sum','cumulative sum'
       
       j=1
       IF(global_pka_sums(i,j).NE.0) WRITE(results_unit,'(5ES11.4)') &
               energy_min,global_pka_recoil_energies_master(j), &
                 global_pka_sums(i,j), &
                   global_pka_sums(i,j)/SUM(global_pka_sums(i,1:global_num_pka_recoil_points_master)), &
                   SUM(global_pka_sums(i,1:j))/SUM(global_pka_sums(i,1:global_num_pka_recoil_points_master))       
       DO j=2,global_num_pka_recoil_points_master
        IF(global_pka_sums(i,j).NE.0) WRITE(results_unit,'(5ES11.4)') &
                   global_pka_recoil_energies_master(j-1),global_pka_recoil_energies_master(j), &
                 global_pka_sums(i,j), &
                   global_pka_sums(i,j)/SUM(global_pka_sums(i,1:global_num_pka_recoil_points_master)), &
                   SUM(global_pka_sums(i,1:j))/SUM(global_pka_sums(i,1:global_num_pka_recoil_points_master))
       END DO  
       END IF !do_tdam
           
   
       ! compute average pka energy
       pka_ave=0.5_DBL*(energy_min+global_pka_recoil_energies_master(1))*global_pka_sums(i,1)
       DO j=2,global_num_pka_recoil_points_master
        pka_ave=pka_ave+0.5_DBL*(global_pka_recoil_energies_master(j)+&
                                 global_pka_recoil_energies_master(j-1))*global_pka_sums(i,j)
       END DO
       pka_ave=pka_ave*1000000._DBL
       
       WRITE(results_unit,'(a1,20x,A,ES11.4,A5)') '#','AVERAGE PKA ENERGY = ',pka_ave/SUM(global_pka_sums(i,:)),' (eV)'
   IF (do_tdam) then
   
    WRITE(results_unit,'(a1,20x,A,ES11.4)') '#',&
             'displacement energy eV/s = ',global_disp_sums(i)*1e6_DBL
    WRITE(results_unit,'(a1,20x,A,ES11.4,A7,F4.1,a)') '#','equivalent NRT dpa/s = ',&
              0.8_DBL*global_disp_sums(i)/(2._DBL*assumed_ed*1e-6_DBL),&
            ' (E_d=',assumed_ed,' eV)'    
   

   END IF          
       WRITE(results_unit,*)
       WRITE(results_unit,*)
       file_index=file_index+1
    END IF !non zero check
   END DO ! loop over global recoils
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   !10/10/2013 - and per element sums
   
   !define order based on position in master_elements
   
   ALLOCATE(order(number_global_recoil_elements),positions(number_global_recoil_elements))
   order=0
   
   ! define element position in list of elements
 DO i=1,number_global_recoil_elements   
   positions(i)=1
   found=.false.
   DO WHILE ((.not.found).AND.(positions(i).LE.num_master_elements))
    IF(TRIM(ADJUSTL(global_elements(i)))==TRIM(ADJUSTL(master_elements(positions(i)))) ) THEN
     found=.true.
    ELSE
     positions(i)=positions(i)+1
    END IF
   END DO 
 END DO
   
   DO i=1,number_global_recoil_elements
      j=minloc(positions(1:number_global_recoil_elements),1)
      order(i)=j
      positions(j)=num_master_elements*10
      !PRINT *,order(i)
   END DO
     
   
   DO k=1,number_global_recoil_elements
    i=order(k)
    IF(sum(global_pka_sums_element(i,1:global_num_pka_recoil_points_master)).NE.0) THEN
       WRITE(results_unit,*) '### index ',file_index,' ##### ( totals )'
       WRITE(results_unit,*) '#  elemental recoil matrix of '//TRIM(ADJUSTL(global_elements(i))) 
    WRITE(index_summary,*) file_index,' elemental recoil matrix of '//TRIM(ADJUSTL(global_elements(i)))  
       
   ! WRITE(results_unit,*) '#DIFF SPEC AVG. PKA RECOIL DISTRIBUTIONS - per target atom'
   WRITE(results_unit,*) '#PKA RECOIL DISTRIBUTIONS - per target atom'  
       
       
       
     IF(do_tdam) THEN

       WRITE(results_unit,fmt2) '#RECOIL energy (MeV low & high)','PKAs/s','norm_sum',&
        'cumulative_sum','tdam-pkas','disp_energy_(eV/s)','NRT_dpa/s'
       WRITE(results_unit,'(1x,a50)') '#(or T-dam energy low+high for tdam-pkas/disp/dpa)' 

        
       j=1
       IF((global_pka_sums_element(i,j).NE.0).or.(global_pka_sums_element_tdam(i,j).NE.0)) &
            WRITE(results_unit,fmt1) &
               energy_min,global_pka_recoil_energies_master(j), &
                 global_pka_sums_element(i,j), &
                   global_pka_sums_element(i,j)/SUM(global_pka_sums_element(i,1:global_num_pka_recoil_points_master)), &
              SUM(global_pka_sums_element(i,1:j))/SUM(global_pka_sums_element(i,1:global_num_pka_recoil_points_master)), &
              global_pka_sums_element_tdam(i,j), &
                   (global_pka_sums_element_tdam(i,j))*1e6_DBL*((global_pka_recoil_energies_master(j)+&
                   global_pka_recoil_energies_master(j)/2._DBL)/2._DBL), &
                   (global_pka_sums_element_tdam(i,j))*1e6_DBL*(0.8_DBL/(2._DBL*assumed_ed))*&
                    ((global_pka_recoil_energies_master(j)+&
                   global_pka_recoil_energies_master(j)/2._DBL)/2._DBL)
       DO j=2,global_num_pka_recoil_points_master
        IF((global_pka_sums_element(i,j).NE.0).or.(global_pka_sums_element_tdam(i,j).NE.0)) &
              WRITE(results_unit,fmt1) &
                   global_pka_recoil_energies_master(j-1),global_pka_recoil_energies_master(j), &
                 global_pka_sums_element(i,j), &
                   global_pka_sums_element(i,j)/SUM(global_pka_sums_element(i,1:global_num_pka_recoil_points_master)), &
      SUM(global_pka_sums_element(i,1:j))/SUM(global_pka_sums_element(i,1:global_num_pka_recoil_points_master)), &
              global_pka_sums_element_tdam(i,j), &
                   (global_pka_sums_element_tdam(i,j))*1e6_DBL*((global_pka_recoil_energies_master(j)+&
                   global_pka_recoil_energies_master(j-1))/2._DBL), &
                   (global_pka_sums_element_tdam(i,j))*1e6_DBL*(0.8_DBL/(2._DBL*assumed_ed))*&
                    ((global_pka_recoil_energies_master(j)+&
                   global_pka_recoil_energies_master(j-1))/2._DBL) 
       END DO 
      
     ELSE
       WRITE(results_unit,fmt3) '#RECOIL energy (MeV low & high)','PKAs/s','norm_sum','cumulative sum'
        
       j=1
       IF(global_pka_sums_element(i,j).NE.0) WRITE(results_unit,'(5ES11.4)') &
               energy_min,global_pka_recoil_energies_master(j), &
                 global_pka_sums_element(i,j), &
                   global_pka_sums_element(i,j)/SUM(global_pka_sums_element(i,1:global_num_pka_recoil_points_master)), &
                  SUM(global_pka_sums_element(i,1:j))/SUM(global_pka_sums_element(i,1:global_num_pka_recoil_points_master))       
       DO j=2,global_num_pka_recoil_points_master
        IF(global_pka_sums_element(i,j).NE.0) WRITE(results_unit,'(5ES11.4)') &
                   global_pka_recoil_energies_master(j-1),global_pka_recoil_energies_master(j), &
                 global_pka_sums_element(i,j), &
                   global_pka_sums_element(i,j)/SUM(global_pka_sums_element(i,1:global_num_pka_recoil_points_master)), &
                  SUM(global_pka_sums_element(i,1:j))/SUM(global_pka_sums_element(i,1:global_num_pka_recoil_points_master)) 
       END DO 
     END IF
  
   
       ! compute average pka energy
       pka_ave=0.5_DBL*(energy_min+global_pka_recoil_energies_master(1))*global_pka_sums_element(i,1)
       DO j=2,global_num_pka_recoil_points_master
        pka_ave=pka_ave+0.5_DBL*(global_pka_recoil_energies_master(j)+&
                                 global_pka_recoil_energies_master(j-1))*global_pka_sums_element(i,j)
       END DO
       pka_ave=pka_ave*1000000._DBL
       
       WRITE(results_unit,'(a1,20x,A,ES11.4,A5)') '#','AVERAGE PKA ENERGY = ',pka_ave/SUM(global_pka_sums_element(i,:)),' (eV)'
   IF (do_tdam) then
   
    WRITE(results_unit,'(a1,20x,A,ES11.4)') '#',&
           'displacement energy eV/s = ',global_disp_sums_element(i)*1e6_DBL
    WRITE(results_unit,'(a1,20x,A,ES11.4,A7,F4.1,a)') '#','equivalent NRT dpa/s = ',&
              0.8_DBL*global_disp_sums_element(i)/(2._DBL*assumed_ed*1e-6_DBL),&
            ' (E_d=',assumed_ed,' eV)'    
   

   END IF       
       WRITE(results_unit,*)
       WRITE(results_unit,*)
       file_index=file_index+1
    END IF !non zero check
   END DO ! loop over global recoils 
   
   DEALLOCATE(order,positions)
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   ! 20/5/2014
   !total recoil spectrum
   
       !23/4/2018 - new lower bin bound for first bin - 
       ! if using user grid then set to zero - otherwise half of lowest bin bound.
       IF(do_user_output_energy_grid) THEN
        energy_min=0_DBL
       ELSE
        energy_min=totalglobal_pka_recoil_energies_master(1)/2._DBL
       END IF   
    IF(do_outputs) WRITE(log_unit,*) 'total sum',SIZE(total_pka_sum),totalglobal_num_pka_recoil_points_master, &
       sum(total_pka_sum(1:totalglobal_num_pka_recoil_points_master))
    !IF(sum(total_pka_sum(1:totalglobal_num_pka_recoil_points_master)).NE.0) THEN   
       WRITE(results_unit,*) '### index ',file_index,' ##### ( totals )'
       IF((do_exclude_light_from_total).AND.(do_exclude_unknown_from_total)) THEN
        WRITE(results_unit,*) '#  total recoil matrix (He+H+unknowns excluded)'
        WRITE(index_summary,*) file_index,' total recoil matrix (He+H+unknowns excluded)'       
       ELSE IF (do_exclude_light_from_total) THEN
        WRITE(results_unit,*) '#  total recoil matrix (He+H excluded)'
        WRITE(index_summary,*) file_index,' total recoil matrix (He+H excluded)'  
       ELSE IF (do_exclude_unknown_from_total) THEN
        WRITE(results_unit,*) '#  total recoil matrix (unknowns excluded)'
        WRITE(index_summary,*) file_index,' total recoil matrix (unknowns excluded)'        
       ELSE
        WRITE(results_unit,*) '#  total recoil matrix '
        WRITE(index_summary,*) file_index,' total recoil matrix'  
       END IF
       
   ! WRITE(results_unit,*) '#DIFF SPEC AVG. PKA RECOIL DISTRIBUTIONS - per target atom'
   WRITE(results_unit,*) '#PKA RECOIL DISTRIBUTIONS - per target atom'       
       
      IF(do_tdam) THEN

       WRITE(results_unit,fmt2) '#RECOIL energy (MeV low & high)','PKAs/s','norm_sum',&
        'cumulative_sum','tdam-pkas','disp_energy_(eV/s)','NRT_dpa/s'
       WRITE(results_unit,'(1x,a50)') '#(or T-dam energy low+high for tdam-pkas/disp/dpa)' 
        
       j=1
       IF((total_pka_sum(j).NE.0).or.(total_pka_sum_tdam(j).NE.0)) &
            WRITE(results_unit,fmt1) &
               energy_min,totalglobal_pka_recoil_energies_master(j), &
                 total_pka_sum(j), &
                   total_pka_sum(j)/SUM(total_pka_sum(1:totalglobal_num_pka_recoil_points_master)) , &
                   SUM(total_pka_sum(1:j))/SUM(total_pka_sum(1:totalglobal_num_pka_recoil_points_master)), &
                   total_pka_sum_tdam(j), &
                   (total_pka_sum_tdam(j))*1e6_DBL*((totalglobal_pka_recoil_energies_master(j)+&
                   totalglobal_pka_recoil_energies_master(j)/2._DBL)/2._DBL), &
                   (total_pka_sum_tdam(j))*1e6_DBL*(0.8_DBL/(2._DBL*assumed_ed))*&
                    ((totalglobal_pka_recoil_energies_master(j)+&
                   totalglobal_pka_recoil_energies_master(j)/2._DBL)/2._DBL)
       DO j=2,totalglobal_num_pka_recoil_points_master
        IF((total_pka_sum(j).NE.0).or.(total_pka_sum_tdam(j).NE.0)) &
            WRITE(results_unit,fmt1) &
                   totalglobal_pka_recoil_energies_master(j-1),totalglobal_pka_recoil_energies_master(j), &
                 total_pka_sum(j), &
                   total_pka_sum(j)/SUM(total_pka_sum(1:totalglobal_num_pka_recoil_points_master)), &
                   SUM(total_pka_sum(1:j))/SUM(total_pka_sum(1:totalglobal_num_pka_recoil_points_master)), &
                   total_pka_sum_tdam(j), &
                   (total_pka_sum_tdam(j))*1e6_DBL*((totalglobal_pka_recoil_energies_master(j)+&
                   totalglobal_pka_recoil_energies_master(j-1))/2._DBL), &
                   (total_pka_sum_tdam(j))*1e6_DBL*(0.8_DBL/(2._DBL*assumed_ed))*&
                    ((totalglobal_pka_recoil_energies_master(j)+&
                   totalglobal_pka_recoil_energies_master(j-1))/2._DBL)
                   
       END DO 

      ELSE
       WRITE(results_unit,fmt3) '#RECOIL energy (MeV low & high)','PKAs/s','norm_sum','cumulative sum'
        
       j=1
       IF(total_pka_sum(j).NE.0) WRITE(results_unit,'(5ES11.4)') &
               energy_min,totalglobal_pka_recoil_energies_master(j), &
                 total_pka_sum(j), &
                   total_pka_sum(j)/SUM(total_pka_sum(1:totalglobal_num_pka_recoil_points_master)) , &
                   SUM(total_pka_sum(1:j))/SUM(total_pka_sum(1:totalglobal_num_pka_recoil_points_master))
       DO j=2,totalglobal_num_pka_recoil_points_master
        IF(total_pka_sum(j).NE.0) WRITE(results_unit,'(5ES11.4)') &
                   totalglobal_pka_recoil_energies_master(j-1),totalglobal_pka_recoil_energies_master(j), &
                 total_pka_sum(j), &
                   total_pka_sum(j)/SUM(total_pka_sum(1:totalglobal_num_pka_recoil_points_master)), &
                   SUM(total_pka_sum(1:j))/SUM(total_pka_sum(1:totalglobal_num_pka_recoil_points_master))
                   
       END DO 
      END IF
  
   
       ! compute average pka energy
       pka_ave=0.5_DBL*(energy_min+totalglobal_pka_recoil_energies_master(1))*total_pka_sum(1)
       DO j=2,totalglobal_num_pka_recoil_points_master
        pka_ave=pka_ave+0.5_DBL*(totalglobal_pka_recoil_energies_master(j)+&
                                 totalglobal_pka_recoil_energies_master(j-1))*total_pka_sum(j)
       END DO
       pka_ave=pka_ave*1000000._DBL
       
       WRITE(results_unit,'(a1,20x,A,ES11.4,A5)') '#','AVERAGE PKA ENERGY = ',pka_ave/SUM(total_pka_sum(:)),' (eV)'
   IF (do_tdam) then
    WRITE(results_unit,'(a1,20x,A,ES11.4)') '#',&
               'displacement energy eV/s = ',total_disp*1e6_DBL
    WRITE(results_unit,'(a1,20x,A,ES11.4,A7,F4.1,a)') '#','equivalent NRT dpa/s = ',&
              0.8_DBL*total_disp/(2._DBL*assumed_ed*1e-6_DBL),&
            ' (E_d=',assumed_ed,' eV)'   
   END IF         
       WRITE(results_unit,*)
       WRITE(results_unit,*)
       file_index=file_index+1
    !END IF !non zero check 
    
    

   
  END SUBROUTINE output_global_sums
