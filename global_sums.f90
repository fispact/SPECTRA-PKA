!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk,                                                     *
!******************************************************************************
!end rubric

  SUBROUTINE add_to_globals(input_pka,recoil_points,tdam_vector)
   use globals
   IMPLICIT NONE
   INTEGER, INTENT (in) :: recoil_points
   REAL (KIND=DBL), INTENT(in) :: input_pka(global_num_pka_recoil_points_master)
   REAL (KIND=DBL), INTENT(in) :: tdam_vector(global_num_pka_recoil_points_master)
   REAL (KIND=DBL), ALLOCATABLE :: tdam_pka_temp(:)
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
     STOP
    END IF
    global_daughter_eles(ii)=daughter_ele
    global_daughter_nums(ii)=daughter_num
    PRINT *,'number global nuclide recoils: ',number_global_recoils
   END IF   
   
   
   IF(jj.GT.number_global_recoil_elements) THEN
    number_global_recoil_elements=jj
    IF(number_global_recoil_elements.GT.max_global_recoils) THEN
     PRINT *,'max number of global recoil elements exceeded- STOP'
     STOP
    END IF
    global_elements(jj)=daughter_ele
    PRINT *,'number global element recoils: ',number_global_recoils
   END IF
   

   

   !PRINT *,num_pka_recoil_points,global_num_pka_recoil_points_master
   !PRINT *,'one'
  IF(recoil_points.NE.global_num_pka_recoil_points_master) THEN
      PRINT *,'grid of recoil energies do not match while global summing over mt numbers'
      io_quit=1   
  ELSE

   global_pka_sums(ii,:)=global_pka_sums(ii,:)+input_pka(:)*pka_ratios(filenum)
   global_pka_sums_element(jj,:)=global_pka_sums_element(jj,:)+input_pka(:)*pka_ratios(filenum)
   IF(do_tdam) THEN
    global_disp_sums(ii)=global_disp_sums(ii)+displacements*pka_ratios(filenum)
    global_disp_sums_element(jj)=global_disp_sums_element(jj)+displacements*pka_ratios(filenum)
    
    
    ! 29/6/2016 collapse pka rates onto global tdam energy spectrum using tdam energy spectrum for
    ! current channel
    ALLOCATE(tdam_pka_temp(global_num_pka_recoil_points_master))
    tdam_pka_temp=input_pka(:)*pka_ratios(filenum)
    !PRINT *,'in',SUM(tdam_pka_temp)
    !print *,tdam_vector
    !print *,tdam_pka_temp
    !print *,global_tdam_energies_master
    CALL collapse_fluxes(tdam_pka_temp, &
                  global_num_pka_recoil_points_master,&
                  global_num_pka_recoil_points_master,&
                  global_tdam_energies_master,&
                  tdam_vector)
   !PRINT *,'out',SUM(tdam_pka_temp)  
   !print *,tdam_pka_temp
   global_pka_sums_tdam(ii,:)=global_pka_sums_tdam(ii,:)+&
        tdam_pka_temp(:)
   global_pka_sums_element_tdam(jj,:)=global_pka_sums_element_tdam(jj,:)+&
        tdam_pka_temp(:)
    
    
   END IF
   
   !20/5/2014 - add to total recoil spectra
   IF((do_exclude_light_from_total).AND.&
    ((TRIM(ADJUSTL(daughter_ele))=='He').OR.(TRIM(ADJUSTL(daughter_ele))=='H'))) THEN
    PRINT *,'omitting gas particle from total pkas'
    !total_pka_sum(:)=total_pka_sum(:)+input_pka(:)*pka_ratios(filenum)
   ELSE IF((do_exclude_unknown_from_total).AND.&
    (TRIM(ADJUSTL(daughter_ele))=='unknown')) THEN
    PRINT *,'omitting unknown daughter from total pkas'
   ELSE 
    total_pka_sum(:)=total_pka_sum(:)+input_pka(:)*pka_ratios(filenum)
    IF(do_tdam) THEN
      total_disp=total_disp+displacements*pka_ratios(filenum)
    ! 29/6/2016 collapse pka rates onto global tdam energy spectrum using tdam energy spectrum for
    ! current channel
     total_pka_sum_tdam=total_pka_sum_tdam+&
        tdam_pka_temp
           
    END IF
   END IF
  END IF 
  
  
  IF(do_tdam) DEALLOCATE(tdam_pka_temp)

  END SUBROUTINE add_to_globals
 
 
 
   SUBROUTINE output_global_sums()
   use globals
   IMPLICIT NONE
   
   INTEGER, allocatable :: order(:),positions(:)
   LOGICAL :: found
   CHARACTER (LEN=100) :: fmt1='(4ES11.4,2x,ES11.4,3x,ES11.4,3x,ES11.4,5x,ES11.4)'
   CHARACTER (LEN=100) :: fmt2='(1x,a24,a8,a11,a15,1x,a11,a19,a11)' 
   CHARACTER (LEN=100) :: fmt3='(1x,a24,a8,a11,a15)'
    
    DO i=1,number_global_recoils
    IF(sum(global_pka_sums(i,1:global_num_pka_recoil_points_master)).NE.0) THEN
       WRITE(results_unit,*) '### index ',file_index,' ##### ( totals )'
       WRITE(number_string,'(I5)') global_daughter_nums(i)
       WRITE(results_unit,*) '#  recoil matrix of '//TRIM(ADJUSTL(global_daughter_eles(i))) &
        //TRIM(ADJUSTL(number_string))
    WRITE(index_summary,*) file_index,' recoil matrix of '//TRIM(ADJUSTL(global_daughter_eles(i))) &
        //TRIM(ADJUSTL(number_string))  
       
       WRITE(results_unit,*) '#DIFF SPEC AVG. PKA RECOIL DISTRIBUTIONS - per target atom'
       
       


       
       IF (do_tdam) THEN
       
       WRITE(results_unit,fmt2) '#ERECOIL(MeV_low_&_high))','PKAs/s','norm_sum',&
        'cumulative_sum','tdam-pkas','disp_energy_(eV/s)','NRT_dpa/s'
       
       
       j=1
       IF((global_pka_sums(i,j).NE.0).or.(global_pka_sums_tdam(i,j).NE.0)) &
              WRITE(results_unit,fmt1) &
               global_pka_recoil_energies_master(j)/2._DBL,global_pka_recoil_energies_master(j), &
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
       
       WRITE(results_unit,fmt3) '#ERECOIL(MeV low & high))','PKAs/s','norm_sum','cumulative sum'
       
       j=1
       IF(global_pka_sums(i,j).NE.0) WRITE(results_unit,'(5ES11.4)') &
               global_pka_recoil_energies_master(j)/2._DBL,global_pka_recoil_energies_master(j), &
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
       pka_ave=0.5_DBL*global_pka_recoil_energies_master(1)*global_pka_sums(i,1)
       DO j=1,global_num_pka_recoil_points_master
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
   !PRINT *,positions(i),global_elements(i)
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
       
       WRITE(results_unit,*) '#DIFF SPEC AVG. PKA RECOIL DISTRIBUTIONS - per target atom'
       
       
       
     IF(do_tdam) THEN

       WRITE(results_unit,fmt2) '#ERECOIL(MeV_low_&_high))','PKAs/s','norm_sum',&
        'cumulative_sum','tdam-pkas','disp_energy_(eV/s)','NRT_dpa/s'
        
       j=1
       IF((global_pka_sums_element(i,j).NE.0).or.(global_pka_sums_element_tdam(i,j).NE.0)) &
            WRITE(results_unit,fmt1) &
               global_pka_recoil_energies_master(j)/2._DBL,global_pka_recoil_energies_master(j), &
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
       WRITE(results_unit,fmt3) '#ERECOIL(MeV low & high)','PKAs/s','norm_sum','cumulative sum'
        
       j=1
       IF(global_pka_sums_element(i,j).NE.0) WRITE(results_unit,'(5ES11.4)') &
               global_pka_recoil_energies_master(j)/2._DBL,global_pka_recoil_energies_master(j), &
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
       pka_ave=0.5_DBL*global_pka_recoil_energies_master(1)*global_pka_sums_element(i,1)
       DO j=1,global_num_pka_recoil_points_master
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
    IF(sum(total_pka_sum(1:global_num_pka_recoil_points_master)).NE.0) THEN   
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
       
       WRITE(results_unit,*) '#DIFF SPEC AVG. PKA RECOIL DISTRIBUTIONS - per target atom'
       
       
      IF(do_tdam) THEN

       WRITE(results_unit,fmt2) '#ERECOIL(MeV_low_&_high))','PKAs/s','norm_sum',&
        'cumulative_sum','tdam-pkas','disp_energy_(eV/s)','NRT_dpa/s'
        
       j=1
       IF((total_pka_sum(j).NE.0).or.(total_pka_sum_tdam(j).NE.0)) &
            WRITE(results_unit,fmt1) &
               global_pka_recoil_energies_master(j)/2._DBL,global_pka_recoil_energies_master(j), &
                 total_pka_sum(j), &
                   total_pka_sum(j)/SUM(total_pka_sum(1:global_num_pka_recoil_points_master)) , &
                   SUM(total_pka_sum(1:j))/SUM(total_pka_sum(1:global_num_pka_recoil_points_master)), &
                   total_pka_sum_tdam(j), &
                   (total_pka_sum_tdam(j))*1e6_DBL*((global_pka_recoil_energies_master(j)+&
                   global_pka_recoil_energies_master(j)/2._DBL)/2._DBL), &
                   (total_pka_sum_tdam(j))*1e6_DBL*(0.8_DBL/(2._DBL*assumed_ed))*&
                    ((global_pka_recoil_energies_master(j)+&
                   global_pka_recoil_energies_master(j)/2._DBL)/2._DBL)
       DO j=2,global_num_pka_recoil_points_master
        IF((total_pka_sum(j).NE.0).or.(total_pka_sum_tdam(j).NE.0)) &
            WRITE(results_unit,fmt1) &
                   global_pka_recoil_energies_master(j-1),global_pka_recoil_energies_master(j), &
                 total_pka_sum(j), &
                   total_pka_sum(j)/SUM(total_pka_sum(1:global_num_pka_recoil_points_master)), &
                   SUM(total_pka_sum(1:j))/SUM(total_pka_sum(1:global_num_pka_recoil_points_master)), &
                   total_pka_sum_tdam(j), &
                   (total_pka_sum_tdam(j))*1e6_DBL*((global_pka_recoil_energies_master(j)+&
                   global_pka_recoil_energies_master(j-1))/2._DBL), &
                   (total_pka_sum_tdam(j))*1e6_DBL*(0.8_DBL/(2._DBL*assumed_ed))*&
                    ((global_pka_recoil_energies_master(j)+&
                   global_pka_recoil_energies_master(j-1))/2._DBL)
                   
       END DO 

      ELSE
       WRITE(results_unit,fmt3) '#ERECOIL(MeV low & high)','PKAs/s','norm_sum','cumulative sum'
        
       j=1
       IF(total_pka_sum(j).NE.0) WRITE(results_unit,'(5ES11.4)') &
               global_pka_recoil_energies_master(j)/2._DBL,global_pka_recoil_energies_master(j), &
                 total_pka_sum(j), &
                   total_pka_sum(j)/SUM(total_pka_sum(1:global_num_pka_recoil_points_master)) , &
                   SUM(total_pka_sum(1:j))/SUM(total_pka_sum(1:global_num_pka_recoil_points_master))
       DO j=2,global_num_pka_recoil_points_master
        IF(total_pka_sum(j).NE.0) WRITE(results_unit,'(5ES11.4)') &
                   global_pka_recoil_energies_master(j-1),global_pka_recoil_energies_master(j), &
                 total_pka_sum(j), &
                   total_pka_sum(j)/SUM(total_pka_sum(1:global_num_pka_recoil_points_master)), &
                   SUM(total_pka_sum(1:j))/SUM(total_pka_sum(1:global_num_pka_recoil_points_master))
                   
       END DO 
      END IF
  
   
       ! compute average pka energy
       pka_ave=0.5_DBL*global_pka_recoil_energies_master(1)*total_pka_sum(1)
       DO j=1,global_num_pka_recoil_points_master
        pka_ave=pka_ave+0.5_DBL*(global_pka_recoil_energies_master(j)+&
                                 global_pka_recoil_energies_master(j-1))*total_pka_sum(j)
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
    END IF !non zero check 
   
  END SUBROUTINE output_global_sums
