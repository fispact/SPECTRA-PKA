!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk,                                                     *
!******************************************************************************
!end rubric

program specter
use globals
 IMPLICIT NONE

 
 io_quit=0


 
 !read_input
 CALL read_input()
 
 
 ! read flux file
 
 IF(io_quit==0) CALL read_flux()
 

 
 !omit unnecessary compound data read
 
 ! open results file
 IF(io_quit==0) THEN
  !3/5/16
  IF(LEN(TRIM(results_stub))==0) THEN
   ! default to original filename method (default for results_stub is empty)
   OPEN(results_unit,FILE=TRIM(ADJUSTL(results_filename)),STATUS='REPLACE',IOSTAT=io_open)
  ELSE 
   OPEN(results_unit,FILE=TRIM(ADJUSTL(results_stub))//'.out',STATUS='REPLACE',IOSTAT=io_open)
  END IF
  IF(io_open.NE.0) THEN
   PRINT *,'error opening results file'
   io_quit=1
  END IF
 END IF
 
 file_index=0
 IF(io_quit==0) THEN
  !3/5/16
  IF(LEN(TRIM(results_stub))==0) THEN
   ! default to original filename method (default for results_stub is empty)
   OPEN(index_summary,FILE='index_summary.dat',STATUS='REPLACE',IOSTAT=io_open)
  ELSE 
   OPEN(index_summary,FILE=TRIM(ADJUSTL(results_stub))//'.indexes',STATUS='REPLACE',IOSTAT=io_open)
  END IF 
 
 IF(io_open.NE.0) THEN
  PRINT *,'error opening index summary file'
  io_quit=1
 END IF
 END IF
 
 IF(io_quit==0) THEN
   !output spectrum info
   WRITE(results_unit,*) '### index ',file_index,' #####'
   WRITE(results_unit,*) '#'//TRIM(ADJUSTL(flux_title))
   WRITE(results_unit,*) '#INPUT SPECTRUM - E(MEV) lower+upper bound  VS.  FLUX'
   DO i=1,number_flux_groups
    WRITE(results_unit,'(1x,4ES12.4)') flux_energies(i),flux_energies(i+1),fluxes(i)
   END DO                              
   
   ! normalise flux distribution
   rtemp=0._DBL
   DO i=1,number_flux_groups
    flux_ebin_widths(i)=flux_energies(i+1)-flux_energies(i) !fluxes in n/c^2/s/MeV units
    IF(igroup==2) flux_ebin_widths(i)=1._DBL   ! fluxes in n/cm^2/s units
    rtemp=rtemp+fluxes(i)*flux_ebin_widths(i)
   END DO
   total_fluence=rtemp*time/acnm
   norm_fluence=1._DBL/rtemp
   IF(flux_norm_type==2) THEN
    fluxes_norm=fluxes*1e-24_DBL  ! 1/cm^2 --> 1/barns
   ELSEIF(flux_norm_type==3) THEN
    ! do nothing
    fluxes_norm=fluxes
   ELSE
    ! 1/11/2013 - even here we want to convert between cm2 and barns
    ! 1/11/2013 - and don't want to normalise
    fluxes_norm=fluxes*flux_ebin_widths*1e-24_DBL
   END IF
   WRITE(results_unit,'(2(A,ES11.4))') '# TOTAL FLUX = ',rtemp,'  FLUENCE = ',total_fluence
   
   WRITE(results_unit,*)
   WRITE(results_unit,*)
   WRITE(index_summary,*) file_index,' input_spectrum'
   file_index=file_index+1
   

END IF


 filenum=1
 first_non_empty=.false.
 DO WHILE ((io_quit==0).AND.(filenum.LE.number_pka_files))
! read and process pka data
 at_end=.false.
 num_pka_elements=0
 alpha_sum_flag=.true.
 proton_sum_flag=.true.
 !open pka file
 OPEN(unit=pka_unit,FILE=TRIM(ADJUSTL(pka_filename(filenum))),IOSTAT=io_open,STATUS='OLD')
 
  
 IF(io_open==0) THEN
  first_read=.true.
  PRINT *,'reading ',TRIM(ADJUSTL(pka_filename(filenum)))
  io_read=0
  DO WHILE(.not.at_end.AND.(io_read==0).AND.(io_quit==0))
   ! read next set of data from pka file
   
   CALL read_pka()
   IF(at_end) cycle
   IF(io_read.NE.0) cycle
   IF(io_quit.NE.0) cycle
   PRINT *,'pka data read in for: '//TRIM(ADJUSTL(pka_element))
   !num_pka_elements=num_pka_elements+1
   !IF(num_pka_elements==69) at_end=.true.
   !write(*,*) 'got here',num_pka_elements,ALLOCATED(pka_sums)
   !ALLOCATE(pka_fluxes(num_pka_recoil_points,MAX(number_flux_ebins,num_pka_incident_energies)))
   
   !IF(io_read==0) CALL read_sigma




   !6/3/2014 - now based on cross section read
   IF(do_ngamma_estimate.AND.(mtd==102).AND.(index(pka_element,"cross").NE.0)) THEN
     CALL process_ng_estimate
     !! and that is all we do for this version, so cycle
     !DEALLOCATE(recoil_kermas,pka_recoil_energies,pka_incident_energies)
     !CYCLE
   END IF
   
   
   num_pka_elements=num_pka_elements+1
   
  IF(index(pka_element,"matrix").NE.0) THEN
   !PRINT *,recoil_kermas(1,:)
   !PRINT *,recoil_kermas(2,:)
   !PRINT *,recoil_kermas(3,:)
   
       !IF (i==570) THEN
       !j=615
       !    WRITE(number_string,'(I5)') daughter_num
       !    WRITE(number_string2,'(I5)') parent_num(filenum)
       !    WRITE(600,*) '#  '//pka_element//' [ '//TRIM(ADJUSTL(daughter_ele)) &
       !       //TRIM(ADJUSTL(number_string))//' ] from [ '//TRIM(ADJUSTL(parent_ele(filenum))) &
       !       //TRIM(ADJUSTL(number_string2)) &
       !//' ]'
       !non_zero=.false.
       ! DO i=1,num_pka_points
       !  IF(recoil_kermas(j,i)>0d0) WRITE(600,*) 'in',pka_incident_energies(i),pka_incident_energies(i+1),recoil_kermas(j,i)
       !   !non_zero=.true.
       !  !END IF
       ! END DO
       ! WRITE(600,*)
       ! WRITE(600,*)
     ! END IF
   DO i=1,num_pka_recoil_points
    !pka_fluxes(i,1:num_pka_incident_energies)
    ! collapse input_pka_energy_spectrum onto flux spectrum

    CALL collapse_xs2(recoil_kermas(i,:),num_pka_points,number_flux_groups,flux_energies,pka_incident_energies)
    
    WHERE(recoil_kermas(i,:).LT.0._DBL) recoil_kermas(i,:)=0._DBL
   END DO
   !PRINT *,recoil_kermas(1,:)
   !PRINT *,recoil_kermas(2,:)
   !PRINT *,recoil_kermas(3,:)
   
       !j=615
       !    WRITE(600,*) '#  '//pka_element//' [ '//TRIM(ADJUSTL(daughter_ele)) &
       !       //TRIM(ADJUSTL(number_string))//' ] from [ '//TRIM(ADJUSTL(parent_ele(filenum))) &
       !       //TRIM(ADJUSTL(number_string2)) &
       !//' ]'       
       ! DO i=1,number_flux_groups
       !  IF(recoil_kermas(j,i)>0d0) WRITE(600,*) 'out',flux_energies(i),flux_energies(i+1),recoil_kermas(j,i)
       ! END DO
       ! WRITE(600,*)
       ! WRITE(600,*)   
   
  END IF

   
   ALLOCATE(pka(11,num_pka_recoil_points),epka(num_pka_recoil_points))
   IF((num_pka_elements==1).AND.(do_mtd_sums)) THEN
    num_pka_recoil_points_master=num_pka_recoil_points

    ALLOCATE(pka_sums(10,num_pka_recoil_points_master),pka_recoil_energies_master(num_pka_recoil_points_master))
    IF(do_tdam) THEN
     ALLOCATE(mtd_disp_sums(10))
     mtd_disp_sums=0._DBL
     
      ALLOCATE(tdam_energies_master(num_pka_recoil_points_master), &
               pka_sums_tdam(10,num_pka_recoil_points_master))
      !29/6/2016 can only  have tdams for varying parent/daughter combinations
      ! using a fixed vector with interpolation
      tdam_energies_master=pka_recoil_energies
      pka_sums_tdam=0._DBL
     
    END IF
    pka_sums=0._DBL
    
    pka_recoil_energies_master=pka_recoil_energies
   END IF
   
   ! for new style input files - just one energy array per file - fixed incident energy spectrum (and recoil)
   IF(energies_once_perfile.and.first_read.and.(.not.first_non_empty)) THEN
    num_pka_incident_energies_master=num_pka_incident_energies
    ALLOCATE(pka_incident_energies_master(num_pka_incident_energies_master))
    
    pka_incident_energies_master=pka_incident_energies
    IF(.not.do_mtd_sums) THEN
     num_pka_recoil_points_master=num_pka_recoil_points
     ALLOCATE(pka_recoil_energies_master(num_pka_recoil_points_master))
     pka_recoil_energies_master=pka_recoil_energies

    END IF
    
   END IF
   
   
   first_read=.false.
   ! set up arrays for global isotope and element sums
   IF((num_pka_elements==1).AND.(.not.first_non_empty).AND.(do_global_sums)) THEN
    global_num_pka_recoil_points_master=num_pka_recoil_points
    ALLOCATE(global_pka_sums(max_global_recoils,global_num_pka_recoil_points_master), &
       global_pka_recoil_energies_master(global_num_pka_recoil_points_master), &
       global_daughter_eles(max_global_recoils),global_daughter_nums(max_global_recoils))
    IF(do_tdam) THEN
      ALLOCATE(global_disp_sums(max_global_recoils))
      global_disp_sums=0._DBL
      total_disp=0._DBL
    END IF
    global_pka_sums=0._DBL
    global_pka_recoil_energies_master=pka_recoil_energies
    IF(do_tdam) THEN
     ALLOCATE(global_tdam_energies_master(global_num_pka_recoil_points_master),&
         global_pka_sums_tdam(max_global_recoils,global_num_pka_recoil_points_master) &
         )
      ! can only have tdams for globals with varying parent/daughter combinations
      ! using a fixed vector with interpolation
     global_tdam_energies_master=pka_recoil_energies
     global_pka_sums_tdam=0._DBL
    END IF
    number_global_recoils=0
    
    
    
    !10/10/2013 element sums
    ALLOCATE(global_pka_sums_element(max_global_recoils,global_num_pka_recoil_points_master), &
             global_elements(max_global_recoils))
    IF(do_tdam) THEN
     ALLOCATE(global_disp_sums_element(max_global_recoils))
     global_disp_sums_element=0._DBL
    END IF
    number_global_recoil_elements=0
    global_pka_sums_element=0._DBL


    
    !20/5/2014
    ALLOCATE(total_pka_sum(global_num_pka_recoil_points_master))
    total_pka_sum=0._DBL
    
    IF(do_tdam) THEN
     ALLOCATE(global_pka_sums_element_tdam(max_global_recoils,&
                                           global_num_pka_recoil_points_master), &
              total_pka_sum_tdam(global_num_pka_recoil_points_master))
     total_pka_sum_tdam=0._DBL
     global_pka_sums_element_tdam=0._DBL
    END IF    
    
   END IF
   
  IF(index(pka_element,"matrix").NE.0) THEN 
   DO i=1,num_pka_recoil_points
    pka(1,i)=SUM(fluxes_norm(1:number_flux_ebins)*recoil_kermas(i,1:number_flux_ebins)) 
    epka(i)=0._DBL
   IF(ksail.GE.0) THEN
    DO j=1,number_flux_ebins
     epka(i)=epka(i)+fluxes_norm(j)*recoil_kermas(i,j)*SUM( &
                   fluxes_norm(1:number_flux_ebins)*recoil_kermas(i,1:number_flux_ebins)*&
                  flux_covariances(j,1:number_flux_ebins) &
                                                          )
    END DO   
    !print *,epka(i),pka(1,i)
    IF(pka(1,i).GT.0._DBL) THEN
     epka(i)=100._DBL*SQRT(epka(i) &
                          )/pka(1,i)
    ELSE
     epka(i)=0._DBL
    END IF
   END IF ! ksail
   END DO

   SELECT CASE(pka_filetype) 
    CASE(3)
    CASE DEFAULT
     !print *,mtd
     IF(INDEX(pka_element,"recoil").NE.0) THEN
      CALL define_daughter(.true.)
     ELSE
      CALL define_daughter(.false.) ! light particle recoil
     END IF
    END SELECT
    print *,'daughter: ',daughter_ele,daughter_num,daughter_z
    
    
    IF(do_tdam) THEN
        IF (daughter_num==-1) THEN
         !unknown - assume parent masses
         CALL calc_tdam(num_pka_recoil_points,pka_recoil_energies,tdam_energies, &
         parent_num(filenum),parent_z,parent_num(filenum),parent_z)
        ELSE
         
         CALL calc_tdam(num_pka_recoil_points,pka_recoil_energies,tdam_energies, &
         daughter_num,daughter_z,parent_num(filenum),parent_z)
        END IF
        !tdam_energies=tdam_energies/1e6
         
         !4/4/2016 calculate overall displacements estimate for channel.
         displacements=0._DBL
         i=1
         ! skip NRT part until end
         !displacements=(0.8_DBL*pka(1,i)/(2._DBL*assumed_ed*1e-6_DBL))*&
         !    (tdam_energies(i)/2._DBL)
         !do i=2,num_pka_recoil_points
         !  displacements=displacements+(0.8_DBL*pka(1,i)/(2._DBL*assumed_ed*1e-6_DBL))*&
         !    ((tdam_energies(i-1)+tdam_energies(i))/2._DBL)
         !end do

         displacements=(pka(1,i))*&
             (tdam_energies(i)/2._DBL)

         do i=2,num_pka_recoil_points
           displacements=displacements+(pka(1,i))*&
             ((tdam_energies(i-1)+tdam_energies(i))/2._DBL)
         end do         
         
    END IF

    
    !10/10/2013
    ! needs to come before sum_pkas - so that we can check alpha/proton sum flags.
    IF(do_global_sums) THEN
    print *,'global sums'
     IF( ((.not.alpha_sum_flag).AND.(mtd.GE.800).AND.(mtd.LE.849)).OR. &
         ((.not.proton_sum_flag).AND.(mtd.GE.600).AND.(mtd.LE.649))) THEN
        ! do not add because main 107/103 z,a/z,p reaction has already been added
        ! (assumed to come first)
        ! and reaction falls within limits of excited state z,a/z,p mtd numbers
     ELSE
      CALL add_to_globals(pka(1,:),num_pka_recoil_points,tdam_energies)
     END IF
    END IF
  
   
   pka(3,:)=pka(1,:)
   IF(do_mtd_sums) THEN
    print *,'mtd sums'
    CALL sum_pkas()
   END IF
   !print *,'output',sum(pka(1,1:num_pka_recoil_points)),io_quit

   IF(io_quit.NE.0) cycle

    
    ! and write out
    !20/6/2013 - but only if non-zero
   IF(sum(pka(1,1:num_pka_recoil_points)).NE.0) THEN
    
    WRITE(results_unit,*) '### index ',file_index,' ##### ',TRIM(ADJUSTL(pka_filename(filenum)))
    WRITE(number_string,'(I5)') daughter_num
    WRITE(number_string2,'(I5)') parent_num(filenum)
    WRITE(results_unit,*) '#  '//pka_element//' [ '//TRIM(ADJUSTL(daughter_ele)) &
       //TRIM(ADJUSTL(number_string))//' ] from [ '//TRIM(ADJUSTL(parent_ele(filenum))) &
       //TRIM(ADJUSTL(number_string2)) &
       //' ]'
    WRITE(index_summary,*) file_index,' ',TRIM(ADJUSTL(pka_element))//' [ '//TRIM(ADJUSTL(daughter_ele)) &
       //TRIM(ADJUSTL(number_string))//' ] from [ '//TRIM(ADJUSTL(parent_ele(filenum))) &
       //TRIM(ADJUSTL(number_string2)) &
       //' ]'//TRIM(ADJUSTL(pka_filename(filenum)))  
    WRITE(index_summary,*) file_index+1,' normalised ',TRIM(ADJUSTL(pka_element))//' [ '//TRIM(ADJUSTL(daughter_ele)) &
       //TRIM(ADJUSTL(number_string))//' ] from [ '//TRIM(ADJUSTL(parent_ele(filenum))) &
       //TRIM(ADJUSTL(number_string2)) &
       //' ]'//TRIM(ADJUSTL(pka_filename(filenum)))   
   ! WRITE(results_unit,*) '#DIFF SPEC AVG. PKA RECOIL DISTRIBUTIONS'
   WRITE(results_unit,*) '#PKA RECOIL DISTRIBUTIONS'
   IF(do_tdam) THEN


    IF(ksail.GE.0) THEN
     WRITE(results_unit,'(1x,a32,a8,3x,a7,a11,a28,A22,A12)') &
           '#RECOIL energy (MeV low & high)','SUM','ERROR(%)',&
              'norm_sum','   T_dam (MeV low & high)','disp_energy (eV/s)','dpa/s'
     i=1
     ! add in displacement energy and dpa
     IF(pka(1,i).NE.0) WRITE(results_unit,'(2ES16.4,ES11.4,3x,F7.3,2ES11.4,2ES20.4)') &
               pka_recoil_energies(i)/2._DBL,pka_recoil_energies(i), &
                   pka(1,i),epka(i),pka(1,i)/SUM(pka(1,1:num_pka_recoil_points)), &
                   tdam_energies(i)/2._DBL,tdam_energies(i),(pka(1,i))*&
             (tdam_energies(i)/2._DBL)*1e6_DBL, &
             0.8_DBL*(pka(1,i))*&
             (tdam_energies(i)/2._DBL)/(2._DBL*assumed_ed*1e-6_DBL)
     DO i=2,num_pka_recoil_points
      IF(pka(1,i).NE.0) WRITE(results_unit,'(2ES16.4,ES11.4,3x,F7.3,2ES11.4,2ES20.4)') &
               pka_recoil_energies(i-1),pka_recoil_energies(i), &
                   pka(1,i),epka(i),pka(1,i)/SUM(pka(1,1:num_pka_recoil_points)), &
                   tdam_energies(i-1),tdam_energies(i),(pka(1,i))*&
             ((tdam_energies(i-1)+tdam_energies(i))/2._DBL)*1e6_DBL, &
             0.8_DBL*(pka(1,i))*&
             ((tdam_energies(i-1)+tdam_energies(i))/2._DBL)/(2._DBL*assumed_ed*1e-6_DBL)
      
     END DO
    ELSE
     WRITE(results_unit,'(1x,a32,a8,a11,a28,A22,A12)') '#RECOIL energy (MeV low & high)',&
            'SUM','norm_sum', &
            '   T_dam (MeV low & high)','disp_energy (eV/s)','dpa/s'
     i=1
     IF(pka(1,i).NE.0) WRITE(results_unit,'(2ES16.4,4ES11.4,2ES20.4)') &
               pka_recoil_energies(i)/2._DBL,pka_recoil_energies(i), &
                   pka(1,i),pka(1,i)/SUM(pka(1,1:num_pka_recoil_points)), &
                  tdam_energies(i)/2._DBL,tdam_energies(i),(pka(1,i))*&
             (tdam_energies(i)/2._DBL)*1e6_DBL, &
             0.8_DBL*(pka(1,i))*&
             (tdam_energies(i)/2._DBL)/(2._DBL*assumed_ed*1e-6_DBL) 
     DO i=2,num_pka_recoil_points
      IF(pka(1,i).NE.0) WRITE(results_unit,'(2ES16.4,4ES11.4,2ES20.4)') &
               pka_recoil_energies(i-1),pka_recoil_energies(i), &     
                 pka(1,i),pka(1,i)/SUM(pka(1,1:num_pka_recoil_points)), &
                   tdam_energies(i-1),tdam_energies(i),(pka(1,i))*&
             ((tdam_energies(i-1)+tdam_energies(i))/2._DBL)*1e6_DBL, &
             0.8_DBL*(pka(1,i))*&
             ((tdam_energies(i-1)+tdam_energies(i))/2._DBL)/(2._DBL*assumed_ed*1e-6_DBL)
     END DO    
    
    END IF

    
   ELSE
    IF(ksail.GE.0) THEN
     WRITE(results_unit,'(1x,a32,a8,3x,a7,a11)') &
           '#RECOIL energy (MeV low & high)','SUM','ERROR(%)','norm_sum'
     i=1
     IF(pka(1,i).NE.0) WRITE(results_unit,'(2ES16.4,ES11.4,3x,F7.3,ES11.4)') &
               pka_recoil_energies(i)/2._DBL,pka_recoil_energies(i), &
                   pka(1,i),epka(i),pka(1,i)/SUM(pka(1,1:num_pka_recoil_points))
     DO i=2,num_pka_recoil_points
      IF(pka(1,i).NE.0) WRITE(results_unit,'(2ES16.4,ES11.4,3x,F7.3,ES11.4)') &
               pka_recoil_energies(i-1),pka_recoil_energies(i), &
                   pka(1,i),epka(i),pka(1,i)/SUM(pka(1,1:num_pka_recoil_points))
      
     END DO
    ELSE
     WRITE(results_unit,'(1x,a32,a8,a11)') '#RECOIL energy (MeV low & high)','SUM','norm_sum'
     i=1
     IF(pka(1,i).NE.0) WRITE(results_unit,'(2ES16.4,2ES11.4)') &
               pka_recoil_energies(i)/2._DBL,pka_recoil_energies(i), &
                   pka(1,i),pka(1,i)/SUM(pka(1,1:num_pka_recoil_points))
     DO i=2,num_pka_recoil_points
      IF(pka(1,i).NE.0) WRITE(results_unit,'(2ES16.4,2ES11.4)') &
               pka_recoil_energies(i-1),pka_recoil_energies(i), &     
                 pka(1,i),pka(1,i)/SUM(pka(1,1:num_pka_recoil_points))
     END DO    
    
    END IF
   END IF
   IF (do_tdam) then
    WRITE(results_unit,'(a1,20x,A,ES11.4)') '#',&
               'displacement energy eV/s = ',displacements*1e6_DBL
    WRITE(results_unit,'(a1,20x,A,ES11.4,A7,F4.1,a)') '#','equivalent NRT dpa/s = ',&
              0.8_DBL*displacements/(2._DBL*assumed_ed*1e-6_DBL),&
            ' (E_d=',assumed_ed,' eV)'
   END IF
    ! 11/3/2018 - unused normalised output omitted - remove blank lines here - will 
    ! be added when kept average energy is output below
    !WRITE(results_unit,*)
    !WRITE(results_unit,*)
    !file_index=file_index+1
    
    ! normalise pka distributions to unity
    DO i=1,1! 1 only   ,3,2 ! 1 and 3 only
     IF(i==1) sum_temp=pka_recoil_energies(1)*pka(i,1)
     sum_temp2=0._DBL
     DO j=1,num_pka_recoil_points-1
      !psv=pka(i,j)
      IF((pka(i,j).NE.0._DBL).AND.(pka(i,j+1).NE.0._DBL)) THEN
       pka(i,j)=exp(0.5_DBL*(LOG(pka(i,j))+LOG(pka(i,j+1))))*(pka_recoil_energies(j+1)-pka_recoil_energies(j))
      ELSE
       pka(i,j)=0.5_DBL*(pka(i,j)+pka(i,j+1))*(pka_recoil_energies(j+1)-pka_recoil_energies(i))
      END IF
      IF(i==1) THEN
       sum_temp=sum_temp+pka(i,j)
      ELSE
       !mqq=0
       !IF(pka_recoil_energies(j+1).LE.
       ! unknown damage limits -- EDX
      END IF
     END DO !j
     pka(i,num_pka_recoil_points)=0._DBL
     IF(sum_temp.NE.0._DBL) THEN
      IF(i==1) THEN      
       pka(i,:)=pka(i,:)/sum_temp
      ELSE
       pka(i,:)=pka(i,:)/sum_temp2
      END IF
     END IF
    END DO ! i
    pka(2,1)=pka(1,1)
    pka(4,1)=0._DBL
    DO i=2,num_pka_recoil_points
     pka(4,i)=pka(4,i-1)+pka(3,i)
     pka(2,i)=pka(2,i-1)+pka(1,i)
    END DO
    ! 11/3/2018 - unused normalised output - remove
    !WRITE(results_unit,*) '### index ',file_index,' ##### ',TRIM(ADJUSTL(pka_filename(filenum)))
    !WRITE(results_unit,'(a,ES11.4)') '# NORMALIZED SPECTRAL AND GROUP AVERAGED RECOIL DIST'// &
    ! 'RIBUTIONS  SUM = ',sum_temp
    !WRITE(results_unit,*) '#  ERECOIL  SUM-0-INTEGRAL  '
    !DO i=1,num_pka_recoil_points
    ! IF(pka(1,i)==0._DBL) cycle
    ! WRITE(results_unit,'(1x,3ES10.3,I10)') pka_recoil_energies(i),(pka(j,i),j=1,2),i
    !END DO
    ! compute average pka energy
    pka_ave=0.5_DBL*pka_recoil_energies(1)*pka(1,1)
    DO i=2,num_pka_recoil_points
     pka_ave=pka_ave+0.5_DBL*(pka_recoil_energies(i)+pka_recoil_energies(i-1))*pka(1,i)
    END DO
    pka_ave=pka_ave*1000000._DBL
    !11/3/2018 - but keep this part - still usefil
    WRITE(results_unit,'(a1,20x,A,ES11.4,A5)') '#','AVERAGE PKA ENERGY = ',pka_ave,' (eV)'
    
    WRITE(results_unit,*)
    WRITE(results_unit,*)
    file_index=file_index+1
    
   ELSE
     print *,'zero recoils for this channel'
   END IF !non zero check
   
  END IF ! cross sectin skip check
   
   DEALLOCATE(recoil_kermas,pka_recoil_energies,pka_incident_energies,pka,epka)
   IF(do_tdam) DEALLOCATE(tdam_energies)
   
   !16/10/2017 - if we get here then the file contained something to be used
   ! set first flag to true to prevent future attempts to allocate global arrays
   ! this prevents exception caused by empty file 1
    first_non_empty=.true.
   
  END DO ! over pka matrices in file
  !STOP
  IF((do_mtd_sums).AND.ALLOCATED(pka_sums)) THEN
    CALL output_sum_pkas()
    DEALLOCATE(pka_recoil_energies_master,pka_sums)
    IF(do_tdam) DEALLOCATE(mtd_disp_sums,pka_sums_tdam,tdam_energies_master)
  END IF
  
  !print *,do_ngamma_estimate,ALLOCATED(ng_pka)
  IF(do_ngamma_estimate.AND.ALLOCATED(ng_pka)) THEN
   CALL output_ng_estimate
   DEALLOCATE(ng_pka,ng_pka_recoil_energies)
   IF(do_tdam) DEALLOCATE(ng_tdam_energies)
  END IF
  
  CLOSE(pka_unit)

  
  
  
 ELSE !pka file io_open
   PRINT *,'error opening pka or cross section files'
   io_quit=1
 END IF

 filenum=filenum+1
   
 END DO ! io_quit  and loop over filenames for pka section
 
 IF(do_global_sums) THEN
   CALL output_global_sums()
   IF(ALLOCATED(global_pka_sums)) DEALLOCATE(global_pka_sums)
   IF(ALLOCATED(global_pka_sums_element)) DEALLOCATE(global_pka_sums_element)
   IF(do_tdam.AND.ALLOCATED(global_tdam_energies_master)) &
       DEALLOCATE(global_tdam_energies_master,global_pka_sums_tdam, &
      global_pka_sums_element_tdam)
 END IF
 IF(io_quit.NE.0) THEN
  PRINT *,'quiting'
 END IF
 CLOSE(results_unit)
 
 CONTAINS
 

 end program specter
 
 
 
 
 
 





