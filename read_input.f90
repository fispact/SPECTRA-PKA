!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk,                                                     *
!******************************************************************************
!end rubric

      SUBROUTINE read_input()
      use globals
      IMPLICIT NONE
      INTEGER :: equal_position,iipka
      character*1000 :: input_line,input_label,input_value
      logical :: processed
      CHARACTER (LEN=1000), ALLOCATABLE :: column_data(:,:)
      
      
    CALL GETARG(1,input_line)
    IF(LEN(TRIM(ADJUSTL(input_line)))==0) THEN
     open(unit=unitin,file='specter.input',IOSTAT=io_open) 
    ELSE
     open(unit=unitin,file=TRIM(ADJUSTL(input_line)),IOSTAT=io_open,STATUS='OLD')
    END IF
    IF(io_open.NE.0) THEN
     write(*,*) 'error opening parameter input file - quiting'
     stop
     io_quit=1
    END IF

      ! first run to identify number of pka files and number of columns
      ! for column data entry
          io_read=0
          number_pka_files=1
          num_columns=1
          DO 
          
            read (unitin,end=100,fmt='(a1000)') input_line
            input_line=ADJUSTL(input_line)
            if (input_line(1:1)=='#') cycle
            if (len_trim(input_line)==0) cycle
            equal_position=index(input_line,'=')
            input_label=input_line(1:equal_position-1)
            input_value=input_line(equal_position+1:len(input_line))         
            if (trim(adjustl(input_label))=='number_pka_files') then
             read (input_value,*) number_pka_files
             write (*,*) 'changing number_pka_files to      :',number_pka_files
            end if
            if (trim(adjustl(input_label))=='num_columns') then
             read (input_value,*) num_columns
             write (*,*) 'changing num_columns to      :',num_columns
            end if            
          END DO
          

100       REWIND(unitin)

          ALLOCATE(pka_filename(number_pka_files),pka_ratios(number_pka_files), &
                   parent_ele(number_pka_files),parent_num(number_pka_files), &
                   ngamma_daughter_mass(number_pka_files), &
                   ngamma_parent_mass(number_pka_files))
          
      !
      ! Set default input parameters
      !
          call set_defaults
 
          io_read=0
          DO 
          
            read (unitin,end=200,fmt='(a1000)') input_line
            input_line=ADJUSTL(input_line)
            if (input_line(1:1)=='#') cycle
            if (len_trim(input_line)==0) cycle
            processed=.false.
            equal_position=index(input_line,'=')
            input_label=input_line(1:equal_position-1)
            input_value=input_line(equal_position+1:len(input_line))   
            
            
            if (trim(adjustl(input_label))=='number_pka_files') then
             write (*,*) 'skipping over number_pka_files in this second pass'
             processed=.true.
            end if  
            if (trim(adjustl(input_label))=='num_columns') then
             write (*,*) 'skipping over num_columns in this second pass'
             processed=.true.
            end if              
            
            
           if (trim(adjustl(input_label))=='columns') then
            ALLOCATE (column_data(number_pka_files,num_columns))
            iipka=1
            DO WHILE ((iipka.LE.number_pka_files).AND.(io_read==0))
             READ(unitin,*,IOSTAT=io_read) column_data(iipka,1:num_columns)
             iipka=iipka+1
            END DO
            IF(io_read==0) THEN
             CALL process_columns(num_columns,column_data,input_value,processed)
            ELSE
             WRITE(*,*) 'error reading raw column data'
            END IF
            IF((io_read==0).AND.(processed)) write (*,*) 'processed_column_data'
            !processed=.true.
           end if             
            
            
            
            if (trim(adjustl(input_label))=='flux_filename') then
             read (input_value,*) flux_filename
             write (*,*) 'changing flux_filename to      :',flux_filename
             processed=.true.
            end if                        
            if (trim(adjustl(input_label))=='results_filename') then
             read (input_value,*) results_filename
             write (*,*) 'changing results_filename to      :',results_filename
             processed=.true.
            end if  
            if (trim(adjustl(input_label))=='kerma_filename') then
             read (input_value,*) kerma_filename
             write (*,*) 'changing kerma_filename to      :',kerma_filename
             processed=.true.
            end if   
            if (trim(adjustl(input_label))=='pka_filename') then
             BACKSPACE(unitin) !go back one line and read pka_filename= into a dummy string
             ! 110/10/2013 allows multiple long filenames to be split on multiple lines
             ! must be a space between = and first filename
             read (unitin,*) input_value,(pka_filename(i),i=1,number_pka_files)
             write (*,*) 'changing pka_filename to      :',pka_filename
             processed=.true.
            end if     
            if (trim(adjustl(input_label))=='pka_ratios') then
             read (input_value,*) (pka_ratios(i),i=1,number_pka_files)
             write (*,*) 'changing pka_ratios to      :',pka_ratios
             processed=.true.
             IF(SUM(pka_ratios(1:number_pka_files))>1._DBL) THEN
              write(*,*) 'ratios for summing results for different pka matrices greater than 1'
              write(*,*) 'normalising ratios to 1'
              pka_ratios=pka_ratios/SUM(pka_ratios(1:number_pka_files))
             END IF
            end if              
            if (trim(adjustl(input_label))=='gas_filename') then
             read (input_value,*) gas_filename
             write (*,*) 'changing gas_filename to      :',gas_filename
             processed=.true.
            end if
            
            if (trim(adjustl(input_label))=='sigma_filename') then
 	     read (input_value,*) sigma_filename
 	     write (*,*) 'changing sigma_filename to      :',sigma_filename
 	     processed=.true.
            end if            
            !if (trim(adjustl(input_label))=='nkerma_groups_input') then
            ! read (input_value,*) nkerma_groups_input
            ! write (*,*) 'changing nkerma_groups_input to      :',nkerma_groups_input
            ! processed=.true.
            !end if           
            !if (trim(adjustl(input_label))=='max_kerma_elements') then
            ! read (input_value,*) max_kerma_elements
            ! write (*,*) 'changing max_kerma_elements to      :',max_kerma_elements
            ! processed=.true.
            !end if
            
             
                       
            if (trim(adjustl(input_label))=='pka_filetype') then
 	     read (input_value,*) pka_filetype
 	     write (*,*) 'changing pka_filetype to      :',pka_filetype
 	     processed=.true.
            end if  

            if (trim(adjustl(input_label))=='flux_energy_rescale_value') then
 	     read (input_value,*) flux_energy_rescale_value
 	     write (*,*) 'changing flux_energy_rescale_value to      :',flux_energy_rescale_value
 	     processed=.true.
            end if 

            if (trim(adjustl(input_label))=='flux_norm_type') then
 	     read (input_value,*) flux_norm_type
 	     write (*,*) 'changing flux_norm_type to      :',flux_norm_type
 	     processed=.true.
            end if 

            if (trim(adjustl(input_label))=='do_mtd_sums') then
 	     read (input_value,*) do_mtd_sums
 	     write (*,*) 'changing do_mtd_sums to      :',do_mtd_sums
 	     processed=.true.
            end if 

            if (trim(adjustl(input_label))=='do_ngamma_estimate') then
 	     read (input_value,*) do_ngamma_estimate
 	     write (*,*) 'changing do_ngamma_estimate to      :',do_ngamma_estimate
 	     processed=.true.
            end if 
            
          

            if (trim(adjustl(input_label))=='parent_ele') then
             read (input_value,*) (parent_ele(i),i=1,number_pka_files)
             write (*,*) 'changing parent_ele to      :',parent_ele
             processed=.true.
            end if
            
            if (trim(adjustl(input_label))=='parent_num') then
             read (input_value,*) (parent_num(i),i=1,number_pka_files)
             write (*,*) 'changing parent_num to      :',parent_num
             processed=.true.
            end if 

            if (trim(adjustl(input_label))=='ngamma_daughter_mass') then
             read (input_value,*) (ngamma_daughter_mass(i),i=1,number_pka_files)
             write (*,*) 'changing ngamma_daughter_mass to      :',ngamma_daughter_mass
             processed=.true.
            end if 
            
            if (trim(adjustl(input_label))=='ngamma_parent_mass') then
             read (input_value,*) (ngamma_parent_mass(i),i=1,number_pka_files)
             write (*,*) 'changing ngamma_parent_mass to      :',ngamma_parent_mass
             processed=.true.
            end if             
            if (trim(adjustl(input_label))=='do_exclude_light_from_total') then
             read (input_value,*) do_exclude_light_from_total
             write (*,*) 'changing do_exclude_light_from_total to      :',do_exclude_light_from_total
             processed=.true.
            end if   
            
            if (trim(adjustl(input_label))=='do_global_sums') then
 	     read (input_value,*) do_global_sums
 	     write (*,*) 'changing do_global_sums to      :',do_global_sums
 	     processed=.true.
            end if  
            if (trim(adjustl(input_label))=='max_global_recoils') then
 	     read (input_value,*) max_global_recoils
 	     write (*,*) 'changing max_global_recoils to      :',max_global_recoils
 	     processed=.true.
            end if              
            
            if (trim(adjustl(input_label))=='flux_rescale_value') then
 	     read (input_value,*) flux_rescale_value
 	     write (*,*) 'changing flux_rescale_value to      :',flux_rescale_value
 	     processed=.true.
            end if  
            
            if (trim(adjustl(input_label))=='energies_once_perfile') then
 	     read (input_value,*) energies_once_perfile
 	     write (*,*) 'changing energies_once_perfile to      :',energies_once_perfile
 	     processed=.true.
            end if  

            if (trim(adjustl(input_label))=='do_exclude_unknown_from_total') then
 	     read (input_value,*) do_exclude_unknown_from_total
 	     write (*,*) 'changing do_exclude_unknown_from_total to      :',do_exclude_unknown_from_total
 	     processed=.true.
            end if            

            if (trim(adjustl(input_label))=='do_tdam') then
 	     read (input_value,*) do_tdam
 	     write (*,*) 'changing do_tdam to      :',do_tdam
 	     processed=.true.
            end if  
            if (trim(adjustl(input_label))=='assumed_ed') then
 	     read (input_value,*) assumed_ed
 	     write (*,*) 'changing assumed_ed to      :',assumed_ed
 	     processed=.true.
            end if 
            if (trim(adjustl(input_label))=='tdam_method') then
 	     read (input_value,*) tdam_method
 	     write (*,*) 'changing tdam_method to      :',tdam_method
 	     processed=.true.
            end if   
            
            if (trim(adjustl(input_label))=='incident_particle') then
 	     read (input_value,*) incident_particle
 	     write (*,*) 'changing incident_particle to      :',incident_particle
 	     processed=.true.
            end if            

            if (trim(adjustl(input_label))=='results_stub') then
 	     read (input_value,*) results_stub
 	     write (*,*) 'changing results_stub to      :',results_stub
 	     processed=.true.
            end if            

            
            if (.not.processed) then
                 write (*,*) 'variable ',trim(adjustl(input_label)),' not recognised'
                 STOP
            end if
          END DO 
 200 close (unitin) 
 
 ! check for conflicts
 ! 20/6/2013 pka_filetype and do_mtd_sums
 IF(do_mtd_sums.and.(pka_filetype.ne.2)) THEN
  PRINT *,' pka summing by mt number requested bu not reading from suitable pka matrix file'
  io_quit=1
 END If
 
 
 
      
      END SUBROUTINE read_input
      
      SUBROUTINE set_defaults
      use globals
      IMPLICIT NONE
      !REAL (KIND=DBL) :: 
      
 
 flux_filename='HF32D.DAT'  !expected default is MeV
 flux_energy_rescale_value=1._DBL ! scale factor for flux energy bins - to convert to MeV if not already.

 flux_rescale_value=1._DBL   !24/7/2015 - applied in read flux - would be undone if normalisation requested
                            ! afterwards

 results_filename='results.out'
 results_stub=''
 pka_ratios=1._DBL
! nkerma_groups_input=171
! max_kerma_elements=100
! kerma_filename='MACKLIB.DAT'
 pka_filename='PKAMIN.DAT'
 pka_filetype=2 ! .ne.2 = original format - 2 new format equivalent to groupr output from NJOY
                !7/9/2017 - 3= new geant-generated format from M. Fleming
! gas_filename='GAS157.DAT'
! sigma_filename='SIGD.DAT'
 flux_norm_type=2 ! .ne.2.and.ne.3 = normalise flux before merging with pka - 3= no normalisation
                  !- 2= convert to barns
                  ! 1 normalise spectrum before collapse based on igroup
 
 do_mtd_sums=.false. ! if true (and if pka_filetype==2) then sum certain mtd groups  
 do_ngamma_estimate=.false. ! if true then estimate n,gamma recoil matrix using (n,g) gamma recoils
 !Feb 2014 - alternate method for ng estimate requires mass of nuclide
 ngamma_daughter_mass=56.935392841_DBL ! 26/2/2014 - array of ngamma daughter masses (carbon 12 scale)
                       ! converted to gram mass using avogadro
 ngamma_parent_mass=55.934936326_DBL !6/3/14 - equivalent array of ngamma daughter masses.
 !october 2013 - daughter definitions
 
 parent_ele='Fe'  !used to define daughters
 parent_num=56
 do_global_sums=.false. ! if true then sum species over multiple files
 do_exclude_light_from_total=.true. ! 20/5/2014 if false then include He+H in total recoil spectra.
 do_exclude_unknown_from_total=.true. !4/3/2016
 max_global_recoils=200
 
 energies_once_perfile=.false.  ! 5/10/2015 if false then old style energy structure for every channel+4/3/2016
 do_tdam=.false.  
 assumed_ed=40._DBL  ! 4/4/2016 assumed E_d for NRT displacements per second per channel
 tdam_method=2 ! 1 from NRT paper, 2 from njoy heatr subroutine.
 
 !21/4/2016 need to identify incident particle to get correct daughter
 incident_particle='n' !can be n or p for now
 
 
     END SUBROUTINE set_defaults
     
     
     
SUBROUTINE process_columns(num_col,string_data,headerstring,complete)

 USE globals
 IMPLICIT NONE
 INTEGER, intent(in) :: num_col
 CHARACTER (LEN=*), INTENT(in) :: string_data(number_pka_files,num_col)
 CHARACTER (LEN=*), INTENT(in) :: headerstring
 CHARACTER (LEN=100) :: headers(num_col)
 LOGICAL, INTENT(out) :: complete
 LOGICAL        :: processed
 INTEGER :: ii,jj
 

!allowed data - pka_filename, pka_ratios, parent_ele, parent_num, 
! ngamma_parent_mass, ngamma_daughter_mass

  READ(headerstring,*,IOSTAT=io_read) headers(:)
  IF(io_read.NE.0) WRITE(*,*) 'error reading column headers'

  ii=1
  DO WHILE ((ii.LE.num_col).AND.(io_read==0))
     PRINT *,ii,headers(ii)
     if (trim(adjustl(headers(ii)))=='pka_filename') then
       jj=1
       DO WHILE ((jj.LE.number_pka_files).AND.(io_read==0))
            !read (string_data(jj,ii),'(a)',IOSTAT=io_read) pka_filename(jj)
            pka_filename(jj)=string_data(jj,ii)
            write (*,'(a,I4,A,A)') 'changing pka_filename for nuclide',jj,' to:',pka_filename(jj) 
            jj=jj+1
       END DO       
       processed=.true.
     end if   

     if (trim(adjustl(headers(ii)))=='pka_ratios') then
       jj=1
       DO WHILE ((jj.LE.number_pka_files).AND.(io_read==0))
            read (string_data(jj,ii),*,IOSTAT=io_read) pka_ratios(jj)
            write (*,'(a,I4,A,E15.10)') 'changing pka_ratio for nuclide',jj,' to:',pka_ratios(jj) 
            jj=jj+1
       END DO   
             IF(SUM(pka_ratios(1:number_pka_files))>1._DBL) THEN
              write(*,*) 'ratios for summing results for different pka matrices greater than 1'
              write(*,*) 'normalising ratios to 1'
              pka_ratios=pka_ratios/SUM(pka_ratios(1:number_pka_files))
             END IF       
       processed=.true.
     end if 
     
     if (trim(adjustl(headers(ii)))=='parent_ele') then
       jj=1
       DO WHILE ((jj.LE.number_pka_files).AND.(io_read==0))
            read (string_data(jj,ii),*,IOSTAT=io_read) parent_ele(jj)
            write (*,'(a,I4,A,A)') 'changing parent_ele for nuclide',jj,' to:',parent_ele(jj) 
            jj=jj+1
       END DO       
       processed=.true.
     end if 
     
     if (trim(adjustl(headers(ii)))=='parent_num') then
       jj=1
       DO WHILE ((jj.LE.number_pka_files).AND.(io_read==0))
            read (string_data(jj,ii),*,IOSTAT=io_read) parent_num(jj)
            write (*,'(a,I4,A,I5)') 'changing parent_num for nuclide',jj,' to:',parent_num(jj) 
            jj=jj+1
       END DO       
       processed=.true.
     end if 
     
     if (trim(adjustl(headers(ii)))=='ngamma_parent_mass') then
       jj=1
       DO WHILE ((jj.LE.number_pka_files).AND.(io_read==0))
            read (string_data(jj,ii),*,IOSTAT=io_read) ngamma_parent_mass(jj)
            write (*,'(a,I4,A,E15.10)') 'changing ngamma_parent_mass for nuclide', &
                                         jj,' to:',ngamma_parent_mass(jj) 
            jj=jj+1
       END DO       
       processed=.true.
     end if 
     
     if (trim(adjustl(headers(ii)))=='ngamma_daughter_mass') then
       jj=1
       DO WHILE ((jj.LE.number_pka_files).AND.(io_read==0))
            read (string_data(jj,ii),*,IOSTAT=io_read) ngamma_daughter_mass(jj)
            write (*,'(a,I4,A,E15.10)') 'changing ngamma_daughter_mass for nuclide', &
                                    jj,' to:',ngamma_daughter_mass(jj) 
            jj=jj+1
       END DO       
       processed=.true.
     end if      
     

     
     if (.not.processed) then
         write (*,*) 'column variable ',trim(adjustl(headers(ii))),' not recognised'
         io_quit=1
         io_read=1
     end if 
     IF(io_read.NE.0) THEN
         write(*,*) 'error processing column data'
         io_quit=1
     END IF     
     ii=ii+1
  END DO
  IF(io_read==0) THEN
    complete=.true.
  ELSE
    complete=.false.
  END IF
  
END SUBROUTINE process_columns       
     