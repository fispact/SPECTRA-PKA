!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk,                                                     *
!******************************************************************************
!end rubric

  SUBROUTINE read_pka()
  use globals
  IMPLICIT NONE
  REAL(KIND=DBL) :: rdummy,rdummy2
  INTEGER :: idummy,jdummy,last_i,last_j,flag1,flag2
  CHARACTER (LEN=100) :: input_line,temp_string
  LOGICAL :: redo !5/2/14 - handle case where there are no non-zero matrix elements 
                  !   file type 2 only (from njoy)
     
     redo=.true.
     DO WHILE(redo)
      SELECT CASE(pka_filetype)
      case(1)
       READ(pka_unit,'(a23,a70)',IOSTAT=io_read) pka_element,input_line
       IF(io_read==0) READ(input_line,*,IOSTAT=io_read) rdummy,num_pka_recoil_points,num_pka_points
      CASE(2)
       READ(pka_unit,'(a30,I5,a70)',IOSTAT=io_read) pka_element,mtd,input_line
       
       IF(io_read==0) write (log_unit,'(a,I10,x,a)') TRIM(ADJUSTL(pka_element)),mtd,TRIM(ADJUSTL(input_line))
       
       !6/3/2014 - allow for (n,g) cross section to be read-in
       IF(io_read==0) THEN
        IF(INDEX(pka_element,'matrix')==0) THEN
         ! not pka matrix
          READ(input_line,*,IOSTAT=io_read) rdummy,num_pka_incident_energies
          num_pka_points=num_pka_incident_energies-1
          num_pka_recoil_points=num_pka_incident_energies
        ELSE
         READ(input_line,*,IOSTAT=io_read) rdummy,num_pka_recoil_points,num_pka_points,&
                      flag1,flag2 ! 25/2/2014 - flags to handle non-legendre format from njoy
        END IF
       END IF
       !stop
      CASE(3) !7/9/2017 - GEANT output format
       READ(pka_unit,'(a70)',IOSTAT=io_read) input_line
       
       !IF(io_read==0) print *,TRIM(ADJUSTL(input_line))
       IF(io_read==0) READ(input_line(1:index(input_line,"A:")-LEN("A:")),'(a)',IOSTAT=io_read) pka_element
       IF(io_read==0) THEN
         READ(input_line(index(input_line,"A:")+LEN("A:"):),*,IOSTAT=io_read) daughter_num
         READ(input_line(index(input_line,"Z:")+LEN("Z:"):),*,IOSTAT=io_read) daughter_z
         READ(input_line(index(input_line,"Ein:")+LEN("Ein:"):),*,IOSTAT=io_read) num_pka_incident_energies
         READ(input_line(index(input_line,"Eout:")+LEN("Eout:"):),*,IOSTAT=io_read) num_pka_recoil_points
       END IF
       WRITE(log_unit,*) io_read,daughter_num,daughter_z,num_pka_incident_energies,num_pka_recoil_points
       num_pka_points=num_pka_incident_energies-1
       daughter_ele=master_elements(daughter_z)
       WRITE(log_unit,*) TRIM(ADJUSTL(pka_element))
       !stop
       pka_element=TRIM(ADJUSTL(pka_element))!//" matrix"
       IF((daughter_num==4).AND.(daughter_z==2)) THEN
        pka_element=TRIM(ADJUSTL(pka_element))//" alpha"
       ELSEIF((daughter_num==1).AND.(daughter_z==1)) THEN
        pka_element=TRIM(ADJUSTL(pka_element))//" proton"
       ELSE
        pka_element=TRIM(ADJUSTL(pka_element))//" recoil"
       END IF
       pka_element=TRIM(ADJUSTL(pka_element))//" matrix"
       !STOP

       
      CASE DEFAULT
       READ(pka_unit,'(a23,a70)',IOSTAT=io_read) pka_element,input_line
       IF(io_read==0) READ(input_line,*,IOSTAT=io_read) rdummy,num_pka_recoil_points,num_pka_points
      END SELECT
       
       !WRITE(log_unit,*)  pka_element,rdummy,num_pka_recoil_points,num_pka_points,io_read,input_line
       IF(io_read==0) THEN
        num_pka_incident_energies=num_pka_points+1
        ALLOCATE(pka_recoil_energies(num_pka_recoil_points),pka_incident_energies(num_pka_incident_energies),&
                        recoil_kermas(num_pka_recoil_points,MAX(number_flux_ebins,num_pka_incident_energies)))
        IF(do_tdam) ALLOCATE(tdam_energies(num_pka_recoil_points))
               !5/10/2015 only read energies once if indicated
        IF(INDEX(pka_element,'matrix')==0) THEN
        
         IF((energies_once_perfile.and.first_read).OR.(.not.energies_once_perfile)) THEN
          READ(pka_unit,*,IOSTAT=io_read) (pka_incident_energies(j),j=1,num_pka_incident_energies)
          pka_recoil_energies=pka_incident_energies

         END IF        
        

        ELSE
         SELECT CASE(pka_filetype)
         CASE (3)
          !7/9/2017 - new format has incident spectrum first.
          IF((energies_once_perfile.and.first_read).OR.(.not.energies_once_perfile)) THEN
           READ(pka_unit,*,IOSTAT=io_read) (pka_incident_energies(i),i=1,num_pka_incident_energies),&
                                (pka_recoil_energies(j),j=1,num_pka_recoil_points)

          END IF          
         CASE DEFAULT
          IF((energies_once_perfile.and.first_read).OR.(.not.energies_once_perfile)) THEN
           READ(pka_unit,*,IOSTAT=io_read) (pka_recoil_energies(i),i=1,num_pka_recoil_points),&
                                (pka_incident_energies(j),j=1,num_pka_incident_energies)

          END IF
         END SELECT
        END IF
        
        IF((energies_once_perfile)) THEN
         IF((.not.first_read)) THEN
          IF((num_pka_incident_energies.NE.num_pka_incident_energies_filemaster).OR.&
            (num_pka_recoil_points.NE.num_pka_recoil_points_filemaster)) THEN
            PRINT *,'Energy matrix only read once, but subseqent recoil matrices have different grid'
    WRITE(log_unit,*) 'Energy matrix only read once, but subseqent recoil matrices have different grid' 
	    STOP
          ELSE
           pka_recoil_energies=pka_recoil_energies_filemaster
           pka_incident_energies=pka_incident_energies_filemaster          
          END IF
         ELSE
          num_pka_recoil_points_filemaster=num_pka_recoil_points
          num_pka_incident_energies_filemaster=num_pka_incident_energies
          ALLOCATE(pka_recoil_energies_filemaster(num_pka_recoil_points_filemaster),&
            pka_incident_energies_filemaster(num_pka_incident_energies_filemaster))
           pka_recoil_energies_filemaster=pka_recoil_energies
           pka_incident_energies_filemaster=pka_incident_energies  
          first_read=.false.  
         END IF
        END IF
        SELECT CASE(pka_filetype)
        CASE(1) ! original format
         IF(io_read==0) READ(pka_unit,'(8(E10.3))',IOSTAT=io_read) &
                  ((recoil_kermas(j,i),j=1,num_pka_recoil_points),i=1,num_pka_points)
         redo=.false.
        CASE(2) ! new format from NJOY
         ! in this case NJOY gives the kerma for recoils in group j due to incident neutron in group i
         ! 24/2/2014 there is no overflow to recoil energies below the minimum group
         ! thus there are only num_pka_recoil_points-1 recoil groups and recoil_kermas(1,:) is always zero
         ! for consistency we should assign the jdummy,idummy entry to jdummy+1,idummy
         ! similar for incident energies, but this is the same as previously and so already covered
         last_i=0
         last_j=0
         recoil_kermas=0._DBL
        IF(INDEX(pka_element,'matrix')==0) THEN
         !6/3/2014 not pka matrix - will read and then skip if not n,g
         DO WHILE (io_read==0)
          ! zero and first order legendre polynomials only
          READ(pka_unit,*,IOSTAT=io_read) idummy,rdummy
          IF(io_read==0) THEN
           recoil_kermas(1,idummy)=rdummy
           last_i=idummy
           last_j=1
          ELSE
           ! error or reached end of pka matrix and have read next line
           backspace(pka_unit)
          END IF
         END DO
        ELSEIF((flag1==1).AND.(flag2==1)) THEN !24.2.2014 - read multiple entries per line
         DO WHILE (io_read==0) 
           read(pka_unit,'(A100)',IOSTAT=io_read) input_line
           read(input_line,*,IOSTAT=io_read) idummy,jdummy
           IF(io_read==0) THEN
            ! remove two integers at start of string
            READ(input_line,*) temp_string
            input_line=TRIM(ADJUSTL(input_line(INDEX(input_line,TRIM(ADJUSTL(temp_string)))+&
                   LEN(TRIM(ADJUSTL(temp_string))):LEN(input_line))))
            READ(input_line,*) temp_string
            input_line=TRIM(ADJUSTL(input_line(INDEX(input_line,TRIM(ADJUSTL(temp_string)))+&
                   LEN(TRIM(ADJUSTL(temp_string))):LEN(input_line))))
                   
            DO WHILE((LEN(TRIM(ADJUSTL(input_line))).NE.0).AND.(io_read==0))
             read(input_line,*,IOSTAT=io_read) temp_string
             IF(io_read==0) THEN
              read(temp_string,*) rdummy
              recoil_kermas(jdummy+1,idummy)=rdummy
              !IF(mtd==102) WRITE(107,*) jdummy+1,idummy,rdummy
              jdummy=jdummy+1 ! next one read from line will be for this group
              input_line=TRIM(ADJUSTL(input_line(INDEX(input_line,TRIM(ADJUSTL(temp_string)))+&
                   LEN(TRIM(ADJUSTL(temp_string))):LEN(input_line)))) 
              
             END IF
            END DO
            last_i=idummy
            last_j=jdummy-1
           ELSE
            !error or reached end of pka matrix
            backspace(pka_unit)
           END IF           
         END DO  
        ELSE ! normal legendre polynomial read
         DO WHILE (io_read==0)
          ! zero and first order legendre polynomials only
          READ(pka_unit,*,IOSTAT=io_read) idummy,jdummy,rdummy,rdummy2
          IF(io_read==0) THEN
           !7/4/2016 - checks of NJOY output reveal that the sum of the P0 Legendre 
           ! constants equals the total cross section for a particular incident energy group
           ! thus, when reading the PKA results into spectra-pka, we need only include
           ! the first value for each i,j pair.
           recoil_kermas(jdummy+1,idummy)=rdummy!+rdummy2
           last_i=idummy
           last_j=jdummy
          ELSE
           ! error or reached end of pka matrix and have read next line
           backspace(pka_unit)
          END IF
         END DO
         
        END IF !flag test 24/2/2014
         
         !6/3/2014 - check cross section case
         ! and skip if not n,g and/or not doing n,g estimate
         IF(INDEX(pka_element,'matrix')==0) THEN
          IF(.not.do_ngamma_estimate) THEN
           last_i=0
           last_j=0
          END IF
          IF(mtd.NE.102) THEN
           last_i=0
           last_j=0
          END IF
         END IF
         

         IF((last_i+last_j).GT.0) THEN
          redo=.false.
         ELSE
          !deallocate arrays and try again
           DEALLOCATE(pka_recoil_energies,pka_incident_energies,recoil_kermas)
           IF(do_tdam) DEALLOCATE(tdam_energies)
           redo=.true.
         END IF
         io_read=0

         
        CASE(3) !GEANT approach

         last_i=0
         last_j=0
         recoil_kermas=0._DBL
        !IF(INDEX(pka_element,'matrix')==0) THEN
         !7/9/2017 not an option
        !ELSE ! normal  read

         DO WHILE (io_read==0)
          ! zero and first order legendre polynomials only
          READ(pka_unit,*,IOSTAT=io_read) idummy,jdummy,rdummy
          IF(io_read==0) THEN
           recoil_kermas(jdummy+1,idummy)=rdummy
           last_i=idummy
           last_j=jdummy
          ELSE
           ! error or reached end of pka matrix and have read next line
           backspace(pka_unit)
          END IF
         END DO
        !END IF      
         
         IF((last_i+last_j).GT.0) THEN
          redo=.false.
         ELSE
          !deallocate arrays and try again
           DEALLOCATE(pka_recoil_energies,pka_incident_energies,recoil_kermas)
           IF(do_tdam) DEALLOCATE(tdam_energies)
           redo=.true.
         END IF
         io_read=0
       
        CASE DEFAULT
         IF(io_read==0) READ(pka_unit,'(8(E10.3))',IOSTAT=io_read) &
                  ((recoil_kermas(j,i),j=1,num_pka_recoil_points),i=1,num_pka_points)        
         redo=.false.
        END SELECT

        IF(io_read.NE.0) THEN
         PRINT *,'error reading pka file'
	 WRITE(log_unit,*)  'error reading pka file'
         io_quit=1
         at_end=.true.
        END IF
        
       ELSE
        PRINT *,'problem reading more of pka file'
        PRINT *,'either we have reached the end (likely) or there is an error'
        PRINT *,'skipping to next pka file if required'
	
        WRITE(log_unit,*) 'problem reading more of pka file'
        WRITE(log_unit,*) 'either we have reached the end (likely) or there is an error'
        WRITE(log_unit,*) 'skipping to next pka file if required'	
	
        !io_quit=1
        at_end=.true.
        redo=.false.
       END IF
     END DO



    

   
       !DO i=1,num_pka_incident_energies
       ! write(202,*) i,pka_incident_energies(i)
       !END DO
       !DO i=1,num_pka_recoil_points
       ! write(201,*) i,pka_recoil_energies(i)
       !END DO
       !STOP
 END SUBROUTINE read_pka
