!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk
!* copyright, 2015 (first version), 2018 (first git repository), UKAEA                                                     *
!******************************************************************************
!end rubric

 SUBROUTINE read_flux()
 use globals
 IMPLICIT NONE
 CHARACTER (LEN=200) :: ddstr
 
 OPEN(unit=flux_unit,FILE=TRIM(ADJUSTL(flux_filename)),STATUS='OLD',IOSTAT=io_open)
  !print *,TRIM(ADJUSTL(flux_filename)),io_open
  IF(io_open==0) THEN
  
    READ(flux_unit,'(A100)',IOSTAT=io_read) flux_title
        IF(io_read==0) READ(flux_unit,*,IOSTAT=io_read) itype,isigma,igroup,ipka,acnm,irrtime
    
    IF(io_read==0) THEN
     ! just consider cases for example files for now
     if((itype.GT.1).AND.(io_read==0)) THEN
      READ(flux_unit,*,IOSTAT=io_read) number_flux_groups,ksail      
     END IF
     IF(io_read==0) THEN
       number_flux_ebins=number_flux_groups+1
       !12/3/2018 - over allocated flux arrays (last bin never becomes non-zero)
       ALLOCATE(flux_energies(number_flux_ebins),fluxes(number_flux_ebins),&
              flux_ebin_widths(number_flux_ebins),fluxes_norm(number_flux_ebins))
       fluxes=0._DBL
       READ(flux_unit,*,IOSTAT=io_read) (flux_energies(i),i=1,number_flux_ebins)
       IF(io_read==0) &
        READ(flux_unit,*,IOSTAT=io_read) (fluxes(i),i=1,number_flux_groups) 
     END IF
    END IF
     IF(io_read==0) THEN     
      IF(ksail.GT.0) THEN
       ALLOCATE(flux_covariances(number_flux_groups+1,number_flux_groups+1))
       i=1
       DO WHILE((io_read==0).AND.(i.LE.number_flux_groups))
        READ(flux_unit,*,IOSTAT=io_read) (flux_covariances(i,j),j=i,number_flux_groups)
        i=i+1
       END DO
       DO i=2,number_flux_groups
        DO j=1,i-1
         flux_covariances(i,j)=flux_covariances(j,i)
        END DO
       END DO
      ELSE IF(ksail.LT.0) THEN
       !no covariances
      END IF
     END IF
    IF(io_read.NE.0) THEN
     PRINT *,'error reading from flux input file'
     io_quit=1
    END IF
    CLOSE(flux_unit,STATUS='KEEP')
    flux_energies=flux_energies*flux_energy_rescale_value
    fluxes=fluxes*flux_rescale_value
 ELSE
   PRINT *,'unable to open flux input file'
   io_quit=1
 END IF
 
 END SUBROUTINE read_flux
