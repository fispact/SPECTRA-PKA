!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk,                                                     *
!******************************************************************************
!end rubric

 


SUBROUTINE write_cfg_start(filename,natoms,lxx,io_quit,do_velocity)
 use accuracy
 use file_units
 IMPLICIT NONE
 
 INTEGER, intent(in) :: natoms
 REAL(KIND=DBL), intent(in) :: lxx(3)
 CHARACTER (LEN=500), intent(in) :: filename
 INTEGER, intent(inout) :: io_quit
 LOGICAL :: do_velocity
 CHARACTER (LEN=300) outstr 
 
 io_quit=0


 print *,TRIM(ADJUSTL(filename))
 OPEN(unit=cfg_file,file=TRIM(ADJUSTL(filename)),STATUS='REPLACE',IOSTAT=io_quit)
 IF(io_quit==0) THEN
 
    WRITE (cfg_file,'(A22,I10)') 'Number of particles = ',natoms
    !WRITE (cfg_file,*) 'A = 1.0 Angstrom (basic length scale)'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(1,1) = ',lxx(1),'A'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(1,2) = ',0._DBL,'A'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(1,3) = ',0._DBL,'A'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(2,1) = ',0._DBL,'A'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(2,2) = ',lxx(2),'A'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(2,3) = ',0._DBL,'A'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(3,1) = ',0._DBL,'A'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(3,2) = ',0._DBL,'A'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(3,3) = ',lxx(3),'A'
    IF(.not.do_velocity) THEN
     WRITE (cfg_file,'(A13)') '.NO_VELOCITY.'   
     WRITE (cfg_file,'(A18)') 'entry_count =    3'
    ELSE
     WRITE (cfg_file,'(A18)') 'entry_count =    6'
    END IF
 

  END IF
 
  

END SUBROUTINE write_cfg_start

SUBROUTINE write_lammps_start(filename,natoms,lxx,io_quit)
 use accuracy
 use file_units
 IMPLICIT NONE
 
 INTEGER, intent(in) :: natoms
 REAL(KIND=DBL), intent(in) :: lxx(3)
 CHARACTER (LEN=500), intent(in) :: filename
 INTEGER, intent(inout) :: io_quit
 

 INTEGER :: noutput

 
 CHARACTER (LEN=100) outstr 
 
 io_quit=0

 
 OPEN(unit=lammps_file,file=TRIM(ADJUSTL(filename)),STATUS='REPLACE',IOSTAT=io_quit)
 IF(io_quit==0) THEN


    WRITE (lammps_file,'(A)') 'lammps configuration output file'
    WRITE (lammps_file,'(I10,A10)') natoms,'atoms'
    WRITE(lammps_file,'(2e16.8,A15)') 0d0,lxx(1),' xlo xhi'
    WRITE(lammps_file,'(2e16.8,A15)') 0d0,lxx(2),' ylo yhi'
    WRITE(lammps_file,'(2e16.8,A15)') 0d0,lxx(3),' zlo zhi'
    WRITE(lammps_file,'(3e16.8,A15)') 0d0,0d0,0d0,'xy xz yz'
    
    

  END IF
 
 

END SUBROUTINE write_lammps_start



SUBROUTINE write_cfg(filename,xx,natoms,lxx,io_quit,ele_type,ele_mass,xee)
 use accuracy
 use file_units
 IMPLICIT NONE
 
 INTEGER, intent(in) :: natoms
 REAL(KIND=DBL), intent(in) :: lxx(3),xx(natoms,3),xee(natoms,3)
 CHARACTER (LEN=500), intent(in) :: filename
 INTEGER, intent(inout) :: io_quit
 LOGICAL :: output(natoms)
 CHARACTER (LEN=2), intent(in) :: ele_type(natoms)
 INTEGER, intent(in) :: ele_mass(natoms)
 INTEGER :: noutput
 CHARACTER (LEN=2) :: current_ele
 LOGICAL :: newblock
 INTEGER :: i,currentistart
 
 CHARACTER (LEN=300) outstr 
 
 io_quit=0
 output=.false.
 !DO i=1,natoms
 ! PRINT *,i,ele_type(i),ele_mass(i)
 !END DO
 print *,TRIM(ADJUSTL(filename))
 OPEN(unit=cfg_file,file=TRIM(ADJUSTL(filename)),STATUS='REPLACE',IOSTAT=io_quit)
 IF(io_quit==0) THEN
 
    WRITE (cfg_file,'(A22,I10)') 'Number of particles = ',natoms
    !WRITE (cfg_file,*) 'A = 1.0 Angstrom (basic length scale)'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(1,1) = ',lxx(1),'A'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(1,2) = ',0._DBL,'A'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(1,3) = ',0._DBL,'A'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(2,1) = ',0._DBL,'A'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(2,2) = ',lxx(2),'A'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(2,3) = ',0._DBL,'A'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(3,1) = ',0._DBL,'A'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(3,2) = ',0._DBL,'A'
    WRITE (cfg_file,'(A10,ES25.10,A2)') 'H0(3,3) = ',lxx(3),'A'
    ! assume velocities will be given
    !WRITE (cfg_file,'(A13)') '.NO_VELOCITY.'   
    WRITE (cfg_file,'(A18)') 'entry_count =    6'
    noutput=0
    DO 
     IF(noutput.GE. natoms) EXIT
     IF(io_quit.NE.0) EXIT
     newblock=.true.
     DO i=1,natoms
       IF(output(i)) cycle
       IF(io_quit.NE.0) EXIT
       IF(newblock) THEN
         WRITE(cfg_file,'(I3)') ele_mass(i)
	 WRITE(cfg_file,'(A2)') ele_type(i)
	 newblock=.false.
	 currentistart=i
       END IF
       !PRINT *,i,currentistart,ele_type(i),ele_type(currentistart),ele_mass(i),ele_mass(currentistart)
       IF((ele_mass(i)==ele_mass(currentistart)).AND. &
          (TRIM(ADJUSTL(ele_type(i)))==TRIM(ADJUSTL(ele_type(currentistart))))) THEN
	  WRITE(outstr,'(1x,6(1x,e16.8))') &
	    xx(i,1)/lxx(1),xx(i,2)/lxx(2),xx(i,3)/lxx(3),&
	     xee(i,1)/lxx(1),xee(i,2)/lxx(2),xee(i,3)/lxx(3)
	  WRITE(cfg_file,'(a)',IOSTAT=io_quit) TRIM(ADJUSTL(outstr))
	  output(i)=.true.
	  noutput=noutput+1
       END IF
     END DO
     
    END DO
    IF(io_quit==0) CLOSE(cfg_file)
  END IF
 
  

END SUBROUTINE write_cfg
     

SUBROUTINE write_lammps(filename,xx,natoms,lxx,io_quit,ele_type,ele_mass,xee)
 use accuracy
 use file_units
 IMPLICIT NONE
 
 INTEGER, intent(in) :: natoms
 REAL(KIND=DBL), intent(in) :: lxx(3),xx(natoms,3),xee(natoms,3)
 CHARACTER (LEN=500), intent(in) :: filename
 INTEGER, intent(inout) :: io_quit
 
 INTEGER :: output(natoms),to_output(natoms)
 CHARACTER (LEN=2), intent(in) :: ele_type(natoms)
 INTEGER, intent(in) :: ele_mass(natoms)
 INTEGER :: noutput
 CHARACTER (LEN=2) :: current_ele
 LOGICAL :: newblock
 INTEGER :: i,currentistart,atom_types,j
 
 CHARACTER (LEN=100) outstr 
 
 io_quit=0
 !DO i=1,natoms
 ! PRINT *,i,ele_type(i),ele_mass(i)
 !END DO
 
 OPEN(unit=lammps_file,file=TRIM(ADJUSTL(filename)),STATUS='REPLACE',IOSTAT=io_quit)
 IF(io_quit==0) THEN


    WRITE (lammps_file,'(A)') 'lammps configuration output file'
    WRITE (lammps_file,'(I10,A10)') natoms,'atoms'
    WRITE(lammps_file,'(2e16.8,A15)') 0d0,lxx(1),' xlo xhi'
    WRITE(lammps_file,'(2e16.8,A15)') 0d0,lxx(2),' ylo yhi'
    WRITE(lammps_file,'(2e16.8,A15)') 0d0,lxx(3),' zlo zhi'
    WRITE(lammps_file,'(3e16.8,A15)') 0d0,0d0,0d0,'xy xz yz'
    
    
! count atom types
atom_types=1
to_output(1)=1
outer: DO i=1,natoms
        DO j=1,atom_types
         IF((ele_mass(i)==ele_mass(to_output(j))).AND. &
          (TRIM(ADJUSTL(ele_type(i)))==TRIM(ADJUSTL(ele_type(to_output(j)))))) THEN
          output(i)=j
          cycle outer
         END IF
        END DO
        ! if we get here then we haven't found a match
          atom_types=atom_types+1
          to_output(atom_types)=i
          output(i)=atom_types
      END DO outer


    WRITE(lammps_file,'(I5,A15)') atom_types,' atom types'
    WRITE(lammps_file,'(A10)') 'Masses'
    WRITE(lammps_file,*)
    DO i=1,atom_types
     WRITE(lammps_file,'(I5,I10)') i,ele_mass(to_output(i))
    END DO
    
    WRITE(lammps_file,'(A)') 'Atoms'
    WRITE(lammps_file,*)
    do i=1,natoms
     WRITE(lammps_file,'(I10,I5,3e16.8,A15)') i,output(i),xx(i,1),xx(i,2),xx(i,3)
    end do    
    
    WRITE(lammps_file,'(A)') 'Velocities'
    WRITE(lammps_file,*)
    do i=1,natoms
     WRITE(lammps_file,'(I10,3e16.8)') i,xee(i,1),xee(i,2),xee(i,3)
    end do  

  END IF
 
IF(io_quit==0) CLOSE(lammps_file,IOSTAT=io_quit)  

END SUBROUTINE write_lammps
     
