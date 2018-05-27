!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk,                                                     *
!******************************************************************************
!end rubric

      SUBROUTINE create_bcc(xx,islx,natoms,latt,lxx)
      use accuracy
      IMPLICIT NONE
      REAL(KIND=DBL), allocatable, intent(inout) :: xx(:,:)
      integer, intent(in) :: islx
      integer :: isly,islz
      integer, intent(out) :: natoms
      REAL(KIND=DBL), intent(out) :: lxx(3)
      REAL(KIND=DBL), intent(in) :: latt
      
      INTEGER, allocatable :: id(:)
      integer :: i,j,k,ijk

      
      isly=islx
      islz=islx

    lxx(1)=REAL(islx,DBL)*latt
    lxx(2)=REAL(isly,DBL)*latt
    lxx(3)=REAL(islz,DBL)*latt
    natoms=2*islx*isly*islz 
    allocate(id(natoms),xx(natoms,3))
    ijk=0
    do i=0,islx-1
        do j=0,isly-1
            do k=0,islz-1
                ijk=ijk+1
                id(ijk)=ijk
                xx(ijk,1)=latt*REAL(i,DBL)
                xx(ijk,2)=latt*REAL(j,DBL)
                xx(ijk,3)=latt*REAL(k,DBL)
           
! in a bcc lattice there is only one more atom in each unit cell:
! the body-centre atom at (0.5,0.5,0.5)
                ijk=ijk+1
                id(ijk)=ijk 
                xx(ijk,1)=latt*(REAL(i,DBL)+0.5_DBL)
                xx(ijk,2)=latt*(REAL(j,DBL)+0.5_DBL)
                xx(ijk,3)=latt*(REAL(k,DBL)+0.5_DBL)
              
            end do
        end do
    end do
  
END SUBROUTINE create_bcc 

SUBROUTINE write_cfg(filename,xx,natoms,lxx,io_quit,ele_type,ele_mass)
 use accuracy
 IMPLICIT NONE
 
 REAL(KIND=DBL), intent(in) :: lxx(3),xx(natoms,3)
 CHARACTER (LEN=500) :: intent(in) :: filename
 INTEGER, intent(inout) :: io_quit
 INTEGER, PARAMETER :: cfg_file=1111
 LOGICAL :: output(natoms)
 CHARACTER (LEN=2), intent(in) :: ele_type(natoms)
 INTEGER, intent(in) :: ele_mass(natoms)
 INTEGER :: noutput
 CHARACTER (LEN=2) :: current_ele
 LOGICAL :: newblock
 INTEGER :: i,currentistart
 CHARACTER (LEND=100) outstr 
 
 io_quit=0
 output=.false.
 
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
    WRITE (cfg_file,'(A13)') '.NO_VELOCITY.'   
 
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
       IF((ele_mass(i)==ele_mass(currentistart)).AND. &
          (ele_type(i)==ele_type(currentistart))) THEN
	  WRITE(outstr,'(1x,3(1x,e16.8))') &
	    xx(i,1)/lxx(1),xx(i,2)/lxx(2),xx(i,3)/lxx(3)
	  WRITE(cfg_file,'(a)',IOSTAT=io_quit)
	  output(i)=.true.
	  noutput=noutput+1
       END IF
     END DO
     IF(io_quit==0) CLOSE(cfg_file,IOSTAT=io_quit)
  END IF
 
  

END SUBROUTINE write_cfg()
     
     
