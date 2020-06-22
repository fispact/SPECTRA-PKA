!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk
!* copyright, 2015 (first version), 2018 (first git repository), UKAEA                                                     *
!******************************************************************************
!end rubric

      SUBROUTINE create_bcc()
      use configs
      IMPLICIT NONE
      integer :: islx
      integer :: isly,islz
      
      INTEGER, allocatable :: id(:)
      integer :: ijk
      REAL (KIND=DBL) :: ww(3)

      islx=box_nunits
      isly=islx
      islz=islx
    !WRITE(*,*) islx,latt,lx
    lx(1)=REAL(islx,DBL)*latt
    lx(2)=REAL(isly,DBL)*latt
    lx(3)=REAL(islz,DBL)*latt
    natoms=2*islx*isly*islz 
    !WRITE(*,*) natoms
    !WRITE(*,*) islx,latt,lx

    allocate(id(natoms),x(natoms,3),xe(natoms,3))

    ijk=0
    do i=0,islx-1
        do j=0,isly-1
            do k=0,islz-1
                ijk=ijk+1
                id(ijk)=ijk
                x(ijk,1)=latt*REAL(i,DBL)
                x(ijk,2)=latt*REAL(j,DBL)
                x(ijk,3)=latt*REAL(k,DBL)
                ! tests
                !CALL define_atom_position_bcc(ijk,ww)
                !write(111,*) ijk,i,j,k,x(ijk,:)
           
! in a bcc lattice there is only one more atom in each unit cell:
! the body-centre atom at (0.5,0.5,0.5)
                ijk=ijk+1
                id(ijk)=ijk 
                x(ijk,1)=latt*(REAL(i,DBL)+0.5_DBL)
                x(ijk,2)=latt*(REAL(j,DBL)+0.5_DBL)
                x(ijk,3)=latt*(REAL(k,DBL)+0.5_DBL)
                !CALL define_atom_position_bcc(ijk,ww)
                !write(111,*) ijk,i,j,k,x(ijk,:)
            end do
        end do
    end do

END SUBROUTINE create_bcc 


      SUBROUTINE define_atom_position_bcc(nn,xx)
      use configs
      IMPLICIT NONE
      integer, INTENT(in) :: nn ! atom number
      integer :: isly,islz,islx
      REAL (KIND=DBL), intent(out) :: xx(3) 
      

      integer :: ijk

      islx=box_nunits
      isly=islx
      islz=islx
 
    IF (nn>natoms) THEN
     PRINT *,'atom id outside of bounds, quiting',nn,natoms
     STOP
    END IF
    !WRITE(*,*) natoms
    !WRITE(*,*) islx,latt,lx
    !natoms=2*islx*isly*islz

    ijk=nn
    IF(mod(nn,2)==0) ijk=nn-1
    i=INT(REAL(ijk,DBL)/(2._DBL*REAL(isly,DBL)*REAL(islz,DBL)))
    j=INT((REAL(ijk,DBL)-REAL(i,DBL)*2._DBL*REAL(isly,DBL)*REAL(islz,DBL))/ &
              (2._DBL*REAL(islz,DBL)))
    k=INT((REAL(ijk,DBL)-REAL(i,DBL)*2._DBL*REAL(isly,DBL)*REAL(islz,DBL)-&
        REAL(j,DBL)*2._DBL*REAL(islz,DBL))/2._DBL)
    xx(1)=latt*(REAL(i,DBL)+REAL(mod(nn-1,2),DBL)*0.5_DBL)
    xx(2)=latt*(REAL(j,DBL)+REAL(mod(nn-1,2),DBL)*0.5_DBL)
    xx(3)=latt*(REAL(k,DBL)+REAL(mod(nn-1,2),DBL)*0.5_DBL)
    !write(111,*) nn,i,j,k,xx


END SUBROUTINE define_atom_position_bcc


