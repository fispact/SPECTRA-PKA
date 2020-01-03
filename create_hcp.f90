!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk,                                                     *
!******************************************************************************
!end rubric

      SUBROUTINE create_hcp()
      use configs
      IMPLICIT NONE
      integer :: islx
      integer :: isly,islz
      
      INTEGER, allocatable :: id(:)
      integer :: ijk
      REAL (KIND=DBL) :: ww(3),pi
      
      pi=atan2(0._DBL,-1._DBL)

      !assume same number of units in each direction
      ! i.e. not square in hcp case
      islx=box_nunits
      isly=islx
      islz=islx
    !WRITE(*,*) islx,latt,lx
    lx(1)=REAL(islx,DBL)*latt
    lx(2)=REAL(isly,DBL)*latt*cos(pi/6.0_DBL)
    lx(3)=REAL(islz,DBL)*latt*sqrt(3.0_DBL)
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
                x(ijk,1)=latt*REAL(i,DBL)-latt*(REAL(j,DBL))*sin(pi/6._DBL)
                x(ijk,2)=latt*REAL(j,DBL)*cos(pi/6._DBL)
                x(ijk,3)=latt*sqrt(3.0_DBL)*REAL(k,DBL)
                ! tests
                !CALL define_atom_position_hcp(ijk,ww)
                !write(111,*) ijk,i,j,k,x(ijk,:)
           

! in a hcp lattice there is only one more atom in each unit cell:
!  at (2/3,1/3,0.5)
                ijk=ijk+1
                id(ijk)=ijk 
                x(ijk,1)=latt*(REAL(i,DBL)+2._DBL/3._DBL)-&
		    latt*(REAL(j,DBL)+1._DBL/3._DBL)*sin(pi/6._DBL)
                x(ijk,2)=latt*(REAL(j,DBL)+1._DBL/3._DBL)*cos(pi/6._DBL)
                x(ijk,3)=latt*sqrt(3.0_DBL)*(REAL(k,DBL)+0.5_DBL)
                !CALL define_atom_position_hcp(ijk,ww)
                !write(111,*) ijk,i,j,k,x(ijk,:)
            end do
        end do
    end do

END SUBROUTINE create_hcp 


      SUBROUTINE define_atom_position_hcp(nn,xx)
      use configs
      IMPLICIT NONE
      integer, INTENT(in) :: nn ! atom number
      integer :: isly,islz,islx
      REAL (KIND=DBL), intent(out) :: xx(3) 
      REAL (KIND=DBL) :: pi
      integer :: ijk
      
      pi=atan2(0._DBL,-1._DBL)

      

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
    xx(1)=latt*(REAL(i,DBL)+REAL(mod(nn-1,2),DBL)*2._DBL/3._DBL)-&
		    latt*(REAL(j,DBL)+REAL(mod(nn-1,2),DBL)*1_DBL/3._DBL)*sin(pi/6._DBL)
		    
    xx(2)=latt*(REAL(j,DBL)+REAL(mod(nn-1,2),DBL)*1._DBL/3._DBL)*cos(pi/6._DBL)
    
    xx(3)=latt*sqrt(3.0_DBL)*(REAL(k,DBL)+REAL(mod(nn-1,2),DBL)*0.5_DBL)
    
    !write(111,*) nn,i,j,k,xx


END SUBROUTINE define_atom_position_hcp
