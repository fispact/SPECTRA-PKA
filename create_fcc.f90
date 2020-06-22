!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk
!* copyright, 2015 (first version), 2018 (first git repository), UKAEA                                                     *
!******************************************************************************
!end rubric

      SUBROUTINE create_fcc()
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
    natoms=4*islx*isly*islz 
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
                !CALL define_atom_position_fcc(ijk,ww)
                !write(111,*) ijk,i,j,k,x(ijk,:)
           
! in a fcc lattice there arw 3 more atoms in each unit cell:
! the body-centre atom at (0.5,0.5,0),(0,0.5,0.5),(0.5,0,0.5)
                ijk=ijk+1
                id(ijk)=ijk 
                x(ijk,1)=latt*(REAL(i,DBL)+0.5_DBL)
                x(ijk,2)=latt*(REAL(j,DBL)+0.5_DBL)
                x(ijk,3)=latt*(REAL(k,DBL))
                !CALL define_atom_position_fcc(ijk,ww)
                !write(111,*) ijk,i,j,k,x(ijk,:)

                ijk=ijk+1
                id(ijk)=ijk 
                x(ijk,1)=latt*(REAL(i,DBL))
                x(ijk,2)=latt*(REAL(j,DBL)+0.5_DBL)
                x(ijk,3)=latt*(REAL(k,DBL)+0.5_DBL)
                !CALL define_atom_position_fcc(ijk,ww)
                !write(111,*) ijk,i,j,k,x(ijk,:)
		
                ijk=ijk+1
                id(ijk)=ijk 
                x(ijk,1)=latt*(REAL(i,DBL)+0.5_DBL)
                x(ijk,2)=latt*(REAL(j,DBL))
                x(ijk,3)=latt*(REAL(k,DBL)+0.5_DBL)
                !CALL define_atom_position_fcc(ijk,ww)
                !write(111,*) ijk,i,j,k,x(ijk,:)				
		
            end do
        end do
    end do

END SUBROUTINE create_fcc 


      SUBROUTINE define_atom_position_fcc(nn,xx)
      use configs
      IMPLICIT NONE
      integer, INTENT(in) :: nn ! atom number
      integer :: isly,islz,islx,mm
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
    !natoms=4*islx*isly*islz

    ijk=nn
    IF(mod(nn,4)==0) ijk=nn-1
    
    i=INT(REAL(ijk,DBL)/(4._DBL*REAL(isly,DBL)*REAL(islz,DBL)))
    j=INT((REAL(ijk,DBL)-REAL(i,DBL)*4._DBL*REAL(isly,DBL)*REAL(islz,DBL))/ &
              (4._DBL*REAL(islz,DBL)))
    k=INT((REAL(ijk,DBL)-REAL(i,DBL)*4._DBL*REAL(isly,DBL)*REAL(islz,DBL)-&
        REAL(j,DBL)*4._DBL*REAL(islz,DBL))/4._DBL)
    
    mm=mod(nn-1,4)
    SELECT CASE(mm)
     CASE(0)
     xx(1)=latt*(REAL(i,DBL))
     xx(2)=latt*(REAL(j,DBL))
     xx(3)=latt*(REAL(k,DBL))
     CASE(1)
     xx(1)=latt*(REAL(i,DBL)+REAL(mm,DBL)*0.5_DBL)
     xx(2)=latt*(REAL(j,DBL)+REAL(mm,DBL)*0.5_DBL)
     xx(3)=latt*(REAL(k,DBL))
     CASE(2)
     xx(1)=latt*(REAL(i,DBL))
     xx(2)=latt*(REAL(j,DBL)+REAL(mm,DBL)*0.5_DBL)
     xx(3)=latt*(REAL(k,DBL)+REAL(mm,DBL)*0.5_DBL)
     CASE(3)
     xx(1)=latt*(REAL(i,DBL)+REAL(mm,DBL)*0.5_DBL)
     xx(2)=latt*(REAL(j,DBL))
     xx(3)=latt*(REAL(k,DBL)+REAL(mm,DBL)*0.5_DBL)
    END SELECT               
    !write(111,*) nn,i,j,k,xx



END SUBROUTINE define_atom_position_fcc



     
