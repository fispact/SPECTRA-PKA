SUBROUTINE create_configs
  use globals
  implicit none
  INTEGER :: io_quit_config,natoms
  
  REAL(KIND=DBL), allocatable :: x(:,:)
  REAL(KIND=DBL) :: lx(3),ireal
  INTEGER, ALLOCATABLE :: atom_mass(:)
  CHARACTER, ALLOCATABLE :: atom_element(:)
  integer, allocatable, dimension(:) :: seed_array
  INTEGER :: jcount,kcount,choice
  
  io_quit_config=0
  x=0._DBL
  
  PRINT *,'creating evolving atomic composition configurations'

    call random_seed(size=getsize)
    allocate (seed_array(getsize))
    seed_array=time() 
    call random_seed(put=seed_array) 
    deallocate (seed_array)
  
  IF(box_type==1) THEN
   CALL create_bcc(x,box_nunits,natoms,latt,lx)
   
  ELSEIF (box_type==2) THEN
   !CALL create_fcc(x,natoms)
  ELSE
    print *,'not a valid lattice type, skipping cfg creation'
    io_quit_config=1
  END IF
  
  ALLOCATE(atom_mass(natoms),atom_element(natoms))
  atom_mass=0
  ! allocate elements by input composition ratios
  
  kcount=0
  print *,'Starting configuration:'
  DO i=1,number_pka_files
   !parent_ele(1,numner_pka_files)
   !parent_num
   !pka_ratios
   jcount=0
   DO
    IF(jcount.GE.NINT(REAL(natoms,DBL)*pka_ratios(i))) exit
    IF(kcount.GE.natoms) exit
    CALL random_number(ireal)
    choice=NINT(ireal*REAL(natoms,DBL)
    IF((choice.LE.0).OR.(choice.GT.natoms)) cycle ! shouldn't happen
    IF(atom_mass(choice).NE.0) cycle ! already assigned
    atom_mass(choice)=parent_num(i)
    atom_element(choice)=parent_ele(i)
    jcount=jcount+1
    kcount=kcount+1
   END DO
   print *,jcount,'  ',parent_ele(i),parent_num(i),' atoms'
  END DO
  print *,kcount,' out of ',natoms,' assigned'
  
  
  isteps=0
  time=0._DBL
  DO
   IF(isteps.GT.nsteps) EXIT
   
 ! how to do this:
 ! call make config and intial composition at the start of code
 ! then after reading each channel, evolve a sequence of compositions according to it????
 ! since interactions can only happen on atoms of the original selection
 ! this shouldn't create bias towards later channels read
 ! in any case, this is only a short-timescale, low transmutation rate approximation
 ! for longer timescales or if transmutation is high, we should loop over configurations
 ! provided by FISPACT-II
 ! something for the future.  
 ! actually base on isotope (or element) sums at end of execution. - but then won't now parent to daughter conversion?
 ! PKAs based on total rate, probability of energy based on total PKA distribution
 ! and then type based on ratio of different isotopes at each 
 
 ! probability per timestep:
 ! total PKAS per channel=m in PKAs/s/box [number of atoms]
 ! number PKAs/timestep=mDelta{T}=n
 ! INT(n) occur per timestep + n-INT(n) occurrences based on fractional random number
 ! [energy based on random selection from a normalised distribution]
 ! loop over all time steps for each channel
 ! select atoms at random - can we create subsets of each atomic species to speed up selection? [but this would change with time]
   
   
   
   isteps=isteps+1
   time=time+timestep
  
  END DO
  
  
  
  





END create_configs
