MODULE configs
  use globals

  REAL(KIND=DBL), allocatable :: x(:,:),xe(:,:)
  REAL(KIND=DBL) :: lx(3)
  INTEGER, ALLOCATABLE :: atom_mass(:)
  CHARACTER (LEN=2), ALLOCATABLE :: atom_element(:)
  integer, allocatable, dimension(:) :: seed_array
  integer :: io_quit_config,getsize
  integer (kind=8) :: natoms
  INTEGER, ALLOCATABLE :: initial_atom_types_count(:)
  INTEGER, ALLOCATABLE ::pka_events_id(:)
  REAL (KIND=DBL), ALLOCATABLE :: pka_events_vec(:,:),pka_events_ke(:),pka_min_sep(:,:)
  INTEGER, ALLOCATABLE :: pka_events_mass(:),pka_events_step(:)
  CHARACTER (LEN=2), ALLOCATABLE :: pka_events_ele(:)
  integer,dimension(8) :: timevalues
  

END MODULE configs


SUBROUTINE create_configs
  use configs
  implicit none
  

  INTEGER :: jcount,kcount,choice,num_events,num_events_previous
  INTEGER :: isteps,ii,jj,ichannel,num_PKAs_step,ienergy,ipkas,jenergy
  REAL(KIND=DBL) :: cumtime,ireal,remain,ke,jrealvector(3),irealvector(3)
  logical :: found,bca_flag
  CHARACTER (LEN=500) :: outstr,tempstr,tempstr2
  INTEGER :: ievent,ifiles,foundi,jevent
  CHARACTER (LEN=30) :: date
  INTEGER(8) ::itime,num_above_threshold
  REAL(KIND=DBL) :: total_channelpkas
  INTEGER, ALLOCATABLE ::temp_pka_events_id(:)
  REAL (KIND=DBL), ALLOCATABLE :: temp_pka_events_vec(:,:),temp_pka_events_ke(:),temp_pka_min_sep(:,:)
  INTEGER, ALLOCATABLE :: temp_pka_events_mass(:),temp_pka_events_step(:)
  CHARACTER (LEN=2), ALLOCATABLE :: temp_pka_events_ele(:)
  
  
  !16/4/19 - separation tracking
  REAL (KIND=DBL), ALLOCATABLE :: total_sep_dist(:),total_seps(:),event_dist(:)
  REAL (KIND=DBL) :: current_sep
  
  



  
  PRINT *,'Analyzing time dependent distribution of PKAs at the atomic lattice scale'
  PRINT *,'using the following ',config_num_pka_vectors,' dominant reaction channels'
  DO jj=1,config_num_pka_vectors
    PRINT *,jj,config_parent_eles(jj),config_parent_nums(jj),config_daughter_eles(jj),config_daughter_nums(jj),config_pka_strings(jj)
    
  
  END DO
 
 
 IF(do_output_configs) PRINT *,'and creating evolving atomic composition configurations'
 
 ! we need to start be defining some array to hold events (basically hit atoms)
 ! in the limit of low PKA density there won't be many
 ! but the coverage will increase in time.
 ! lets say sqrt(natoms)
 ! need atom ids, atom velocities and atom types/mass
 !15/4/2019 - no lets do this on the fly
 !ii=INT(sqrt(REAL(natoms)))
 !ALLOCATE(pka_events_id(ii),pka_events_vec(ii,3), &
 !         pka_events_ele(ii),pka_events_mass(ii), &
 !         pka_events_step(ii),pka_events_ke(ii))
 !max_events=ii
  
  
  isteps=1
  cumtime=timestep
  num_events=0
  num_events_previous=0
  IF (do_bca) CALL bca_setup()
  
    call random_seed(size=getsize)
    allocate (seed_array(getsize))
    call date_and_time(VALUES=timevalues)
    seed_array=SUM(timevalues)+time() 
    !print *,'time',time(),SUM(timevalues)
    
    call random_seed(put=seed_array) 
    deallocate (seed_array) 
    
    
  OPEN(unit=pka_events,file='config_events.pka',STATUS='REPLACE',IOSTAT=io_open)
  IF(io_open.NE.0) THEN
   PRINT *,'error opening pka events out file, quiting config analysis'
   return
  END IF
  OPEN(unit=pka_analysis,file='config_analysis.pka',STATUS='REPLACE',IOSTAT=io_open)  
  IF(io_open.NE.0) THEN
   PRINT *,'error opening pka events analysis out file, quiting config analysis'
   return
  END IF
  outstr=''
  
  
  
  DO i=1,totalglobal_num_pka_recoil_points_master
    write(tempstr,'(ES15.3)') totalglobal_pka_recoil_energies_master(i)
    write(tempstr2,'(A20)') '  >'//TRIM(ADJUSTL(tempstr))//'(MeV)'
    outstr=TRIM(outstr)//TRIM(tempstr2)
  END DO
  write(pka_analysis,'(A1,A19,4A20,A40,A)') '#','time step num','ave. sep. dist.', &
      'total pairs','ave. min. sep.', 'min. sep.','event pair of min. sep.',TRIM(outstr)
  write(pka_events,'(x,A19,3A20,3A20,A20,4A20)') &
          'PKA atom','PKA-position-x','PKA-position-y','pka-position-z', &
          'pka-velocity-x','pka-velocity-y','pka-velocity-z', &
          'PKA element','PKA mass','time step num',&
          'PKA-energy (MeV)','cumulative time'
     
  ALLOCATE(total_sep_dist(totalglobal_num_pka_recoil_points_master+1), &
           total_seps(totalglobal_num_pka_recoil_points_master+1), &
           event_dist(totalglobal_num_pka_recoil_points_master+1))
  total_sep_dist=0._DBL
  total_seps=0._DBL
  event_dist=0._DBL
  
  
   DO ichannel=1,config_num_pka_vectors
    total_channelpkas=SUM(config_pka_vectors(ichannel,1:config_global_num_pka_recoil_points))
    print *,'Channel ',ichannel,': PKAs per atom per timestep: ',total_channelpkas
   END DO  
  
  DO
   PRINT *,'Time step: ',isteps
   WRITE(log_unit,*) 'Time step: ',isteps 
   num_above_threshold=0
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
 ! actually base on isotope (or element) sums at end of execution. - but then won't know parent to daughter conversion?
 ! PKAs based on total rate, probability of energy based on total PKA distribution
 ! and then type based on ratio of different isotopes at each 
 
 ! probability per timestep:
 ! total PKAS per channel=m in PKAs/s/box [number of atoms]
 ! number PKAs/timestep=mDelta{T}=n
 ! INT(n) occur per timestep + n-INT(n) occurrences based on fractional random number
 ! [energy based on random selection from a normalised distribution]
 ! loop over all time steps for each channel
 ! select atoms at random - can we create subsets of each atomic species to speed up selection? [but this would change with time]




!4/4/4 - might need to do this in abstract from the atomic lattice due to memory constraints.



    
    
    
    ! reset atom energies 
    
   
   DO ichannel=1,config_num_pka_vectors
    !PRINT *,config_daughter_nums(ichannel),config_daughter_eles(ichannel),ichannel
    !15/4/2019 - select correct number of PKAs for channel
    ! and then sample energy dist.
    ! to sort out problem with small numbers
    total_channelpkas=SUM(config_pka_vectors(ichannel,1:config_global_num_pka_recoil_points))
    !print *,'channel ',ichannel,' PKAs per atom per timestep: ',total_channelpkas
    !DO ienergy=1,totalglobal_num_pka_recoil_points_master
    !DO ienergy=11,11
    ! assumes pkas have pkas/s/atom units - will depend on user specifying correct input flux
     ! but we will generally have very small numbers, so lets try total time and divide
     ireal=total_channelpkas*timestep*natoms*nsteps
     num_PKAs_step=INT(ireal/REAL(nsteps))
     remain=ireal/REAL(nsteps)-REAL(INT(ireal/REAL(nsteps)))
!write(700,*) config_pka_vectors(ichannel,:)
!write(700,*) config_global_pka_recoil_energies
     !PRINT *,num_PKAs_step,remain,natoms,timestep
    !stop

     CALL random_number(ireal)
     IF(ireal<remain) THEN ! add an extra one
       num_PKAs_step=num_PKAs_step+1
     END IF
     write(log_unit,*) 'channel ',ichannel,': ',num_PKAs_step,' PKAs required in timestep ',isteps 
     ! increase arrays by required amount
     !print *,allocated(pka_events_id)
!stop
IF(ALLOCATED(pka_events_id))THEN
 IF(ALLOCATED(temp_pka_events_id))THEN
   DEALLOCATE(temp_pka_events_id,temp_pka_events_vec, &
          temp_pka_events_ele,temp_pka_events_mass, &
          temp_pka_events_step,temp_pka_events_ke,temp_pka_min_sep)
 END IF
 ALLOCATE(temp_pka_events_id(num_events),temp_pka_events_vec(num_events,3), &
          temp_pka_events_ele(num_events),temp_pka_events_mass(num_events), &
          temp_pka_events_step(num_events),temp_pka_events_ke(num_events), &
          temp_pka_min_sep(num_events,2*totalglobal_num_pka_recoil_points_master+2))   
 temp_pka_events_id(1:num_events)=pka_events_id(1:num_events)
 temp_pka_events_vec(1:num_events,1:3)=pka_events_vec(1:num_events,1:3)
 temp_pka_events_ele(1:num_events)=pka_events_ele(1:num_events)
 temp_pka_events_mass(1:num_events)=pka_events_mass(1:num_events)
 temp_pka_events_step(1:num_events)=pka_events_step(1:num_events)
 temp_pka_events_ke(1:num_events)=pka_events_ke(1:num_events)
 temp_pka_min_sep(1:num_events,:)=pka_min_sep(1:num_events,:)
 DEALLOCATE(pka_events_id,pka_events_vec, &
          pka_events_ele,pka_events_mass, &
          pka_events_step,pka_events_ke,pka_min_sep)
            
 ALLOCATE(pka_events_id(num_events+num_PKAs_step),pka_events_vec(num_events+num_PKAs_step,3), &
          pka_events_ele(num_events+num_PKAs_step),pka_events_mass(num_events+num_PKAs_step), &
          pka_events_step(num_events+num_PKAs_step),pka_events_ke(num_events+num_PKAs_step), &
          pka_min_sep(num_events+num_PKAs_step,2*totalglobal_num_pka_recoil_points_master+2))
 pka_events_id(1:num_events)=temp_pka_events_id(1:num_events)
 pka_events_vec(1:num_events,1:3)=temp_pka_events_vec(1:num_events,1:3)
 pka_events_ele(1:num_events)=temp_pka_events_ele(1:num_events)
 pka_events_mass(1:num_events)=temp_pka_events_mass(1:num_events)
 pka_events_step(1:num_events)=temp_pka_events_step(1:num_events)
 pka_events_ke(1:num_events)=temp_pka_events_ke(1:num_events)  
 pka_min_sep(1:num_events,:)=temp_pka_min_sep(1:num_events,:)        
ELSE
  ALLOCATE(pka_events_id(num_events+num_PKAs_step),pka_events_vec(num_events+num_PKAs_step,3), &
          pka_events_ele(num_events+num_PKAs_step),pka_events_mass(num_events+num_PKAs_step), &
          pka_events_step(num_events+num_PKAs_step),pka_events_ke(num_events+num_PKAs_step), &
          pka_min_sep(num_events+num_PKAs_step,2*totalglobal_num_pka_recoil_points_master+2)) 
END IF     
     
  !STOP   
     
    

     !cycle
     DO ipkas=1,num_PKAs_step
       found=.false.
       !print *,ipkas
choiceloop:DO WHILE (.not.found) 
        CALL random_number(ireal)
        choice=NINT(ireal*REAL(natoms,DBL))
        !PRINT *,choice,num_events,ipkas,isteps
        ! in a given timestep, we are allowing the same atom to be hit twice
        ! first check if atom has been hit before
        ievent=1
        foundi=num_events+1
        DO 
         IF(ievent.gt.num_events) exit
         !5/4/2019 we are now keeping all events, so atoms may appear more than once
         ! search whole list and take last state.
         !print *,choice,pka_events_id(ievent),ievent,num_events
         IF(choice==pka_events_id(ievent)) foundi=ievent
         ievent=ievent+1
        END DO
        !PRINT *,'here',choice,ievent,pka_events_id(ievent),foundi,num_events
        
        IF(foundi.LE.num_events) THEN
         IF((pka_events_mass(foundi)==config_parent_nums(ichannel)).AND. &
(TRIM(ADJUSTL(pka_events_ele(foundi)))==TRIM(ADJUSTL(config_parent_eles(ichannel))))) THEN 
          found=.true.
         END IF
        ELSE
         ! here we have selected an atom not previously assigned
         ! to avoid conflict we must assume that it can be given the correct parent type
         ! (it has already been selected at random, so this is OK)
         ! provided we have not reached the maximum for the parent atom type
         ! first see if it (the parent) is one of the input types (must be)
         !PRINT *,'here',number_pka_files
         ifiles=1
         DO 
           IF(ifiles.GT.number_pka_files) cycle choiceloop ! shouldn't happen
           !PRINT *,parent_num(i),config_parent_nums(ichannel)
           !PRINT *,TRIM(ADJUSTL(parent_ele(i))),TRIM(ADJUSTL(config_parent_eles(ichannel)))
           IF((parent_num(ifiles)==config_parent_nums(ichannel)).AND. &
(TRIM(ADJUSTL(parent_ele(ifiles)))==TRIM(ADJUSTL(config_parent_eles(ichannel))))) EXIT
           ifiles=ifiles+1
         END DO
         !PRINT *,'here2',ifiles,initial_atom_types_count(ifiles),REAL(natoms,DBL)*pka_ratios(ifiles)
         IF(initial_atom_types_count(ifiles).GE.NINT(REAL(natoms,DBL)*pka_ratios(ifiles))) THEN
          ! cannot have any more target atoms of this type
          ! need to cycle and find one already defined
          cycle choiceloop
         END IF
         initial_atom_types_count(ifiles)=initial_atom_types_count(ifiles)+1
         found=.true.
         !ievent=num_events+1
         !num_events=num_events+1
         !16/4/19 - don't need this for on the fly assignment
         !IF(num_events.GT.max_events) THEN
         ! PRINT *,'maximum number of PKA events exceeded'
         ! PRINT *,'consider reducing timestep'
         ! STOP
         !END IF
         !pka_events_step(ievent)=0 ! so that we know to reset velocities, etc
         ! don't need to assign atom type officially as is about to change
        END IF ! ievent<num_events
         
       END DO choiceloop!found
       
       !15/4/2019
       ! WILL always be new event - check for repeat atoms now only confirms that
       ! it is possible for a previously defined atom to undergo this event
       ievent=num_events+1
       num_events=num_events+1
       
       !PKA found
       pka_events_id(ievent)=choice
       pka_events_mass(ievent)=config_daughter_nums(ichannel)
       pka_events_ele(ievent)=config_daughter_eles(ichannel)
       
       ! now define velocity/energy
       !IF(pka_events_step(ievent).NE.isteps) THEN
        pka_events_vec(ievent,:)=0._DBL ! reset for new timestep (only once)
        pka_events_ke(ievent)=0._DBL
        pka_min_sep(ievent,1)=REAL(box_nunits,DBL)*latt
        pka_min_sep(ievent,2:)=0._DBL
       !END IF
       pka_events_step(ievent)=isteps
       
       ! 15/4/2019 - select from cumulative
       CALL random_number(ireal)
       ienergy=1
       DO
         ! do we need to worry about going beyond config_global_num_pka_recoil_points?no?
         IF(SUM(config_pka_vectors(ichannel,1:ienergy))/total_channelpkas.GE.ireal) THEN
          ! we have found correct energy bin
          exit
         END IF
         ienergy=ienergy+1
       END DO
       !PRINT *,ienergy,ireal,SUM(config_pka_vectors(ichannel,1:ienergy)),total_channelpkas
       
       
       CALL random_number(ireal)
       IF(ienergy==1) THEN
       
       !always use half of lower now
       !IF(do_user_output_energy_grid) THEN
       ! ke=ireal*totalglobal_pka_recoil_energies_master(ienergy)
       !ELSE
        !between half of lowest bin and lowest bin
        ke=(1._DBL+ireal)*config_global_pka_recoil_energies(ienergy)/2._DBL
       !END IF       


       ELSE
        ke=(1._DBL-ireal)*config_global_pka_recoil_energies(ienergy-1)+&
            ireal*config_global_pka_recoil_energies(ienergy)
       END IF
       !PRINT *,ke,ienergy,config_global_pka_recoil_energies(ienergy)
       pka_events_ke(ievent)=pka_events_ke(ievent)+ke
       !to MeV to J
       ke=ke/j_to_mev
       CALL random_number(jrealvector(1))
       CALL random_number(jrealvector(2))
       CALL random_number(jrealvector(3))
       ! 2/8/2019 - don't need to convert to velocity here - do it later
       ! and instead store direction vector.
       ! however this should be between -1 and 1.
       jrealvector=jrealvector*2._DBL-1._DBL
       pka_events_vec(ievent,:)=pka_events_vec(ievent,:)+&
                  (jrealvector(:))
	 
	 
	 
	 
       
       
                  
                       
                  
     END DO !ipkas, num_PKAs_step
    

   
    !END DO ! ienergy
   END DO !ichannel
   
 IF(do_output_configs) THEN 
  ! this will not work unless the atoms are identified in advance
  ! think we need to wait until the end - even for the initial config.
  print *,'output atomic configuration',isteps
  write(tempstr,*) isteps
  outstr=TRIM(ADJUSTL(tempstr))//'_config.cfg'
  !CALL write_cfg(outstr,x,natoms,lx,io_quit_config,atom_element,atom_mass,xe)
  outstr=TRIM(ADJUSTL(tempstr))//'_config.lammmps'
  !CALL write_lammps(outstr,x,natoms,lx,io_quit_config,atom_element,atom_mass,xe)   
 END IF  
 
  !total_seps=REAL(num_events,DBL)*REAL(num_events-1,DBL)/2._DBL
  !i=0
 IF(((num_events-num_events_previous)>0).AND.(num_events>1)) THEN
  ! add to separation totals - cycle over all events in outer to save time
  DO jevent=1,num_events
   SELECT CASE(box_type)
   CASE(1)
    call define_atom_position_bcc(pka_events_id(jevent),jrealvector)
   CASE(2)
    call define_atom_position_fcc(pka_events_id(jevent),jrealvector)
   CASE(3)
    call define_atom_position_hcp(pka_events_id(jevent),jrealvector)    
   CASE DEFAULT
    call define_atom_position_bcc(pka_events_id(jevent),jrealvector)
   END SELECT
   ! locate pka in output grid
   ! want "greater than values", so 
   jenergy=1
   DO
    IF(jenergy>totalglobal_num_pka_recoil_points_master) exit
    IF(pka_events_ke(jevent)<=totalglobal_pka_recoil_energies_master(jenergy)) THEN
     
     exit
    END IF
    jenergy=jenergy+1
   END DO
   ! pka is in the greater than jenergy-1 bin
   IF(jevent>num_events_previous) THEN
    ! add to events binning
    event_dist(jenergy)=event_dist(jenergy)+1
   END IF
    
   DO ievent=MAX(jevent,num_events_previous+1),num_events
    IF(ievent==jevent) cycle
   !PRINT *,jevent,ievent,pka_events_id(jevent),pka_events_id(ievent)
   SELECT CASE(box_type)
   CASE(1)
    call define_atom_position_bcc(pka_events_id(ievent),irealvector)
   CASE(2)
    call define_atom_position_fcc(pka_events_id(ievent),irealvector)
   CASE(3)
    call define_atom_position_hcp(pka_events_id(ievent),irealvector)    
   CASE DEFAULT
    call define_atom_position_bcc(pka_events_id(ievent),irealvector)
   END SELECT   

    
    current_sep=0._DBL
    DO k=1,3
     IF(ABS(irealvector(k)-jrealvector(k))> &
       REAL(box_nunits,DBL)*latt/2._DBL) THEN
       current_sep=current_sep+&
       (ABS(irealvector(k)-jrealvector(k))- &
        REAL(box_nunits,DBL)*latt/2._DBL)**2  
     ELSE
       current_sep=current_sep+&
       (irealvector(k)-jrealvector(k))**2
     END IF  
    END DO
    !total_sep_dist(1)=total_sep_dist(1)+sqrt(current_sep)
    !IF(sqrt(current_sep)<pka_min_sep(ievent,1)) THEN
    !   pka_min_sep(ievent,1)=sqrt(current_sep)
    !   pka_min_sep(ievent,2)=REAL(pka_events_id(jevent),DBL)
    !END IF
    !IF(sqrt(current_sep)<pka_min_sep(jevent,1)) THEN
    !   pka_min_sep(jevent,1)=sqrt(current_sep)
    !   pka_min_sep(jevent,2)=REAL(pka_events_id(ievent),DBL)
    !END IF   
    
   ienergy=1
   DO
    IF(ienergy>totalglobal_num_pka_recoil_points_master) exit
    IF(pka_events_ke(ievent)<=totalglobal_pka_recoil_energies_master(ienergy)) THEN
     
     exit
    END IF
    ienergy=ienergy+1
   END DO    
   !write(701,*) jevent,ievent,ienergy,jenergy,pka_min_sep(jevent,:),pka_min_sep(ievent,:)
   DO ii=1,MIN(ienergy,jenergy)
    ! can contribute to all bins up to lowest PKA energy
    total_sep_dist(ii)=total_sep_dist(ii)+sqrt(current_sep)
    total_seps(ii)=total_seps(ii)+1._DBL
    IF((sqrt(current_sep)<pka_min_sep(ievent,2*ii-1)).OR.&
      (pka_min_sep(ievent,2*ii-1)==0._DBL)) THEN
       pka_min_sep(ievent,2*ii-1)=sqrt(current_sep)
       pka_min_sep(ievent,2*ii)=REAL(pka_events_id(jevent),DBL)
    END IF
    IF((sqrt(current_sep)<pka_min_sep(jevent,2*ii-1)).OR.&
      (pka_min_sep(jevent,2*ii-1)==0._DBL)) THEN
       pka_min_sep(jevent,2*ii-1)=sqrt(current_sep)
       pka_min_sep(jevent,2*ii)=REAL(pka_events_id(ievent),DBL)
    END IF 
   END DO   
   !write(701,*) jevent,ievent,ienergy,jenergy,pka_min_sep(jevent,:),pka_min_sep(ievent,:)
    
     
   END DO ! ievent
  END DO ! jevent   
  outstr=''
  DO i=1,totalglobal_num_pka_recoil_points_master
    ! available events in range needs to be above 1 for there to be a pair
    IF(SUM(event_dist(i+1:totalglobal_num_pka_recoil_points_master+1))>1._DBL) THEN
    !start from 1, but actually second entry pair
    write(tempstr,'(F20.5)') &
        SUM(pka_min_sep(1:num_events,2*i+1))/&
           SUM(event_dist(i+1:totalglobal_num_pka_recoil_points_master+1))
        ! event_dist for pairs greater than i is at i+1
    ELSE 
      write(tempstr,'(A20)') 'none'
    END IF
    outstr=TRIM(outstr)//TRIM(tempstr)
  END DO
  
  write(pka_analysis,'(x,I19,F20.5,I20,2(F20.5),2I20,A)') isteps,&
             total_sep_dist(1)/total_seps(1),NINT(total_seps(1)), &
             SUM(pka_min_sep(1:num_events,1))/REAL(num_events,DBL), &
             MINVAL(pka_min_sep(1:num_events,1)),&
             pka_events_id(MINLOC(pka_min_sep(1:num_events,1))), &
             NINT(pka_min_sep(MINLOC(pka_min_sep(1:num_events,1)),2)), &
               TRIM(outstr)
  
  !PRINT *,pka_min_sep(1:num_events,1)
  !PRINT *,pka_events_id(1:num_events)
  !write(700,*) num_events,event_dist
 ELSE
   jenergy=1
   jevent=1
   DO
    IF(jenergy>totalglobal_num_pka_recoil_points_master) exit
    IF(pka_events_ke(jevent)<=totalglobal_pka_recoil_energies_master(jenergy)) THEN
     
     exit
    END IF
    jenergy=jenergy+1
   END DO  
   event_dist(jenergy)=event_dist(jenergy)+1
 END IF ! num_events
 IF((num_events-num_events_previous)>0) THEN
  !write(pka_events,*) '#  ',isteps
  DO ievent=num_events_previous+1,num_events
   !IF(pka_events_step(ievent)==isteps) THEN
    !i=i+1
    
   SELECT CASE(box_type)
   CASE(1)
    call define_atom_position_bcc(pka_events_id(ievent),jrealvector)
   CASE(2)
    call define_atom_position_fcc(pka_events_id(ievent),jrealvector)
   CASE(3)
    call define_atom_position_hcp(pka_events_id(ievent),jrealvector)    
   CASE DEFAULT
    call define_atom_position_bcc(pka_events_id(ievent),jrealvector)
   END SELECT     
    

    
    
    
    
    
    ! to normalize direction vector to one
    !  and convert mass to kg via mass/(1000*AN) to define velocity vector
    ireal=sqrt(DOT_PRODUCT(pka_events_vec(ievent,1:3),pka_events_vec(ievent,1:3)))* &
         sqrt(2._DBL*ke*1000._DBL/(pka_events_mass(ievent)*avogadro))
 
    
    write(pka_events,'(x,I19,3F20.5,3ES20.5,2x,A15,3x,I15,5x,I20,2ES20.9,I10,3ES20.5)') &
          pka_events_id(ievent),jrealvector,pka_events_vec(ievent,:)/ireal, &
          pka_events_ele(ievent),pka_events_mass(ievent), &
          pka_events_step(ievent),pka_events_ke(ievent),cumtime,ievent,pka_events_vec(ievent,:)
	  
! 2/8/2019 - call to sdtrimsp to generate a cascade of the right size.
! need to filter on threshold energies to avoid unnecessary calls - done in bca routine
   bca_flag=.false.
   if (do_bca) THEN
   SELECT CASE(bca_code)
    CASE(1) ! SDTrimSP
     CALL sdtrim_run(pka_events_ke(ievent)*1e6_DBL, &
           pka_events_vec(ievent,:),&
           ievent,jrealvector,pka_events_ele(ievent),isteps,bca_flag) 
    CASE DEFAULT ! SDTrimSP
     CALL sdtrim_run(pka_events_ke(ievent)*1e6_DBL, &
           pka_events_vec(ievent,:),&
           ievent,jrealvector,pka_events_ele(ievent),isteps,bca_flag) 
   END SELECT 
   end if !do_bca
   if(bca_flag) num_above_threshold=num_above_threshold+1
    
   !END IF
  END DO  !IF(i>0) 
  !write(pka_events,*)
 END IF
  
   isteps=isteps+1
   cumtime=cumtime+timestep   
   PRINT *,'   Total time: ',cumtime, ' s'
   WRITE(log_unit,*) '   Total time: ',cumtime, ' s'
   
   if(do_bca) THEN
   
   WRITE(*,*) '   ',num_events-num_events_previous,' PKAs in timestep (', &
    num_above_threshold,' above threshold and sent to BCA)'    
   WRITE(log_unit,*) '   ',num_events-num_events_previous,' PKAs in timestep (', &
    num_above_threshold,' above threshold and sent to BCA)'   
   
   ELSE
    WRITE(*,*) '   ',num_events-num_events_previous,' PKAs in timestep'    
   WRITE(log_unit,*) '   ',num_events-num_events_previous,' PKAs in timestep'   
     
   END IF
   num_events_previous=num_events
   
  END DO !isteps
  
  
  
CLOSE(pka_events)
CLOSE(pka_analysis)  

if(do_bca) call bca_close()



END SUBROUTINE create_configs


SUBROUTINE define_atom(atmass,atele,numbertype)

 use configs
 implicit none
 !integer, intent(in) :: nn
 REAL (KIND=DBL) :: ireal
 LOGICAL :: assigned
 INTEGER :: ipp
  INTEGER, intent(out) :: atmass,numbertype
  CHARACTER (LEN=2), intent(out) :: atele
 
 ! new approach to define as needed
 
 !loop until assigned
 assigned=.false.
 DO WHILE (.not.assigned)
  CALL random_number(ireal) !progression in pka list
  ipp=1
  DO
   IF(ipp.GT.number_pka_files) EXIT !impossible
   IF(sum(pka_ratios(1:ipp)).GT.ireal) EXIT
   ipp=ipp+1
  END DO
  IF(initial_atom_types_count(ipp).GE.NINT(REAL(natoms,DBL)*pka_ratios(ipp))) THEN
    ! already at max
    CYCLE
  END IF
  atmass=parent_num(ipp)
  atele=parent_ele(ipp)
  numbertype=ipp
  initial_atom_types_count(ipp)=initial_atom_types_count(ipp)+1  
  assigned=.true.
 END DO
 
 
END SUBROUTINE define_atom


SUBROUTINE initial_config
  use configs
  implicit none
  INTEGER :: jcount,kcount,choice,inn,numbertype
  CHARACTER (LEN=500) :: outstr
  REAL (KIND=DBL) :: ireal,xp(3),pi
      
      pi=atan2(0._DBL,-1._DBL)

       
       
       io_quit_config=0
       
  IF(do_outputs) WRITE(*,*) "creating atomic lattice structure"   
  IF(box_type==1) THEN
   ! we don't do this anymore - create on the fly
   !CALL create_bcc(x,box_nunits,natoms,latt,lx)
   natoms=2*box_nunits**3
    lx(1)=REAL(box_nunits,DBL)*latt
    lx(2)=REAL(box_nunits,DBL)*latt
    lx(3)=REAL(box_nunits,DBL)*latt   
  ELSEIF (box_type==2) THEN
   !CALL create_fcc(x,natoms)
   natoms=4*box_nunits**3
    lx(1)=REAL(box_nunits,DBL)*latt
    lx(2)=REAL(box_nunits,DBL)*latt
    lx(3)=REAL(box_nunits,DBL)*latt   
  ELSEIF (box_type==3) THEN
   !CALL create_fcc(x,natoms)
   natoms=2*box_nunits**3   
    lx(1)=REAL(box_nunits,DBL)*latt
    lx(2)=REAL(box_nunits,DBL)*latt*cos(pi/6.0_DBL)
    lx(3)=REAL(box_nunits,DBL)*latt*sqrt(3.0_DBL)   
  ELSE
    print *,'not a valid lattice type, skipping cfg creation'
    io_quit_config=1
  END IF
  ALLOCATE(initial_atom_types_count(number_pka_files))
  initial_atom_types_count=0


    call random_seed(size=getsize)
    allocate (seed_array(getsize))
    seed_array=time() 
    call random_seed(put=seed_array) 
    deallocate (seed_array)
  
  IF(do_output_configs) THEN
   print *,'output initial atomic configuration'
  outstr='initial.cfg'   
   CALL write_cfg_start(outstr,natoms,lx,io_quit_config,.false.)
   outstr='initial.lammps'
   CALL write_lammps_start(outstr,natoms,lx,io_quit_config)
   
    WRITE(lammps_file,'(I5,A15)') number_pka_files,' atom types'
    WRITE(lammps_file,'(A10)') 'Masses'
    WRITE(lammps_file,*)
    DO i=1,number_pka_files
     WRITE(lammps_file,'(I5,I10)') i,parent_num(i)
    END DO
    
    WRITE(lammps_file,'(A)') 'Atoms'
    WRITE(lammps_file,*)   
   
   
  
  
  ALLOCATE(atom_mass(natoms),atom_element(natoms))
  
  atom_mass=0
  ! allocate elements by input composition ratios
  
  kcount=0
  print *,'Starting configuration:'
  initial_atom_types_count=0
  DO inn=1,natoms
   !parent_ele(1,numner_pka_files)
   !parent_num
   !pka_ratios
   
   CALL define_atom(atom_mass(inn),atom_element(inn),numbertype)
   SELECT CASE(box_type)
   CASE(1)
    call define_atom_position_bcc(inn,xp)
   CASE(2)
    call define_atom_position_fcc(inn,xp)
   CASE(3)
    call define_atom_position_hcp(inn,xp)    
   CASE DEFAULT
    call define_atom_position_bcc(inn,xp)
   END SELECT    
   ! and write to output files if required
   IF(do_output_configs) THEN
         WRITE(cfg_file,'(I3)') atom_mass(inn)
	 WRITE(cfg_file,'(A2)') atom_element(inn)
         WRITE(cfg_file,'(1x,6(1x,e16.8))') &
	    xp(1)/lx(1),xp(2)/lx(2),xp(3)/lx(3)	 
    WRITE(lammps_file,'(I10,I5,3e16.8,A15)') inn,numbertype,xp(1),xp(2),xp(3)
   
   END IF
  END DO
  DO i=1,number_pka_files
   print *,initial_atom_types_count(i),'  ',parent_ele(i),parent_num(i),' atoms'
  END DO
  
  
 

  !CALL write_cfg(outstr,x,natoms,lx,io_quit_config,atom_element,atom_mass,xe)
  
  !CALL write_lammps(outstr,x,natoms,lx,io_quit_config,atom_element,atom_mass,xe)
  
  CLOSE(cfg_file)
  CLOSE(lammps_file,IOSTAT=io_quit) 
  
 END IF
 !STOP
  
  END SUBROUTINE initial_config
  
  
  
  SUBROUTINE add_to_config_pkas(input_pka,recoil_points,energy_vector)
   use globals
   IMPLICIT NONE
   INTEGER, INTENT (in) :: recoil_points
   REAL (KIND=DBL), INTENT(in) :: input_pka(recoil_points)
   REAL (KIND=DBL), INTENT(in) :: energy_vector(recoil_points)
   REAL (KIND=DBL)  :: sum_temp
   LOGICAL :: include
   REAL (KIND=DBL), ALLOCATABLE :: tdam_pka_temp(:),pka_temp(:)
   INTEGER :: ii,jj !ii stores isotope position, jj the element position
   

   include=.true.
   
   IF((config_do_exclude_light_pkas).AND.&
    ((TRIM(ADJUSTL(daughter_ele))=='He').OR.(TRIM(ADJUSTL(daughter_ele))=='H'))) THEN
    WRITE(log_unit,*) 'omitting gas particle from config pkas'
    include=.false.
   END IF
   IF(TRIM(ADJUSTL(daughter_ele))=='unknown') THEN
    WRITE(log_unit,*) 'omitting unknown daughter from config pkas'
    include=.false. 
   END IF   
   
IF((sum(input_pka(1:recoil_points)).NE.0).AND.(include)) THEN
    ALLOCATE(pka_temp(MAX(config_global_num_pka_recoil_points,recoil_points)))
    pka_temp(1:recoil_points)=input_pka(:)*pka_ratios(filenum)
    
    
  config_num_pka_vectors=config_num_pka_vectors+1
    
  IF(config_num_pka_vectors.GT.config_max_pka_vectors) THEN
    IF(do_outputs) THEN
      WRITE(log_unit,*) 'maximum number of pka channels to consider in atomic configuration creation exceeded'
      WRITE(log_unit,*) 'all subsequent pka vectors will only be included if they exceed (in total number of pkas)'
      WRITE(log_unit,*) 'one of those already selected (which they will replace)'
    END IF
    ! select minimum
    jj=1
    
    ! 10/12/2019 - should only sum against energies above threshold.
    
    DO ii=2,config_max_pka_vectors
     IF(sum(config_pka_vectors(ii,config_threshold_group:config_global_num_pka_recoil_points)).LT. &
        sum(config_pka_vectors(jj,config_threshold_group:config_global_num_pka_recoil_points))) THEN
        jj=ii
     END IF
    END DO   
    IF(sum(config_pka_vectors(jj,config_threshold_group:config_global_num_pka_recoil_points)).LT. &
      sum(pka_temp(config_threshold_group:recoil_points))) THEN
       include=.true.
    ELSE
      include=.false.
    END IF
    config_num_pka_vectors=config_max_pka_vectors
  ELSE
   jj=config_num_pka_vectors
   include=.true.
  END IF !greater than max

   
 IF(include) THEN
   IF(do_outputs) WRITE(log_unit,*) parent_ele(filenum),parent_num(filenum),daughter_ele,daughter_num,'adding to config pkas'


    CALL collapse_fluxes(pka_temp, &
                  recoil_points,&
                  config_global_num_pka_recoil_points,&
                  config_global_pka_recoil_energies,&
                  energy_vector)  

                  
   
    ! 24/4/2018 - need to compute tdam_pka_temp, pka_temp, etc seperately
    config_pka_vectors(jj,1:config_global_num_pka_recoil_points)=pka_temp(1:config_global_num_pka_recoil_points)
    config_daughter_eles(jj)=daughter_ele
    config_parent_eles(jj)=parent_ele(filenum)
    config_daughter_nums(jj)=daughter_num
    config_parent_nums(jj)=parent_num(filenum)
    config_pka_strings(jj)=pka_element       

 END IF !include

  
  
  DEALLOCATE(pka_temp)
  
END IF ! greater than zero
  
  
  END SUBROUTINE add_to_config_pkas 
