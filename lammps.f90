! started August 2025

MODULE lammps
  use globals
  use file_units
  
  IMPLICIT NONE
  
  !REAL(KIND=DBL), parameter :: min_pka_energy=20.0 ! eV
  !REAL(KIND=DBL), parameter :: xstart=100._DBL
  !INTEGER, ALLOCATABLE :: bca_cells(:,:) 
       ! nlc,4 - cell, number events,timestep & ievent(last occupied)
       ! stores events per cell of box
  !INTEGER :: nlcx,nlcy,nlcz,nlc_max,nlc
  !INTEGER :: previous_itime
  



contains




END module lammps


SUBROUTINE lammps_setup()
 use lammps

 IMPLICIT NONE
 
 
 ! define link cell numbers
 !nlcx=int(lx(1)/bca_cell_size)
 !nlcy=int(lx(2)/bca_cell_size)
 !nlcz=int(lx(3)/bca_cell_size)
 
 ! allocate array to hold populations events in each cell
 !ALLOCATE(bca_cells(0:nlcz*nlcy*nlcx-1))
 !nlc_max=nlcz*nlcy*nlcx
 !bca_cells=0
 !nlc=0
 !previous_itime=0
 

 
END SUBROUTINE lammps_setup



SUBROUTINE lammps_run(ke,direction_vec,ievent,pos_vec,pkaelement,itime,done_lammps)
  use lammps
  implicit none
  
  REAL (KIND=DBL), intent(in) :: ke,direction_vec(3),pos_vec(3)
  INTEGER, intent(in) :: ievent,itime
  CHARACTER (LEN=2), intent(in) :: pkaelement
  CHARACTER (LEN=200) :: sdstr,sdstr2
  CHARACTER (LEN=500) :: commandstring
  LOGICAL, intent(out) :: done_lammps
 
 ! skip unless above min energy=threshold energy
 
 IF(ke.GE.assumed_ed) THEN
  CALL lammpswrite(ke,pkaelement)
 
  !write(commandstring,*) TRIM(ADJUSTL(sdtrim_path))//' >tri.cmd.out'
  !CALL system(TRIM(ADJUSTL(commandstring)))
  
  
  write(sdstr2,*) ievent
  

!CALL read_md(ievent,direction_vec,pos_vec,itime)

  !IF(do_store_lammps_output) THEN
  ! sdstr="mv trajec_stop_p.dat trajec_stop_p_"//TRIM(ADJUSTL(sdstr2))//".dat"
  ! CALL system(sdstr)
  ! sdstr="mv trajec_stop_r.dat trajec_stop_r_"//TRIM(ADJUSTL(sdstr2))//".dat"
  ! CALL system(sdstr) 
   
  ! sdstr="mv partic_stop_p.dat partic_stop_p_"//TRIM(ADJUSTL(sdstr2))//".dat"
  ! CALL system(sdstr)
  ! sdstr="mv partic_stop_r.dat partic_stop_r_"//TRIM(ADJUSTL(sdstr2))//".dat"
  ! CALL system(sdstr) 
  
  
   

   
  
 ! END IF 
  
   done_lammps=.true. 
  
  !STOP
 ELSE
  done_lammps=.false.
 END IF ! above threshold

  
END SUBROUTINE lammps_run




SUBROUTINE lammpswrite(ke.pkaelement)
 IMPLICIT NONE
 use lammps
 INTEGER, PARAMETER :: unitlammpsin=201
 INTEGER :: io_read,io_open
 
 OPEN(unitlammps,file="lammps.in",STATUS='REPLACE',iostat=io_open)
 IF(io_open==0) THEN
 
 


write(unitlammps,*) "# Initialization"
write(unitlammps,*) "units metal"
write(unitlammps,*) "atom_style atomic"
write(unitlammps,*) "boundary p p p"

!IF(do_create_lammps_box)
 ! create_config file will have already created
 !write(unitlammps,) "read_data ",TRIM(ADJUSTL(lammps_config_file))
!ELSE
write(unitlammps,*) "# Create a simple BCC lattice"
write(unitlammps,*) "lattice bcc",latt
! this box dimension could be as specified by box_units
! but that is likely too big. 
write(unitlammps,*) "region box block 0 20 0 20 0 20"
write(unitlammps,*) "create_box 1 box"
write(unitlammps,*) "create_atoms 1 box"
!END IF
 


write(unitlammps,*) "# Define EAM potential for a metal"
write(unitlammps,*) "pair_style eam/fs"
write(unitlammps,*) "pair_coeff * * Fe_mm.eam.fs Fe"

!we may not do this if we are using an existing box
write(unitlammps,*) "# Equilibrate"
write(unitlammps,*) "minimize 1.0e-10 1.0e-10 1000 10000"
write(unitlammps,*) "fix 1 all npt temp 300.0 300.0 100.0 iso 0.0 0.0 1000.0"
!write(unitlammps,*) "dump 1 all custom 1000 dump.equil id type x y z
write(unitlammps,*) "run 1000"

write(unitlammps,*) "# PKA initiation"




write(unitlammps,*) "# Define a small region around the point of interest
write(unitlammps,*) "region my_small_region sphere 0.0 0.0 0.0 1.0 units box
write(unitlammps,*) "# Create a group of atoms within that region
write(unitlammps,*) "group near_point region my_small_region

# Then, perform the same distance calculation and min reduction, but only on this group.
# The variable and compute commands can be limited to the `near_point` group.

create_atoms 1 region mybox
variable px equal 0.0
variable py equal 0.0
variable pz equal 0.0
variable distsq atom "(x-v_px)^2 + (y-v_py)^2 + (z-v_pz)^2"
variable dist atom "sqrt(v_distsq)"
compute mindist near_point reduce min v_dist
variable closest_atom_id
variable closest_atom_x
variable closest_atom_y
variable closest_atom_z
label find_closest_atom
  # Use next to loop through all atoms
  next closest_atom_id loop 1 1 ${N}
  if "$(v_dist[v_closest_atom_id]) == c_mindist" then &
    "variable closest_atom_x equal x[v_closest_atom_id]" &
    "variable closest_atom_y equal y[v_closest_atom_id]" &
    "variable closest_atom_z equal z[v_closest_atom_id]" &
    "jump SELF break" # Exit the loop after finding the first match

print "Closest atom ID: ${closest_atom_id}"
print "Closest atom coordinates: (${closest_atom_x}, ${closest_atom_y}, ${closest_atom_z})"

# To use this in a simulation, you would run minimization or dynamics
# and then add the above logic to your script.

unfix 1
group pka id 8422 #  atom ID is near center
variable pka_energy equal 500 # 50 eV
variable pka_mass equal 56.0 # Ni mass (g/mol)
# divide by 100 to get from m/s to A/ps for metal units
# divide pka_mass (in g/mol) to g per atom by dividing by avogadro (atoms/mol)
# then divide by 1000 to convert from g to kg
# convert pka_energy in eV to J by multiplying by 1.6E-19
variable pka_velocity equal sqrt(2*${pka_energy}*1.60218e-19*6.02214e23*1000/(${pka_mass}))/100 # Convert eV to J, then calculate velocity

#velocity pka set ${pka_velocity} 0 0 units box # Example velocity in x direction

# make it resolved in x,y - want velocity in y to be half that of x
# but resolved must still equal pka_velocity (resolved vector)
# y**2=x**2+(x/2)**2, x=sqrt(4y/5)
variable pka_velocity_x equal 2.0*${pka_velocity}/sqrt(5)
variable pka_velocity_y equal ${pka_velocity_x}/2.0

velocity pka set ${pka_velocity_x} ${pka_velocity_y} 0 units box # Example velocity in x direction

print "pka_mass $(v_pka_mass) "
print "pka_velocity $(v_pka_velocity) "
print "pka_energy $(v_pka_energy) "

# Cascade simulation
fix 2 all nve
#fix 3 all temp/berendsen 300.0 300.0 100.0 region boundary # Define a boundary region
timestep 0.0001 # Small initial timestep 0.1 fs
thermo 100
dump 2 all custom 1 dump.cascade id type x y z vx vy vz
dump_modify 2 sort id
run 1000 # Example run time for cascade

compute        peratome all pe/atom
compute        peratomke all ke/atom

#more time 
timestep 0.001 # timestep 1 fs
thermo 100


run 50 # Example run time for cascade

write_dump  all custom dump.cascade_end id type x y z vx vy vz c_peratome c_peratomke modify sort id 



 ELSE
 
 END IF


END SUBROUTINE
  
