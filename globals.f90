!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk,                                                     *
!******************************************************************************
!end rubric

MODULE globals
 USE accuracy
 
 
 integer, parameter ::flux_unit=21,results_unit=22,kerma_unit=23,pka_unit=24,gas_unit=25, &
                      sigma_unit=26,unitin=20,index_summary=27,user_ebinsunit=28
 CHARACTER (LEN=200) :: flux_filename,results_filename,kerma_filename, &
                        gas_filename,sigma_filename,results_stub
 CHARACTER (LEN=100) :: flux_title
 REAL(KIND=DBL) :: rtemp
 
 INTEGER :: i,j,k,l,m,file_index
 INTEGER :: io_open,io_quit,io_read
 
 
 ! variables in flux input file
 INTEGER :: itype,isigma,igroup,ipka,number_flux_groups,ksail,number_flux_ebins
 INTEGER :: flux_norm_type
 REAL (KIND=DBL) :: acnm,time,total_fluence,norm_fluence,flux_energy_rescale_value, &
                    flux_rescale_value
 REAL (KIND=DBL), ALLOCATABLE :: flux_energies(:),fluxes(:),flux_covariances(:,:),&
              flux_ebin_widths(:),fluxes_norm(:)
 
 
 

 
 
 !pka variables
 INTEGER :: num_pka_elements
 INTEGER :: num_pka_points ! number of PKA bins
 INTEGER :: num_pka_incident_energies ! number of incident energy boundaries
 INTEGER :: num_pka_recoil_points ! number of PKA bin boundaries
 REAL (KIND=DBL), ALLOCATABLE :: pka_incident_energies(:),pka_recoil_energies(:),&
                                 recoil_kermas(:,:),pka(:,:),epka(:)
 CHARACTER (LEN=30) :: pka_element
 LOGICAL :: at_end
 REAL(KIND=DBL) :: sum_temp,sum_temp2,pka_ave
 INTEGER :: pka_filetype
 !20/6/2013 - variables for pka summing
 INTEGER :: mtd,num_pka_recoil_points_master
 LOGICAL :: do_mtd_sums
 REAL(KIND=DBL),ALLOCATABLE :: pka_sums(:,:),pka_recoil_energies_master(:)
 
 !26/6/2013 flags to switch of summing if total already exists
 LOGICAL :: alpha_sum_flag,proton_sum_flag
 
 !17/7/2013 - flag to switch basic energy balance estimate for (n,g) recoils
 LOGICAL :: do_ngamma_estimate
 INTEGER :: ng_num_pka_recoil_points,ng_num_pka_incident_energies,ng_num_pka_points
 REAL(KIND=DBL), ALLOCATABLE :: ng_pka(:,:),ng_pka_recoil_energies(:)
 REAL (KIND=DBL), ALLOCATABLE :: ng_pka_incident_energies(:),&
                                 ng_recoil_kermas(:,:)
 
 !8/10/2013
 ! daughter definitions
 CHARACTER (LEN=2),ALLOCATABLE :: parent_ele(:)
 CHARACTER (LEN=10) :: daughter_ele,number_string,number_string2
 INTEGER,ALLOCATABLE :: parent_num(:)
 INTEGER :: daughter_num
 
 
 !9/10/2013 ! multiple file handling and summing
 REAL (KIND=DBL),ALLOCATABLE ::pka_ratios(:)
 CHARACTER (LEN=200),ALLOCATABLE :: pka_filename(:)
 INTEGER :: number_pka_files,filenum
 LOGICAL :: do_global_sums
 INTEGER :: max_global_recoils
 REAL(KIND=DBL),ALLOCATABLE :: global_pka_sums(:,:),global_pka_recoil_energies_master(:)
 CHARACTER (LEN=10), allocatable :: global_daughter_eles(:)
 INTEGER, allocatable :: global_daughter_nums(:)
 INTEGER :: global_num_pka_recoil_points_master,number_global_recoils
 
 !24/4/2018 - separate total results to allow separate use of user grid
 REAL(KIND=DBL),ALLOCATABLE :: totalglobal_pka_recoil_energies_master(:)
 INTEGER :: totalglobal_num_pka_recoil_points_master


 !10/10/2013 elemental sums
 ! INTEGER, PARAMETER :: max_global_elements=40
 ! 17/8/2016 use same max array size as recoil globals
 REAL(KIND=DBL),ALLOCATABLE :: global_pka_sums_element(:,:)
 CHARACTER (LEN=10), ALLOCATABLE :: global_elements(:)
 INTEGER :: number_global_recoil_elements
 
 !26/2/2014  - alternate ngamma estimate
 REAL (KIND=DBL), PARAMETER :: avogadro=6.022141930e23_DBL, &
                               clight=299792458._DBL, &
                               j_to_mev=1._DBL/(1.602176565e-13_DBL), &
                               neutron_mass=1.008664923_DBL, &
                               proton_mass=1.007825032_DBL
 REAL(KIND=DBL), ALLOCATABLE :: ngamma_daughter_mass(:),ngamma_parent_mass(:)
 
 
 !20/5/2014 - total sum flag to exclude He/H
 LOGICAL :: do_exclude_light_from_total
 REAL(KIND=DBL),ALLOCATABLE :: total_pka_sum(:)
 
 
 !8/4/2015 - handle input data for each pka file in column mode.
 INTEGER :: num_columns
 
 !4/6/2015 master element list is now a global
 ! use to order elemental sum outputs
 INTEGER, PARAMETER :: num_master_elements=109
 CHARACTER (LEN=2), PARAMETER :: master_elements(0:num_master_elements)= &
        (/'n ','H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
          'Na','Mg','Al','Si','P ','S ','Cl','Ar', &
          'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',&
          'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
          'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',&
          'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',&
          'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg', &
         'Bh','Hs','Mt'/)
         
 
 
 !5/10/2015 - way to handle new njoy output with reduced energy info.+4/3/2016
 LOGICAL :: energies_once_perfile,first_read
 REAL(KIND=DBL),ALLOCATABLE :: pka_incident_energies_filemaster(:)
 INTEGER :: num_pka_incident_energies_filemaster
 !24/4/2018 - revised after user grd changes (these are held in the master vectors
 REAL(KIND=DBL),ALLOCATABLE :: pka_recoil_energies_filemaster(:)
 INTEGER :: num_pka_recoil_points_filemaster 
 
 !4/3/2016 exclude unknowns from total - needed for proton matrix problems.
 LOGICAL :: do_exclude_unknown_from_total
 
 
 !21/3/2016 calculate damage energy for each PKA energy & estimate displacements
 LOGICAL :: do_tdam
 REAL (KIND=DBL), ALLOCATABLE :: tdam_energies(:),ng_tdam_energies(:)
 INTEGER :: daughter_z,parent_z,tdam_method
 !4/4/2016 - calculate displacements
 REAL (KIND=DBL) :: displacements,assumed_ed,total_disp
 REAL (KIND=DBL), ALLOCATABLE :: global_disp_sums(:),global_disp_sums_element(:)  ! for species/element sums
 REAL (KIND=DBL), ALLOCATABLE :: mtd_disp_sums(:)   ! for alpha/proton sums
 
 character (LEN=1) :: incident_particle
 
 
 REAL (KIND=DBL), ALLOCATABLE :: global_pka_sums_element_tdam(:,:), &
                                 total_pka_sum_tdam(:), &
                                 global_tdam_energies_master(:), &
                                 global_pka_sums_tdam(:,:),&
                                 tdam_energies_master(:)
 REAL (KIND=DBL), ALLOCATABLE :: pka_sums_tdam(:,:)       
 
 !16/10/2017 - handle exception caused by first file being empty and preventing array allocation
 LOGICAL :: first_non_empty
 
 
 !11/3/2018
 INTEGER :: deallocerr
 
 !23/4/2018
 ! variables to read in user-specified output group structure
 LOGICAL :: do_user_output_energy_grid
 CHARACTER (LEN=500) :: user_energybin_file
 INTEGER :: user_grid_option
 
 !24/4/2018 - extended text outputs flag
 LOGICAL :: do_outputs,doing_ng
 
 !9/5/2018
 LOGICAL :: do_timed_configs
 INTEGER :: box_nunits,box_type,nsteps
 REAL (KIND=DBL) :: timestep,latt
 CHARACTER (LEN=500) :: config_namestub
 
 end module globals
