!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk
!* copyright, 2015 (first version), 2018 (first git repository), UKAEA                                                     *
!******************************************************************************
!end rubric

   SUBROUTINE process_ng_estimate()
   use globals
   IMPLICIT NONE
   logical :: found_pos
   INTEGER :: lower_bin,upper_bin
   REAL(KIND=DBL) :: upper_energy_out,lower_energy_out,extra_momentum_energy, &
                     incident_mass

   
      SELECT CASE(incident_particle)
       case("n")
        incident_mass=neutron_mass
       case("p")
        incident_mass=proton_mass
      
      CASE DEFAULT
       incident_mass=neutron_mass
      END SELECT   
   
    IF(do_outputs) WRITE(log_unit,*) 'performing (n,g) evaluation'
 
    ! see pg 126 bk 5
    ALLOCATE(ng_recoil_kermas(num_pka_recoil_points,MAX(number_flux_ebins,num_pka_incident_energies)), &
    ng_pka_incident_energies(num_pka_incident_energies),ng_pka_recoil_energies(num_pka_recoil_points)) 
    
    !ng_pka_recoil_energies=pka_recoil_energies
    ng_pka_incident_energies=pka_incident_energies
    ng_pka_recoil_energies=pka_recoil_energies
    IF(do_tdam) THEN
       ALLOCATE(ng_tdam_energies(num_pka_recoil_points))
    END IF
    ng_recoil_kermas=0._DBL
    ng_num_pka_incident_energies=num_pka_incident_energies
    ng_num_pka_points=num_pka_points
    ng_num_pka_recoil_points=num_pka_recoil_points
    
        !DO i=2,ng_num_pka_recoil_points
         DO j=1,ng_num_pka_incident_energies-1
         
         
         !24/2/2014  the energy bounds associated with recoil group i are i and i-1, not i and i+1
         !upper_energy_out=ng_pka_incident_energies(j+1)-ng_pka_recoil_energies(i-1) !(i)
         !lower_energy_out=ng_pka_incident_energies(j)-ng_pka_recoil_energies(i) !(i+1)
         
         !26/2/2014 - different estimate - based on momentum conservation of recoil after gamma emission
         ! however, does not account for case when multiple gammas are emitted to create
         ! final daughter
         ! better estimate is to consider energy difference associated with mass difference
         ! between parent and daughter and (conservatively) assume all contributes to extra energy
         ! i.e. is emitted as a single gamma
         ! 6/3/2014 must also include momentum conservation associated with incident neutron
         
         !momentum is mv=m*sqrt(2E/m) - E is MeV, so divide by convertor to turn into J=kgm^2s^-2
         

         
         !conservation of momentum - momentum_before=momentum of neutron
         !velocity of new particle (ngamma_daughter_mass(filenum))=momentum_before/new_mass
         ! extra energy is then j_to_mev*0.5*new_mass*velocity^2
         !see bk 12 - pg46
         
         upper_energy_out=ng_pka_incident_energies(j+1)*incident_mass/ngamma_daughter_mass(filenum)
         
         lower_energy_out=ng_pka_incident_energies(j)*incident_mass/ngamma_daughter_mass(filenum)
         
         
         !WRITE(107,'(I5,3ES15.7)',ADVANCE='NO') j,upper_energy_out,lower_energy_out,ng_pka_incident_energies(j)
         
         
         ! and add in the energy associated with the mass difference
         ! see bk 12 pg 48
         extra_momentum_energy=(j_to_mev)*((ngamma_parent_mass(filenum)+incident_mass-ngamma_daughter_mass(filenum))**2)*clight**2/ &
              (2._DBL*avogadro*1000._DBL*ngamma_daughter_mass(filenum))
         
         
         
         upper_energy_out=upper_energy_out+extra_momentum_energy
         lower_energy_out=lower_energy_out+extra_momentum_energy

! 3/4/2014 actually use original method.         
!         !11/3/2014 - more reaslitic, however, is that the extra_momentum_energy is added
!         ! isotropically in all directions, so that the total added can be between
!         ! +extra and -extra (i.e. gamma generated in same direction as neutron direction + so 
!         ! reverse momentum subtracts from total.
!         !thus:
!         
!         upper_energy_out=upper_energy_out+extra_momentum_energy
!         lower_energy_out=lower_energy_out-extra_momentum_energy
!         ! but now lower_energy_out can be negative
!         ! so make min of zero
!         lower_energy_out=MAX(lower_energy_out,0._DBL)
         
         !WRITE(107,*) upper_energy_out,lower_energy_out,extra_momentum_energy,recoil_kermas(1,j)
         
         !this above defines the range of the recoiling energies of the heavy nucleus
         ! need to merge the recoil_kerma (in fact a cross section) using this onto normal range by adding appropriate fractions
         ! locate positions within recoil spectrum
         found_pos=.false.
         k=1
         DO WHILE ((.not.found_pos).AND.(k.LE.ng_num_pka_recoil_points))
          IF(lower_energy_out.LT.ng_pka_recoil_energies(k)) THEN           
           found_pos=.true.
          ELSE
           !24/2/2014 - k defines lower limit of bin k+1
           k=k+1
          END IF
          
         END DO
         !24/2/2014 - k defines upper limit of bin k
         lower_bin=k  !-1
         !IF(lower_bin==0) lower_bin=1
         found_pos=.false.
         k=ng_num_pka_recoil_points
         DO WHILE ((.not.found_pos).AND.(k.GE.1))
          IF(upper_energy_out.GT.ng_pka_recoil_energies(k)) THEN           
           found_pos=.true.
          ELSE
           
           k=k-1
          END IF
         END DO 
         !24/2/2014, k defines lower limit of bin k+1
         upper_bin=k+1  !k
         
         
         
         IF((upper_bin==1).OR.(lower_bin==ng_num_pka_recoil_points+1)) THEN
           !skip - less than lower bin or above upper bin
         ELSE 
          ! now we have the range of bins into which the kerma must be split

          !3/4/2014 - total width of recoil group is the same
          ! as incident group because calculations are simply multiplication(plus constant)
          ! of lower and upper incident bounds
          ! thus only split is relative sizes of recoil bins
          !divide kerma by full width ratio          
          ! and multiply by appropriate fractions
          ! whole bins or partials      
          
          !11/7/2016
          ! no Kerma is a cross section and hence a probability
          ! there is no proportion of cross section
          ! if a neutron has an energy in the group, then it has the 
          ! xs of that group
          ! so we actually want to add to each ng_recoil_kerma
          ! the kerma of the incident group multiplied by the amount
          ! of overlap with the recoil group
          ! at the end we will then divide by the recoil bin width
          ! check 10/10/2017 - the 2016 method is correct, but results in much smaller
          ! xs than the old approach (i.e. upper-lower is a smaller number in general than the recoil
          ! bin width)
          DO k=lower_bin,upper_bin
          
          ! ng_recoil_kermas(k,j)=ng_recoil_kermas(k,j)+ &
          !  (recoil_kermas(1,j)/(upper_energy_out-lower_energy_out))* &
          !  (MIN(upper_energy_out,ng_pka_recoil_energies(k))- &  !24/2/14 not k+1
          !   MAX(lower_energy_out,ng_pka_recoil_energies(k-1)))  !24/2/14 not k

          !11/7/2016 - xs (probability) times contributing width of bin
           ng_recoil_kermas(k,j)=ng_recoil_kermas(k,j)+ &
            recoil_kermas(1,j)* &
            (MIN(upper_energy_out,ng_pka_recoil_energies(k))- &  !24/2/14 not k+1
             MAX(lower_energy_out,ng_pka_recoil_energies(k-1))) !24/2/14 not k
            
             
          END DO
         END IF
         
         !11/7/2016 - now divide xs in each recoil group by bin width
         !25/4/2018 array handling for k=1
         k=1
         ng_recoil_kermas(k,j)=ng_recoil_kermas(k,j)/ &
             (ng_pka_recoil_energies(k)- &  
             ng_pka_recoil_energies(k)/2._DBL)         
         DO k=2,ng_num_pka_recoil_points
          ng_recoil_kermas(k,j)=ng_recoil_kermas(k,j)/ &
             (ng_pka_recoil_energies(k)- &  
             ng_pka_recoil_energies(k-1)) 
         END DO !k
         
         END DO !j
       ! END DO !i
 
     ALLOCATE(ng_pka(2,ng_num_pka_recoil_points))

   DO i=1,ng_num_pka_recoil_points
    !pka_fluxes(i,1:num_pka_incident_energies)
    ! collapse input_pka_energy_spectrum onto flux spectrum
    CALL collapse_xs2(ng_recoil_kermas(i,:),ng_num_pka_points,number_flux_groups,flux_energies,ng_pka_incident_energies,0._DBL)
    WHERE(ng_recoil_kermas(i,:).LT.0._DBL) ng_recoil_kermas(i,:)=0._DBL
   END DO 
   
   DO i=1,ng_num_pka_recoil_points
    ng_pka(1,i)=SUM(fluxes_norm(1:number_flux_ebins)*ng_recoil_kermas(i,1:number_flux_ebins)) 

   END DO   
   
   
   
   DEALLOCATE(ng_recoil_kermas,ng_pka_incident_energies)
   
   
   
  END SUBROUTINE process_ng_estimate
 
  SUBROUTINE output_ng_estimate()
  use globals
  IMPLICIT NONE

    !estimated (n,g) recoil matrix
    ! ng_pka(1,:)
    
   IF(do_outputs) WRITE(log_unit,*) 'output (n,g) evaluation'
    
   mtd=102
   CALL define_daughter(.true.)
    IF(do_tdam) THEN
         CALL calc_tdam(num_pka_recoil_points,ng_pka_recoil_energies,ng_tdam_energies, &
         daughter_num,daughter_z,parent_num(filenum),parent_z)
         
         !4/4/2016 calculate displacements estimate.
         displacements=0._DBL
         i=1
         ! skip NRT part until end

         
         displacements=(ng_pka(1,i))*&
             (ng_tdam_energies(i)/2._DBL)
         do i=2,ng_num_pka_recoil_points
           displacements=displacements+(ng_pka(1,i))*&
             ((ng_tdam_energies(i-1)+ng_tdam_energies(i))/2._DBL)
         end do          
         
         
    END IF
    IF(do_outputs) WRITE(log_unit,*) 'add ng estimate to globals'
     doing_ng=.true.
   IF(do_global_sums) CALL add_to_globals(ng_pka(1,:),ng_num_pka_recoil_points,ng_tdam_energies,ng_pka_recoil_energies)
   j=1
   
   IF(sum(ng_pka(j,1:ng_num_pka_recoil_points)).NE.0) THEN
      WRITE(results_unit,*) '### index ',file_index,' ##### ',TRIM(ADJUSTL(pka_filename(filenum)))
      WRITE(number_string,'(I5)') daughter_num
      WRITE(number_string2,'(I5)') parent_num(filenum)
      WRITE(results_unit,*) '#  estimated (n,g) recoil matrix [ '//TRIM(ADJUSTL(daughter_ele)) &
       //TRIM(ADJUSTL(number_string))//' ] from [ '//TRIM(ADJUSTL(parent_ele(filenum))) &
       //TRIM(ADJUSTL(number_string2)) &
       //' ]'
    WRITE(index_summary,*) file_index,' estimated (n,g) recoil matrix [ '//TRIM(ADJUSTL(daughter_ele)) &
       //TRIM(ADJUSTL(number_string))//' ] from [ '//TRIM(ADJUSTL(parent_ele(filenum))) &
       //TRIM(ADJUSTL(number_string2)) &
       //' ]'//TRIM(ADJUSTL(pka_filename(filenum)))        
   ! WRITE(results_unit,*) '#DIFF SPEC AVG. PKA RECOIL DISTRIBUTIONS'
   WRITE(results_unit,*) '#PKA RECOIL DISTRIBUTIONS'
       
       
       i=1
      IF(do_tdam) THEN
       WRITE(results_unit,'(1x,a31,a8,2x,a11,a28,2x,A22,A12)') '#RECOIL energy (MeV low & high)','PKAs','norm_sum', &
              '  T_dam (MeV low and high)','disp_energy (eV/s)','NRT_dpa/s'
              
       IF(ng_pka(j,i).NE.0) WRITE(results_unit,'(2ES16.6,4ES11.4,2ES20.4)') &
               ng_pka_recoil_energies(i)/2._DBL,ng_pka_recoil_energies(i), &
                 ng_pka(j,i),ng_pka(j,i)/SUM(ng_pka(j,1:ng_num_pka_recoil_points)), &
                 ng_tdam_energies(i)/2._DBL,ng_tdam_energies(i),(ng_pka(1,i))*&
             (ng_tdam_energies(i)/2._DBL)*1e6_DBL, &
             0.8_DBL*(pka(1,i))*&
             (ng_tdam_energies(i)/2._DBL)/(2._DBL*assumed_ed*1e-6_DBL)
       DO i=2,ng_num_pka_recoil_points
        IF(ng_pka(j,i).NE.0) WRITE(results_unit,'(2ES16.6,4ES11.4,2ES20.4)') &
                  ng_pka_recoil_energies(i-1),ng_pka_recoil_energies(i), &       
                 ng_pka(j,i),ng_pka(j,i)/SUM(ng_pka(j,1:ng_num_pka_recoil_points)), &
                 ng_tdam_energies(i-1),ng_tdam_energies(i),(ng_pka(1,i))*&
             ((ng_tdam_energies(i-1)+ng_tdam_energies(i))/2._DBL)*1e6_DBL, &
             0.8_DBL*(ng_pka(1,i))*&
             ((ng_tdam_energies(i-1)+ng_tdam_energies(i))/2._DBL)/(2._DBL*assumed_ed*1e-6_DBL)
       END DO
      ELSE
       WRITE(results_unit,'(1x,a31,a8,a11)') '#RECOIL energy (MeV low & high)','PKAs','norm_sum'
       IF(ng_pka(j,i).NE.0) WRITE(results_unit,'(2E14.4,5x,2ES11.4)') &
               ng_pka_recoil_energies(i)/2._DBL,ng_pka_recoil_energies(i), &
                 ng_pka(j,i),ng_pka(j,i)/SUM(ng_pka(j,1:ng_num_pka_recoil_points))       
       DO i=2,ng_num_pka_recoil_points
        IF(ng_pka(j,i).NE.0) WRITE(results_unit,'(2E14.4,5x,2ES11.4)') &
                  ng_pka_recoil_energies(i-1),ng_pka_recoil_energies(i), &       
                 ng_pka(j,i),ng_pka(j,i)/SUM(ng_pka(j,1:ng_num_pka_recoil_points))
       END DO      
      END IF



  
      ! compute average pka energy
      pka_ave=0.5_DBL*ng_pka_recoil_energies(1)*ng_pka(j,1)
      DO i=2,ng_num_pka_recoil_points
       pka_ave=pka_ave+0.5_DBL*(ng_pka_recoil_energies(i)+ng_pka_recoil_energies(i-1))*ng_pka(j,i)
      END DO
      pka_ave=pka_ave*1000000._DBL
      
      WRITE(results_unit,'(a1,20x,A,ES11.4,A5)') '#','AVERAGE PKA ENERGY = ',pka_ave/SUM(ng_pka(j,:)),' (eV)'
   IF (do_tdam) then
    WRITE(results_unit,'(a1,20x,A,ES11.4)') '#',&
               'displacement energy eV/s = ',displacements*1e6_DBL 
    WRITE(results_unit,'(a1,20x,A,ES11.4,A7,F4.1,a)') '#','equivalent NRT dpa/s = ',&
              0.8_DBL*displacements/(2._DBL*assumed_ed*1e-6_DBL),&
            ' (E_d=',assumed_ed,' eV)' 
   END IF      
      WRITE(results_unit,*)
      WRITE(results_unit,*)
      file_index=file_index+1
   END IF !non zero check
   IF(do_outputs) WRITE(log_unit,*) 'end of (n,g) output'
  
  END SUBROUTINE output_ng_estimate
