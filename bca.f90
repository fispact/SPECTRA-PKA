MODULE bca
  use globals
  use file_units
  
  IMPLICIT NONE
  
  REAL(KIND=DBL), parameter :: min_pka_energy=20.0 ! eV
  REAL(KIND=DBL), parameter :: xstart=100._DBL
  INTEGER, ALLOCATABLE :: bca_cells(:,:) 
       ! nlc,4 - cell, number events,timestep & ievent(last occupied)
  INTEGER :: nlcx,nlcy,nlcz,nlc_max,nlc
  INTEGER :: previous_itime
  



contains

FUNCTION transformx(x,y,z,x0,y0,z0)
 IMPLICIT NONE
 
 !direction of SDTRIM PKA is (1,0,0)
 REAL(KIND=DBL),PARAMETER :: x1=1.0_DBL
 REAL(KIND=DBL) :: transformx
 REAL(KIND=DBL) :: x,y,z ! position of recoil event
 REAL(KIND=DBL) :: x0,y0,z0 !initial vector of PKA 

transformx=(1._DBL-((x1**2)*(y0**2+z0**2)/(1._DBL+x1*x0)))*x-x1*y0*y-x1*z0*z

END FUNCTION

FUNCTION transformy(x,y,z,x0,y0,z0)
 IMPLICIT NONE
 
 !direction of SDTRIM PKA is (1,0,0)
 REAL(KIND=DBL),PARAMETER :: x1=1.0_DBL
 REAL(KIND=DBL) :: transformy
 REAL(KIND=DBL) :: x,y,z ! position of recoil event
 REAL(KIND=DBL) :: x0,y0,z0 !initial vector of PKA 


transformy=x1*y0*x+(1-(x1**2)*(y0**2)/(1+x1*x0))*y-((x1**2)*y0*z0/(1+x1*x0))*z

END FUNCTION

FUNCTION transformz(x,y,z,x0,y0,z0)
 IMPLICIT NONE
 
 !direction of SDTRIM PKA is (1,0,0)
 REAL(KIND=DBL),PARAMETER :: x1=1.0_DBL
 REAL(KIND=DBL) :: transformz
 REAL(KIND=DBL) :: x,y,z ! position of recoil event
 REAL(KIND=DBL) :: x0,y0,z0 !initial vector of PKA 

transformz=x1*z0*x-((x1**2)*z0*y0/(1+x1*x0))*y+(1-(x1**2)*z0**2/(1+x1*x0))*z

END FUNCTION


END module bca


SUBROUTINE bca_setup()
 use bca
 use configs
 IMPLICIT NONE
 
 
 ! define link cell numbers
 nlcx=int(lx(1)/bca_cell_size)
 nlcy=int(lx(2)/bca_cell_size)
 nlcz=int(lx(3)/bca_cell_size)
 
 ! allocate array to hold populations events in each cell
 !ALLOCATE(bca_cells(0:nlcz*nlcy*nlcx-1))
 nlc_max=nlcz*nlcy*nlcx
 !bca_cells=0
 nlc=0
 previous_itime=0
 
OPEN(bca_output,file="config_bca.dat",STATUS='REPLACE',IOSTAT=io_open)
IF(io_open==0) THEN
  WRITE(bca_output,'(8A20)') "PKA id","Time-step","x","y","z","PKA=1,Recoil=2","recoil-level","recoil-energy-eV"

ELSE
 print *,"error opening bca output file"
 STOP
END IF


OPEN(bca_analysis,file="analysis_bca.dat",STATUS='REPLACE',IOSTAT=io_open)
IF(io_open==0) THEN
  !WRITE(bca_analysis,'(7A20)') "PKA id","Time-step","x","y","z","PKA=1,Recoil=2","recoil-level"

ELSE
 print *,"error opening bca analysis file"
 STOP
END IF

OPEN(bca_overlap,file="bca_overlap.dat",STATUS='REPLACE',IOSTAT=io_open)
IF(io_open==0) THEN
  !WRITE(bca_analysis,'(7A20)') "PKA id","Time-step","x","y","z","PKA=1,Recoil=2","recoil-level"

ELSE
 print *,"error opening bca overlap file"
 STOP
END IF
 
END SUBROUTINE bca_setup

SUBROUTINE bca_close()
use bca
 
CLOSE(bca_output)
CLOSE(bca_analysis)
CLOSE(bca_overlap)

 
END SUBROUTINE bca_close

SUBROUTINE sdtrim_run(ke,direction_vec,ievent,pos_vec,pkaelement,itime,done_bca)
  use bca
  implicit none
  
  REAL (KIND=DBL), intent(in) :: ke,direction_vec(3),pos_vec(3)
  INTEGER, intent(in) :: ievent,itime
  CHARACTER (LEN=2), intent(in) :: pkaelement
  CHARACTER (LEN=200) :: sdstr,sdstr2
  CHARACTER (LEN=500) :: commandstring
  LOGICAL, intent(out) :: done_bca
 
 ! skip unless above min energy=threshold energy
 
 IF(ke.GE.assumed_ed) THEN
  CALL triinpwrite(ke,pkaelement)
 
  write(commandstring,*) TRIM(ADJUSTL(sdtrim_path))//' >tri.cmd.out'
  CALL system(TRIM(ADJUSTL(commandstring)))
  
  
  write(sdstr2,*) ievent
  

CALL read_bca(ievent,direction_vec,pos_vec,itime)

  IF(do_store_bca_output) THEN
   sdstr="mv trajec_stop_p.dat trajec_stop_p_"//TRIM(ADJUSTL(sdstr2))//".dat"
   CALL system(sdstr)
   sdstr="mv trajec_stop_r.dat trajec_stop_r_"//TRIM(ADJUSTL(sdstr2))//".dat"
   CALL system(sdstr) 
   
   sdstr="mv partic_stop_p.dat partic_stop_p_"//TRIM(ADJUSTL(sdstr2))//".dat"
   CALL system(sdstr)
   sdstr="mv partic_stop_r.dat partic_stop_r_"//TRIM(ADJUSTL(sdstr2))//".dat"
   CALL system(sdstr) 
  
  
   

   
  
  END IF 
  
   done_bca=.true. 
  
  !STOP
 ELSE
  done_bca=.false.
 END IF ! above threshold

  
END SUBROUTINE sdtrim_run


SUBROUTINE read_bca(ievent,direction_vec,pos_vec,itime)
 use bca
 use configs
 IMPLICIT NONE
 
 CHARACTER (LEN=100) :: filestr,dummystr,stub_string(2)
 REAL (KIND=DBL), intent(in) :: direction_vec(3),pos_vec(3)
 INTEGER, intent(in) :: ievent,itime
 REAL (KIND=DBL):: dummyarray(10),zpos,ypos,xpos,histene,mag
 INTEGER :: ix,iy,iz,ip,num_ips,bcai,bcaj
 LOGICAL :: fileend,foundc
 INTEGER, ALLOCATABLE :: bca_cells_temp1(:,:),bca_cells_temp2(:,:)
 




 stub_string(1)="trajec_stop_p"
 stub_string(2)="trajec_stop_r"

num_ips=0 
DO i=1,2
 write(dummystr,*) ievent
 !write(filestr,*) TRIM(ADJUSTL(stub_string(i)))//"_"//TRIM(ADJUSTL(dummystr))//".dat"
 write(filestr,*) TRIM(ADJUSTL(stub_string(i)))//".dat"
 OPEN(bca_result,file=TRIM(ADJUSTL(filestr)),STATUS='OLD',IOSTAT=io_open)
 IF(io_open==0) THEN ! file exists and has opened OK
  fileend=.false.
fileread:  DO WHILE(.not.fileend)
 
   
   READ(bca_result,*) dummystr
   IF(TRIM(ADJUSTL(dummystr))=="ende") THEN
    fileend=.true.
    cycle
   END IF
   BACKSPACE(bca_result)
   READ(bca_result,*,IOSTAT=io_read) dummyarray(:)
   !PRINT *,io_read,NINT(dummyarray(8))
  IF((io_read==0).AND.(NINT(dummyarray(8))==1)) THEN
   ! define according to coordinate transform
   
   mag=sqrt(DOT_PRODUCT(direction_vec,direction_vec))
   !mag=sqrt(direction_vec(1)**2+direction_vec(2)**2+direction_vec(3)**2)
   xpos=transformx(dummyarray(2)-xstart,dummyarray(3),dummyarray(4),&
        direction_vec(1)/mag,direction_vec(2)/mag,direction_vec(3)/mag)
   ypos=transformy(dummyarray(2)-xstart,dummyarray(3),dummyarray(4),&
        direction_vec(1)/mag,direction_vec(2)/mag,direction_vec(3)/mag)
   zpos=transformz(dummyarray(2)-xstart,dummyarray(3),dummyarray(4),&
        direction_vec(1)/mag,direction_vec(2)/mag,direction_vec(3)/mag)
   histene=dummyarray(6)
   
 !17/9/2019 - as standard we will skip if below the displacement threshold 
 ! this time it will be the full displacement thresh (not the PKA min)
 ![ update 3/1/2020 - now we don't intiate bca unless above threshold]
 ! atom will not be displaced (and thus won't contribute to damage) if 
 ! below E_d
 IF(histene.LT.assumed_ed) cycle fileread
   !shift according to PKA position
   xpos=xpos+pos_vec(1)
   ypos=ypos+pos_vec(2)
   zpos=zpos+pos_vec(3)
   
   
   ! and shift through (assumed) PBCs as required
   IF(do_bca_pbc) THEN
    IF(xpos>lx(1)) xpos=xpos-lx(1)
    IF(ypos>lx(2)) ypos=ypos-lx(2)
    IF(zpos>lx(3)) zpos=zpos-lx(3)
   
    IF(xpos<=0._DBL) xpos=xpos+lx(1)
    IF(ypos<=0._DBL) ypos=ypos+lx(2)
    IF(zpos<=0._DBL) zpos=zpos+lx(3)  
   ELSE 
    !skip as appropriate
    IF(xpos>lx(1)) cycle fileread
    IF(ypos>lx(2)) cycle fileread
    IF(zpos>lx(3)) cycle fileread
   
    IF(xpos<=0._DBL) cycle fileread
    IF(ypos<=0._DBL) cycle fileread
    IF(zpos<=0._DBL) cycle fileread    
   END IF
   
   
   ! define cell
   
   ix=INT(xpos*REAL(nlcx,DBL)/lx(1))+1
   iy=INT(ypos*REAL(nlcy,DBL)/lx(2))+1
   iz=INT(zpos*REAL(nlcz,DBL)/lx(3))+1
   ip=ix-1+nlcx*(iy-1+nlcy*(iz-1))
   !write(200,*) ip,xpos,ypos,zpos,ix,iy,iz,i,dummyarray(10),ievent,pka_events_step(ievent)
   write(bca_output,'(2I20,3F20.8,2I20,F20.8)') ievent,pka_events_step(ievent),&
       xpos,ypos,zpos,i,NINT(dummyarray(10)) ,histene
   foundc=.false.
   bcai=1
   DO 
    IF(bcai.GT.num_ips) exit
    IF(bca_cells_temp1(bcai,1)==ip) THEN
     foundc=.true.
     exit
    END IF
    bcai=bcai+1
   END DO
   IF(ip.LT.0) THEN
       print *,xpos,ypos,zpos,ix,iy,iz,ievent,itime,dummyarray
       STOP
   END IF
   IF(.not.foundc) THEN
    IF(num_ips==0) THEN
     num_ips=num_ips+1
     ALLOCATE(bca_cells_temp1(num_ips,2))
     bca_cells_temp1(num_ips,1)=ip
     bca_cells_temp1(num_ips,2)=1
    ELSE
     ALLOCATE(bca_cells_temp2(num_ips,2))
     bca_cells_temp2(1:num_ips,:)=bca_cells_temp1(1:num_ips,:)
     DEALLOCATE(bca_cells_temp1)
     num_ips=num_ips+1
     ALLOCATE(bca_cells_temp1(num_ips,2))
     bca_cells_temp1(1:num_ips-1,:)=bca_cells_temp2(1:num_ips-1,:)
     DEALLOCATE(bca_cells_temp2)
     bca_cells_temp1(num_ips,1)=ip
     bca_cells_temp1(num_ips,2)=1
    END IF
   ELSE
    bca_cells_temp1(bcai,2)=bca_cells_temp1(bcai,2)+1
   END IF
    
   !bca_cells(ip)=bca_cells(ip)+1
 
  END IF !dummyarray(8)
  END DO fileread !io_read filend
  CLOSE(bca_result)
 ELSE
   PRINT *,'error opening BCA result file',TRIM(ADJUSTL(filestr))
 END IF !io_open
 
 
END DO ! i
write(bca_output,*)
write(bca_output,*)

IF(num_ips.GT.0) THEN
IF(nlc==0) THEN
 nlc=num_ips
 ALLOCATE(bca_cells(nlc,4))
 bca_cells(1:nlc,1:2)=bca_cells_temp1(1:nlc,1:2)
 bca_cells(:,3)=itime
 bca_cells(:,4)=ievent
 previous_itime=itime
 
ELSE


! output previous if required
IF(previous_itime.LT.itime) THEN ! new time
  call write_bca_analysis
  previous_itime=itime
END IF


!now add to global bca array
ALLOCATE(bca_cells_temp2(nlc,4))
bca_cells_temp2(1:nlc,:)=bca_cells(1:nlc,:)
DEALLOCATE(bca_cells)
ALLOCATE(bca_cells(nlc+num_ips,4)) !max size
bca_cells(1:nlc,:)=bca_cells_temp2(1:nlc,:)
DEALLOCATE(bca_cells_temp2)
DO bcai=1,num_ips
 foundc=.false.
 bcaj=1
 DO
  IF(bcaj.GT.nlc) exit
  IF(bca_cells(bcaj,1)==bca_cells_temp1(bcai,1)) THEN
   foundc=.true.
   
  
   exit
  END IF
  bcaj=bcaj+1
 END DO !bcaj
 IF(foundc) THEN

   ! this our first test for overlap
   ! this will only check if the overlap has occurred in  previous cascade
   ! will not define time. (unless we track time).
   IF(ievent.NE.bca_cells(bcaj,4)) THEN
    WRITE(bca_overlap,*) 'cascade overlap at timestep ',itime, &
    'overlapping with cascade at timestep ',bca_cells(bcaj,3)
    IF(overlap_stop) THEN
      print *,'existing as requested at first cascade overlap'
      STOP
    END IF
   ELSE
    WRITE(bca_overlap,*) 'two cascades in timestep ',itime,' have overlapped'
    WRITE(bca_overlap,*) 'consider reducing the timestep in simulations'
   END IF

  bca_cells(bcaj,2)=bca_cells(bcaj,2)+bca_cells_temp1(bcai,2)
  bca_cells(bcaj,3)=itime
  bca_cells(bcaj,4)=ievent
 ELSE
  !new ip
  nlc=nlc+1
  bca_cells(nlc,1:2)=bca_cells_temp1(bcai,1:2)
  bca_cells(nlc,3)=itime
  bca_cells(nlc,4)=ievent
 
 END IF
END DO !bcai
   
   
END IF !nlc=0   


DEALLOCATE(bca_cells_temp1)
END IF !num_ips>0

 
 
!STOP
 
 
 
END SUBROUTINE read_bca



SUBROUTINE write_bca_analysis()
 use bca
 use configs
 IMPLICIT NONE
 INTEGER :: ix,iy,iz
 
 
! and output - but wait until we reach the end of the timstep

  IF(previous_itime.GT.0) THEN ! not fool proof, may still output blank lines at start
   write(bca_analysis,*)
   write(bca_analysis,*)
  END IF
  ! we are writing the previous state at start of new time.
  write(bca_analysis,*) "# timestep: ",previous_itime
  

DO i=1,nlc
  IF(bca_cells(i,2).GT.0) THEN
   ix=mod(bca_cells(i,1),nlcx)+1
   iy=mod((bca_cells(i,1)-ix+1)/nlcx,nlcy)+1
   iz=(((bca_cells(i,1)-ix+1)/nlcx)-iy+1)/nlcy+1
   write(bca_analysis,*) bca_cells(i,1),bca_cells(i,2),ix,iy,iz,&
     REAL((ix-1),DBL)*lx(1)/REAL(nlcx,DBL), &
     REAL((iy-1),DBL)*lx(2)/REAL(nlcy,DBL), &
     REAL((iz-1),DBL)*lx(3)/REAL(nlcz,DBL),i
     
     
  END IF
END DO 
 
 
END SUBROUTINE write_bca_analysis




SUBROUTINE triinpwrite(ke,pkaelement)
  use bca
  
  implicit none
  real (KIND=DBL), intent(in) :: ke
  character (LEN=200) :: tstr,tstr2
  CHARACTER (LEN=2), intent(in) :: pkaelement
  
  
  
  
  OPEN(triinp,file="tri.inp",STATUS='REPLACE')

! 25/11/2019 - put PKA species first
! so that it is easier to track if parent is PKA


write(triinp,*) "run "  
write(triinp,*) "   &TRI_INP  " 

! number of ncp is total types in material plus one extra for projectile
! start with input composition (later could change)
write(triinp,*) "     ncp = ", number_pka_files+1
tstr=""
tstr=TRIM(ADJUSTL(tstr))//'"'//TRIM(ADJUSTL(pkaelement))//'"'
DO i=1,number_pka_files
 tstr=TRIM(ADJUSTL(tstr))//',"'//TRIM(ADJUSTL(parent_ele(i)))//'"'
END DO
write(triinp,*) '     symbol ='//TRIM(ADJUSTL(tstr))
     
write(triinp,*) "     flc      =   10.0  "
write(triinp,*) "     nh       =   1 "
write(triinp,*) "     idout    =   1 "
write(triinp,*) "     nr_pproj =   1 "     
write(triinp,*) "     idrel = 1 "
write(triinp,*) "     isbv  = 1 "
write(triinp,*) "     ipot  = 1 "



tstr=""
tstr=TRIM(ADJUSTL(tstr))//'1.00'
DO i=1,number_pka_files
 tstr=TRIM(ADJUSTL(tstr))//',0.00'
END DO
write(triinp,*) "     qubeam = "//TRIM(ADJUSTL(tstr))

tstr=""
tstr=TRIM(ADJUSTL(tstr))//'1.00'
DO i=1,number_pka_files
 tstr=TRIM(ADJUSTL(tstr))//',1.00'
END DO
write(triinp,*) "     qumax= "//TRIM(ADJUSTL(tstr))

write(triinp,*) "     case_e0=0 "

tstr=""
write(tstr2,'(F10.2)') ke
tstr=TRIM(ADJUSTL(tstr))//TRIM(ADJUSTL(tstr2))
DO i=1,number_pka_files
 tstr=TRIM(ADJUSTL(tstr))//',0.00'
END DO
write(triinp,*) "     e0= "//TRIM(ADJUSTL(tstr))


tstr=""
write(tstr2,'(F8.1)') xstart
tstr=TRIM(ADJUSTL(tstr))//TRIM(ADJUSTL(tstr2))
DO i=1,number_pka_files
 tstr=TRIM(ADJUSTL(tstr))//',0.00'
END DO

write(triinp,*) "     x0= "//TRIM(ADJUSTL(tstr))

write(triinp,*) "     case_alpha=0 "


tstr=""
tstr=TRIM(ADJUSTL(tstr))//'0.0'
DO i=1,number_pka_files
 tstr=TRIM(ADJUSTL(tstr))//',0.00'
END DO
write(triinp,*) "     alpha0= "//TRIM(ADJUSTL(tstr))

write(triinp,*) "     ttarget = 5000 "
write(triinp,*) "     nqx     =  500 "


tstr=""
tstr=TRIM(ADJUSTL(tstr))//'0.0'
DO i=1,number_pka_files
 write(tstr2,'(F10.2)') pka_ratios(i)
 tstr=TRIM(ADJUSTL(tstr))//','//TRIM(ADJUSTL(tstr2))
END DO
write(triinp,*) "     qu= "//TRIM(ADJUSTL(tstr))

tstr=""
tstr=TRIM(ADJUSTL(tstr))//'20.00'
DO i=1,number_pka_files
 tstr=TRIM(ADJUSTL(tstr))//',20.00'
END DO
write(triinp,*) "     e_cutoff= "//TRIM(ADJUSTL(tstr))

tstr=""
write(tstr2,'(F10.2)') assumed_ed
tstr=TRIM(ADJUSTL(tstr))//TRIM(ADJUSTL(tstr2))
DO i=1,number_pka_files
 tstr=TRIM(ADJUSTL(tstr))//','//TRIM(ADJUSTL(tstr2))
END DO
write(triinp,*) "     e_displ= "//TRIM(ADJUSTL(tstr))


write(triinp,*) "     ltraj_p   = .true.  "  
write(triinp,*) "     ltraj_r   = .true. "   
write(triinp,*) "     numb_hist = 10000 "   
write(triinp,*) "     ioutput_hist = 10000, 10, 0, 10000, 10, 0,"
write(triinp,*) "     lparticle_p=.true."
write(triinp,*) "     lparticle_r=.true."
write(triinp,*) "     ioutput_part=10000,100,0,10000,100,0"
write(triinp,*) "     tableinp = '/work/SDTrimSP/tables' "
write(triinp,*) "/ "


CLOSE(triinp)





end SUBROUTINE triinpwrite
  
  
