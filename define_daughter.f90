!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk,                                                     *
!******************************************************************************
!end rubric

     SUBROUTINE define_daughter(recoil)
     use globals
      IMPLICIT NONE
      INTEGER :: ele,shift
      LOGICAL :: found,recoil
      


      
      
      SELECT CASE(incident_particle)
       case("n")
        shift=0
       case("p")
        shift=1
      
      CASE DEFAULT
       shift=0
      END SELECT
     
      
      ! define parent position in list of elements
      ele=1
      found=.false.
      DO WHILE ((.not.found).AND.(ele.LE.num_master_elements))
       IF(TRIM(ADJUSTL(parent_ele(filenum)))==TRIM(ADJUSTL(master_elements(ele))) ) THEN
        found=.true.
       ELSE
        ele=ele+1
       END IF
      END DO
      !this also defines z number
      parent_z=ele
      
      !by mtd cases
     IF(found) THEN 
      daughter_ele='unknown'
      daughter_num=0     
      IF(recoil) THEN
       SELECT CASE(mtd)
        CASE(2)
         !scattering
         daughter_ele=parent_ele(filenum)
         daughter_num=parent_num(filenum)
         daughter_z=ele
        CASE(4,51:91)
         !neutron emission
         daughter_ele=master_elements(ele+shift)
         daughter_num=parent_num(filenum)
         daughter_z=ele+shift        
        CASE(107,800:849) !(z,a)
         daughter_ele=master_elements(ele-2+shift)
         daughter_num=parent_num(filenum)-3
         daughter_z=ele-2+shift
        CASE(102) !(z,g)
         daughter_ele=master_elements(ele+shift)
         daughter_num=parent_num(filenum)+1
         daughter_z=ele+shift
        CASE(103,600:649) !(z,p)
         daughter_ele=master_elements(ele-1+shift)
         daughter_num=parent_num(filenum)
         daughter_z=ele-1+shift
        CASE(11) !z,2nd
	 daughter_ele=master_elements(ele-1+shift)
         daughter_num=parent_num(filenum)-3  
         daughter_z=ele-1+shift
        CASE(16,875:891) !z,2n
         daughter_ele=master_elements(ele+shift)
         daughter_num=parent_num(filenum)-1
         daughter_z=ele+shift
        CASE(17) !z,3n
	 daughter_ele=master_elements(ele+shift)
         daughter_num=parent_num(filenum)-2
         daughter_z=ele+shift
        CASE(22) !z,na
	 daughter_ele=master_elements(ele-2+shift)
         daughter_num=parent_num(filenum)-4 
         daughter_z=ele-2+shift
        CASE(23) !z,n3a
	 daughter_ele=master_elements(ele-6+shift)
         daughter_num=parent_num(filenum)-12 
         daughter_z=ele-6+shift
        CASE(24) !z,2na
	 daughter_ele=master_elements(ele-2+shift)
         daughter_num=parent_num(filenum)-5
         daughter_z=ele-2+shift
        CASE(25) !z,3na
	 daughter_ele=master_elements(ele-2+shift)
         daughter_num=parent_num(filenum)-6 
         daughter_z=ele-2+shift
        CASE(28) !z,np
	 daughter_ele=master_elements(ele-1+shift)
         daughter_num=parent_num(filenum)-1
         daughter_z=ele-1+shift
        CASE(29) !z,n2a
	 daughter_ele=master_elements(ele-4+shift)
         daughter_num=parent_num(filenum)-8
         daughter_z=ele-4+shift
        CASE(30) !z,2n2a
	 daughter_ele=master_elements(ele-4+shift)
         daughter_num=parent_num(filenum)-9
         daughter_z=ele-4+shift
        CASE(32) !z,nd
	 daughter_ele=master_elements(ele-1+shift)
         daughter_num=parent_num(filenum)-2
         daughter_z=ele-1+shift
        CASE(33) !z,nt
	 daughter_ele=master_elements(ele-1+shift)
         daughter_num=parent_num(filenum)-3 
         daughter_z=ele-1+shift
        CASE(34) !z,nh
	 daughter_ele=master_elements(ele-2+shift)
         daughter_num=parent_num(filenum)-3  
         daughter_z=ele-2+shift
        CASE(35) !z,nd2a
	 daughter_ele=master_elements(ele-5+shift)
         daughter_num=parent_num(filenum)-10
         daughter_z=ele-5+shift
        CASE(36) !z,nt2a
	 daughter_ele=master_elements(ele-5+shift)
         daughter_num=parent_num(filenum)-11
         daughter_z=ele-5+shift
        CASE(37) !z,4n
	 daughter_ele=master_elements(ele+shift)
         daughter_num=parent_num(filenum)-3 
         daughter_z=ele+shift
        CASE(104,650:699) !z,d
	 daughter_ele=master_elements(ele-1+shift)
         daughter_num=parent_num(filenum)-1
         daughter_z=ele-1+shift
        CASE(105,700:749) !z,t
	 daughter_ele=master_elements(ele-1+shift)
         daughter_num=parent_num(filenum)-2  
         daughter_z=ele-1+shift
        CASE(108) !z,2a
	 daughter_ele=master_elements(ele-4+shift)
         daughter_num=parent_num(filenum)-7
         daughter_z=ele-4+shift
        CASE(109) !z,3a
	 daughter_ele=master_elements(ele-6+shift)
         daughter_num=parent_num(filenum)-11   
         daughter_z=ele-6+shift
        CASE(111) !z,2p
	 daughter_ele=master_elements(ele-2+shift)
         daughter_num=parent_num(filenum)-1
         daughter_z=ele-2+shift
        CASE(112) !z,pa
	 daughter_ele=master_elements(ele-3+shift)
         daughter_num=parent_num(filenum)-4
         daughter_z=ele-3+shift
        CASE(41) !z,2np
	 daughter_ele=master_elements(ele-1+shift)
         daughter_num=parent_num(filenum)-2
         daughter_z=ele-1+shift
        CASE(42) !z,3np
	 daughter_ele=master_elements(ele-1+shift)
         daughter_num=parent_num(filenum)-3
         daughter_z=ele-1+shift
        CASE(44) !z,n2p
	 daughter_ele=master_elements(ele-2+shift)
         daughter_num=parent_num(filenum)-2
         daughter_z=ele-2+shift
        CASE(45) !z,npa
	 daughter_ele=master_elements(ele-3+shift)
         daughter_num=parent_num(filenum)-5         
         daughter_z=ele-3+shift


        CASE(106) !z,h
	 daughter_ele=master_elements(ele-2+shift)
         daughter_num=parent_num(filenum)-2
         daughter_z=ele-2+shift
        CASE(113) !z,t2a
	 daughter_ele=master_elements(ele-5+shift)
         daughter_num=parent_num(filenum)-10 
         daughter_z=ele-5+shift
        CASE(114) !z,d2a
	 daughter_ele=master_elements(ele-5+shift)
         daughter_num=parent_num(filenum)-9  
         daughter_z=ele-5+shift
        CASE(115) !z,pd
	 daughter_ele=master_elements(ele-2+shift)
         daughter_num=parent_num(filenum)-2
         daughter_z=ele-2+shift
        CASE(116) !z,pt
	 daughter_ele=master_elements(ele-2+shift)
         daughter_num=parent_num(filenum)-3 
         daughter_z=ele-2+shift
        CASE(117) !z,da
	 daughter_ele=master_elements(ele-3+shift)
         daughter_num=parent_num(filenum)-5
         daughter_z=ele-3+shift
         
     
         
        CASE DEFAULT
        
        
       END SELECT
      ELSE !no recoil
      ! 21/4/2014 case system not required for light particles.
      ! and actually doesn't work anyway for proton case.
       !SELECT CASE(mtd)
       ! CASE(5,22:25,29:30,35:36,45,107:109,112:114,117,800:849) !(z,a)
         IF(INDEX(pka_element,"alpha").NE.0) THEN
          daughter_ele='He'
          daughter_num=4
          daughter_z=2
         END IF
       ! CASE DEFAULT
       !END SELECT
      ! SELECT CASE(mtd)
      !  CASE(102) !(z,g)
        IF(INDEX(pka_element,"gamma").NE.0) THEN
         daughter_ele='gamma'
         daughter_num=0
         daughter_z=0
        END IF
      !  CASE DEFAULT
      ! END SELECT 
       !SELECT CASE(mtd)
       ! CASE(5,28,41:45,103,111:112,115:116,600:649) !protons
         IF(INDEX(pka_element,"proton").NE.0) THEN
          daughter_ele='H '
          daughter_num=1
          daughter_z=1
         END IF
       ! CASE DEFAULT
       !END SELECT         


        
        
      END IF
     ELSE
      PRINT *,'unknown element: ',parent_ele(filenum)
      daughter_ele='unknown'
      daughter_num=-1
     END IF
     
     
     END SUBROUTINE define_daughter