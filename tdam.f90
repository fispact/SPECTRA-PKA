!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk
!* copyright, 2015 (first version), 2018 (first git repository), UKAEA                                                     *
!******************************************************************************
!end rubric

     SUBROUTINE calc_tdam(num_points,Epkas,tdams,arec,zrec,alat,zlat)
      use globals
      IMPLICIT NONE
      INTEGER, INTENT(in) :: num_points
      INTEGER, INTENT(in) :: arec,alat,zrec,zlat
      REAL :: ar,zr,al,zl
      REAL(KIND=DBL), intent(in) :: Epkas(num_points)
      REAL(KIND=DBL), intent(out) :: tdams(num_points)
      REAL(KIND=DBL), parameter :: abohr=52.91772108_DBL !in pm
      REAL(KIND=DBL), parameter :: ec2=1.4399764_DBL !in eVnm
      REAL(KIND=DBL) :: aa,ee,gg,kk,pi,denom,rel,el
      REAL(KIND=DBL),parameter :: twothd=2._DBL/3._DBL
      REAL(KIND=DBL),parameter :: threeq=3._DBL/4._DBL
      REAL(KIND=DBL),parameter :: sixth=1._DBL/6._DBL
      REAL(KIND=DBL),parameter :: onethd=1._DBL/3._DBL
      REAL(KIND=DBL),parameter :: onep5=3._DBL/2._DBL
      REAL(KIND=DBL),parameter :: c1=30.724_DBL
      REAL(KIND=DBL),parameter :: c2=0.07952410617_DBL
      REAL(KIND=DBL),parameter :: c3=3.4008_DBL
      REAL(KIND=DBL),parameter :: c4=.40244_DBL      



     zr=REAL(zrec,DBL)
     zl=REAL(zlat,DBL)
     ar=REAL(arec,DBL)
     al=REAL(alat,DBL)

 pi=atan2(0d0,-1d0)
IF(tdam_method==1) THEN
      
!NRT paper, Nucl Eng. Des. 1975     

aa=((9._DBL*pi**2/128._DBL)**onethd)* &
      abohr/sqrt(zr**twothd+ &
      zl**twothd)
kk=0.1337_DBL*zr**sixth*&
      (zr/ar)**(1.0_DBL/2.0_DBL)
  IF(zrec==0) THEN
  tdams(1:num_points)=0._DBL
  ELSE
DO i=1,num_points
  !convert from MeV to eV in ee formula (and then abohr to nm, for combined 1e3 value), but not in final tdam.
  ee=(al*Epkas(i)*1e3_DBL/(ar+al)) &
      *aa/(zr*zl*ec2)
  gg=3.4008_DBL*ee**sixth+0.40244_DBL*ee**threeq+ee
  
  ! no need for keV conversion here
  tdams(i)=Epkas(i)/(1.0_DBL+kk*gg)
END DO
  END IF



ELSEIF(tdam_method==2) THEN

!njoy coding and also Robinson, JNM, 1994
 ! 12/7/2016
 !c2=0.0793=(32/3*pi)*sqrt(m_e)
 !      =(3.395)*sqrt(9.1093897e-28 g *(A_N=6.0221367e23))
 !      =0.07952410617
 !      = let's use this instead
 
  ! where does the c1=30.724 constant come from?
  ! 
  !(9pi^2/128)^(1/3)*abohr/ec2=32.535..
  !1/32.535=0.030736
  ! but then need convert to bohr to nm 
  ! so becomes 1e3*0.030736=30.736 - very close and probably differences in constants.
  ! let's use real value 
 
el=(ec2/(abohr*1e-3_DBL*(9._DBL*pi**2/128._DBL)**onethd))*zr*zl* &
   sqrt(zr**twothd+ &
        zl**twothd&
       )* &
   (ar+al)/al
rel=1._DBL/el
denom=(zr**twothd+zl**twothd)** &
       threeq*ar**onep5*sqrt(al)
kk=c2*zr**twothd*sqrt(zl)* &
                     (ar+al)**onep5/ &
                     denom
                     
                     
                   
                     
  IF(zrec==0) THEN
  tdams(1:num_points)=0._DBL
  ELSE
DO i=1,num_points
!convert Epkas(i) in MeV to eV (rel is in eV from ec2) to get correct ratio
ee=Epkas(i)*1e6_DBL*rel   
gg=(c3*ee**sixth+c4*ee**threeq+ee)
tdams(i)=Epkas(i)/(1.0_DBL+kk*gg)    

END DO
  END IF
END IF


 
  

     
     END SUBROUTINE calc_tdam
