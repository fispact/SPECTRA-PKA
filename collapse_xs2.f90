!begin rubric
!******************************************************************************                                                                         *
!* Original Author: Mark Gilbert                                                       *
!* mark.gilbert@ukaea.uk,                                                     *
!******************************************************************************
!end rubric

 SUBROUTINE collapse_xs2(xxs,ninput_groups,noutput_groups,output_energies,input_energies)
 use globals
 IMPLICIT NONE
 
 INTEGER, INTENT (in) :: ninput_groups,noutput_groups
 REAL(KIND=DBL), INTENT (inout) :: xxs(MAX(noutput_groups+1,ninput_groups+1))
 REAL(KIND=DBL), INTENT (in) :: output_energies(noutput_groups+1),input_energies(ninput_groups+1)
 INTEGER :: ii,jj,nid(noutput_groups+1),kk
 REAL(KIND=DBL) :: xxs_rounded(MAX(noutput_groups+1,ninput_groups+1))
 integer :: nd1,ni
 logical :: split
 REAL (KIND=DBL) :: frac,frac2,ebottom,due,dud
 
!30/6/2016 modified from CEM's grpconvert routine in FISPACT-II - itself based on routines written by N.P.Taylor.
! converts a cross section vector in one group structure (input_energies) into another (output_energies)
! assumes that the first ninput_energies elements of the input xxs contains the input vector
! and that the input xxs is big enough to hold the cross section in either the input or output structures
! the first noutput_energies of the returned xxs vector will be the converted cross sections.
 
     xxs_rounded=0.0
     nd1=noutput_groups+1
     ni=1
 
     jj=0
     split=.false.
     outer: do kk=ni,nd1
        ii = kk - 1
 
        inner: do
           if(.not.split) then
              jj = jj + 1
              if(jj>ninput_groups) exit outer
              
           end if
           if(input_energies(jj+1)>real(output_energies(ii+1))) exit inner
 
           if(kk/=1) then
              ebottom = input_energies(jj)
              if(split) ebottom = output_energies(ii)           
              due=input_energies(jj+1)-ebottom
              dud=output_energies(ii+1)-output_energies(ii)
              frac2 = due / dud
              xxs_rounded(ii) = xxs_rounded(ii) + xxs(jj) * frac2
             ! print *,frac,jj,ii
           end if
 
           split = .false.
           if(input_energies(jj+1) == real(output_energies(ii+1))) cycle outer
        end do inner
 
        ebottom = input_energies(jj)
        if(split) ebottom = output_energies(ii)
        


        
        ! defines overlap of out bin with inner
        due = output_energies(ii+1)-ebottom
        
        
        ! dealing with probabilities
        ! so denominator is always width of out bin
        ! i.e. to take a probability average
        !dud = input_energies(jj+1) -ebottom   
        dud=output_energies(ii+1)-output_energies(ii)
        
        frac2 = due / dud
 
        if(kk/=1) then
           xxs_rounded(ii) = xxs_rounded(ii) + xxs(jj) * frac2
           !print *,frac2,jj,ii
        end if
 
        split = .true.
    end do outer
 
 

 xxs=0._DBL
 xxs(1:noutput_groups)=xxs_rounded(1:noutput_groups)

 
 END SUBROUTINE collapse_xs2