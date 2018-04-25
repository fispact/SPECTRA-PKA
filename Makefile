###########################################
# Makefile adopted for Linux + mpich
###########################################

FF = gfortran
LD = gfortran
FFLAGS =  -ffree-line-length-none -g 
LDFLAGS = -g 
LIBS =

#FF = ifort
#LD = ifort
#FFLAGS = 
#LDFLAGS = 
#FFLAGS =  -check -debug inline_debug_info -debug-parameters all
#LDFLAGS = -check -debug inline_debug_info -debug-parameters all
#LIBS = 

PROG = spectra-pka

.SUFFIXES : .o .f90

%.o:%.f90
	$(FF) -c $(FFLAGS) $<

OBJ = accuracy.o globals.o read_input.o \
      spectra-pka.o collapse_xs2.o \
      collapse_fluxes.o define_daughter.o \
      ng_estimate.o output_sum_pkas.o \
      sum_pkas.o read_flux.o read_pka.o \
      global_sums.o tdam.o
      
$(PROG) : $(OBJ) 
	  $(LD) -o $(PROG) $(LDFLAGS) $(OBJ)
	
clean:
	rm -f $(OBJ) core *~

distclean:
	rm -f $(PROG) *.o *.mod core *~	
	
	
	
	




