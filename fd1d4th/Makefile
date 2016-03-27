NETCDFF :=/opt/netcdf4/4.2.2-ifort
FC := ifort
FFLAGS := -c -w95 -warn all -zero -module ./ #  -I$(NETCDFF)/include 
FFLAGSMP := -O3 -mp1 -unroll0 -ip-no-inlining $(FFLAGS)
OBJ = fd1dstg.o io_mod.o 
fd : $(OBJ)
	@$(FC)  $^  -o $@  -L$(NETCDFF)/lib -lnetcdff
fd1dstg.o :  fd1dstg.F90 io_mod.o 
	@$(FC) $(FFLAGSMP) $< -o $@
io_mod.o : io_mod.F90
	$(FC) $(FFLAGS) -DGETARG $< -o $@
.PHONY : clean
clean :
	rm -f fd $(OBJ) 
	rm -f *.mod
