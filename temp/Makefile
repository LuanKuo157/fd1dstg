NETCDFF :=/opt/netcdf4/4.2.2-ifort
FC := ifort
OBJ = fd1dstg.o io_mod.o
fd : $(OBJ)
	@$(FC)  $^ -o $@ -L$(NETCDFF)/lib -lnetcdff
io_mod.mod : io_mod.o io_mod.F90
	@$(FC) -c  io_mod.F90 
io_mod.o : io_mod.F90
	@$(FC) -c  $^
fd1dstg.o : io_mod.mod  fd1dstg.F90
	@$(FC) -c fd1dstg.F90 
.PHONY : clean
clean :
	@rm -f fd $(OBJ) 
	rm -f *.mod
