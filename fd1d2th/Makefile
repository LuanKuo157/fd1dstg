FC = ifort 
OBJ = fd1dstg.o read_conf.o read_media.o stf_ricker.o
fd2 : $(OBJ)
	@$(FC)  -O2  $^ -o  $@ 
fd1dstg.o : fd1dstg.F90
	$(FC) -c fd1dstg.F90 
read_conf.o : read_conf.F90
	$(FC) -c  read_conf.F90
read_media.o : read_media.F90
	$(FC) -c read_media.F90
stf_ricker.o : stf_ricker.F90
	$(FC) -c stf_ricker.F90
.PHONY : clean
clean:
	@rm -f $(OBJ)
