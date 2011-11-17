LIBDIR = $(ATSP)/lib
LIBA = $(LIBDIR)/libmpiu.a
#MPI_INC = ${MPI_INCLUDES}

LIBOBJ =  \
    mpi_work_dir_$(FC_MALLOC).o sys_chdir_$(FC_MALLOC).o  \
   sys_getwd_$(FC_MALLOC).o  sys_mkdir_$(FC_MALLOC).o mpi_allr_dp.o mpi_allr_int.o

install: $(LIBA)
$(LIBA) : $(LIBOBJ)
	@echo " Building " $(LIBA)
#	lorder $(LIBOBJ) | tsort | xargs ar cr $(LIBA)
	ar -curvs $(LIBA) $(LIBOBJ)
	ranlib $(LIBA)
	rm -f *.o

clean:
	rm -f *.o core

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC_MPI) -c $(FC_MPIFLAGS) $<  -I . -o $@
.f.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) $< -o $@

