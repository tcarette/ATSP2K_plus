EXED = ${ATSP}/bin
TAR = `date '+%y_%m_%d:%H'`_brci.tar   # tar -file : year_month_date_hour
# executable targets
MSN = $(EXED)/nonh_mpi

# librar
LIBDIR = ${ATSP}/lib

RAD_COM = -lrad  -lcom -lrad
SPARSE = -ldvdson -lnet
# Libraries to be searched
LIBUSE = -lcom -ldvdson -lang  -lrad -lcom -lmpiu

#FC_FLAGS = -g -dp -e0
#linker flags

#objects
OBJECT_SN  = mpi_nonh_${FC_MALLOC}.o mpi_angmom.o genint.o mpi_orbitals.o\
	isearch.o ltriang.o  savels_u_du.o ffgg.o analy_blk.o\
	cfgblk.o 

install : $(MSN)
$(MSN): $(OBJECT_SN)
	$(FC_MPI) -o $(MSN) $(FC_MPILD)  -L$(LIBDIR) \
		$(OBJECT_SN) $(LIBUSE)  
archive:
	tar -cvf $(TAR) *.f Makefile

clean:
	rm -rf *.o core

SOURCE_SNONH = mpi_snonh.f angmomg_u.f genint.f mpi_orbitals.f\
        isearch.f ltriang.f  savels_u_du.f ffgg.f  analy_blk.f\
        cfgblk.f


SNONH.f: $(SOURCE_SNONH)
	cat $(SOURCE_SNONH) > SNONH.f 
	
.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC_MPI) -c $(FC_MPIFLAGS) $<  -I . -o $@
.f.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) $< -o $@

