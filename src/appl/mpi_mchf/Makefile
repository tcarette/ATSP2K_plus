EXED = ${ATSP}/bin

TAR = `date '+%y_%m_%d:%H'`_smchf_mpi.tar   # tar -file : year_month_date_hour

#executable targets
MSM = $(EXED)/mchf_mpi #+output

# libraries
#LIBDIR = /u1/georgio/atsp2K/lib #../../../lib
#LIBDIR = /u2/georgio/SPII/atsp2K/lib #../../../lib
LIBDIR = ${ATSP}/lib

RAD_COM = -lrad  -ldvdson -lcom -lrad -lmpiu  
SPARSE = -ldvdson 
LIBUSE = -lcom -ldvdson -lang  -lrad -lcom 

# objects
USPOBJ= alloc_mem.o de.o diff.o dinit.o e.o eijset.o eig_out.o \
	methd1.o nmrvs.o nodec.o orthog.o output.o quad.o scale.o\
	mpi_sintgrl.o search.o summry.o unpacki.o wavefn.o\
	mpi_mchf_${FC_MALLOC}.o mpi_data.o mpi_diag.o mpi_solve.o \
        mpi_grange.o mpi_scf.o \
	mpi_sum.o mpi_update.o mpi_updatc.o mpi_xch.o mpi_potl.o \
	mpi_dssbmv.o mpi_diag_hmx.o prprty.o mpi_spalcsts_$(FC_MALLOC).o

install: de.o $(MSM)

#de.o : de.f
#	$(FC_MPI) $(FC_FLAGS) -c de.f

$(MSM):  $(USPOBJ) 
	$(FC_MPI) -o $(MSM) ${FC_MPILD}  \
                  -L$(LIBDIR) $(USPOBJ) $(RAD_COM) ${LAPACK_DIR} ${LAPACK_LIB}


clean:
	rm -rf *.o core 
archive:
	tar -cvf $(TAR) *.f Makefile

source_smchf: $(SOURCE_SMCHF)
	cat $(SOURCE_SMCHF) > source_smchf 
	
.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC_MPI) -c $(FC_MPIFLAGS) $<  -I . -o $@
.f.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) $< -o $@

