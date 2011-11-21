
################################################################################
# Makefile for all applications. 
# .For each subdirectory the default target "install" is recursively built 
# ."make clean" will descend into each directory and remove *.o files.
################################################################################


TARGETS= install 
#SUBDIRS= hfs
#SUBDIRS= iso_save
#SUBDIRS= iso
#SUBDIRS= mchf_sort
#SUBDIRS= mchf_iso 
#SUBDIRS= mpi_mchf mpi_nonh
#SUBDIRS= mpi_nonh
#SUBDIRS= mchf mpi_mchf
#SUBDIRS= biomat
#SUBDIRS= lsreduce mpi_mchf
#SUBDIRS= lsreduce
#SUBDIRS= lsreduce hfs mpi_mchf mpi_nonh nonh mchf
#SUBDIRS= mchf
SUBDIRS= nonh
#SUBDIRS= bpci
#SUBDIRS= bp_eiv bp_ang bp_mat mpi_bp_ang mpi_bp_eiv mpi_bp_mat
#SUBDIRS= bp_eiv bp_ang bp_mat bpci hfs iso lsreduce mchf nonh trans biotr \
#         mpi_bp_ang mpi_bp_eiv mpi_bp_mat mpi_mchf \
#         mpi_nonh mchf_iso


$(TARGETS): $(SUBDIRS)
install:
	@for i in $(SUBDIRS) ; \
	do \
                cd $$i ; \
                echo "....Entering: $$i" ; \
                $(MAKE) $@ ; \
                cd .. ; \
                echo "....Leaving: $$i" ; echo; echo;\
	done

backup:
	@for i in $(SUBDIRS); \
	do \
	cd $$i ; \
	echo "Current directory: $$i" ; \
		$(MAKE) $@ ; \
		cd .. ; \
	done

clean:
	@for i in $(SUBDIRS); \
	do \
	cd $$i ; \
	echo "Current directory: $$i" ; \
	        $(MAKE) $@ ; \
	        cd .. ; \
	done

