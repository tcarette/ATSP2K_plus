
################################################################################
# Makefile for all applications. 
# .For each subdirectory the default target "install" is recursively built 
# ."make clean" will descend into each directory and remove *.o files.
################################################################################


TARGETS= install 

SUBDIRS =  biotr_ang biotr_tr 

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

