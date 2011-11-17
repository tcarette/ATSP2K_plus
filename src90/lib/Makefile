# makefile for Applications

TARGETS= install #debug
SUBDIRS= libmod libang90 libcom90 libdvdson90 librad90 mpi

$(TARGETS) : $(SUBDIRS)
install:
	@for i in $(SUBDIRS) ; \
	do \
                cd $$i ; \
                echo "....Entering: $$i" ; \
                $(MAKE) $@ ; \
                cd .. ; \
                echo "....Leaving: $$i" ; echo; echo;\
	done

clean:
	@for i in $(SUBDIRS); \
	do \
	cd $$i ; \
	echo "Current directory: $$i" ; \
		$(MAKE) $@ ; \
	        cd .. ; \
	done

