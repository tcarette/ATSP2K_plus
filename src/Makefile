# makefile for all

TARGETS= install 
#SUBDIRS= lib appl util
#SUBDIRS= appl util
SUBDIRS= lib appl util
#SUBDIRS= util
#SUBDIRS= appl
MKDIR = ${ATSP} ${ATSP}/bin ${ATSP}/lib

$(TARGETS):
	@for D in $(MKDIR) ; \
        do \
           test -d $$D || echo ... creating a new directory $$D... ; \
           test -d $$D || mkdir $$D ; \
        done

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

