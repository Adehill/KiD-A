#
# Makefile for KiD model
#
#
#

TESTING=#Default is off
OBJDIR=../obj
EXECDIR=../bin

ifdef DEF_UM
deps=-v UM_source
endif

ifdef DEF_4A
deps2=-v 4A_code
endif

include preprocess

all:    preprocess
	./mkdepend.py $(deps) $(deps2) -a $(OBJDIR)/
	make -f makefile.inc TESTING=$(TESTING) COMPILER=$(COMPILER) \
		 OBJDIR=$(OBJDIR) EXECDIR=$(EXECDIR) all

clean:
	rm -f $(EXECDIR)/*.exe $(OBJDIR)/*.so $(OBJDIR)/*.o \
	*.mod *~ .includes \
	*/ppdone */o */*.o */*.mod UM_source/*.f90 */*.done 
	rm -rf $(OBJDIR)/[^.]*
