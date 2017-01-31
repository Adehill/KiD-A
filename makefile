#
# Makefile for KiD model
#
#
#

PRG_ENV=default
TESTING=#Default is off

CASES=1D SC_2D CU_2D ISDAC_2D SQUALL_2D WMO_CASE1# OROG_2D
COMPILERS=nagfor ifort pgf90 gfortran

include src/compiler_options.inc

EXECDIR=bin
NMLDIR=namelists
OUTDIR=output
OBJDIR=obj

$(COMPILERS): build_compilers
	x=$@
	OUTDIR=output/"$$x"_$(TESTING); \
	NMLDIR=namelists; \
	EXECDIR=bin_"$$x"_$(TESTING);\
	OBJDIR=obj_"$$x"_$(TESTING);\
	make COMPILER="$$x" EXECDIR="$$EXECDIR" \
	     NMLDIR="$$NMLDIR" OUTDIR="$$OUTDIR" run_cases; \

build_compilers:
	for x in $(COMPILERS); do \
	make clean; \
	OUTDIR=output/"$$x"_$(TESTING); \
	NMLDIR=namelists; \
	EXECDIR=bin_"$$x"_$(TESTING);\
	OBJDIR=obj_"$$x"_$(TESTING);\
	make TESTING=$(TESTING) COMPILER="$$x" EXECDIR="$$EXECDIR" \
	     NMLDIR="$$NMLDIR" OUTDIR="$$OUTDIR" OBJDIR="$$OBJDIR" \
		build_cases; \
	done

$(CASES): %: $(EXECDIR)/KiD_%.exe
	sed -e "s|fileIn=.*|fileIn=\'$(NMLDIR)/$@.nml\'|" \
	 -e "s|fileOut=.*|fileOut=\'$(OUTDIR)/$@_out.nc\'|" namelists/input.nml > $(NMLDIR)/$@_input.nml;
	$<

run_cases:  $(CASES)

build_cases:
	for x in $(CASES); do \
	make TESTING=$(TESTING) COMPILER=$(COMPILER) CASE="$$x" \
		EXECDIR=$(EXECDIR) NMLDIR=$(NMLDIR) OUTDIR=$(OUTDIR) \
		all; \
	done

all: force_look	
	if [ "$(PRG_ENV)" != "default" ]; then . $(PRG_ENV); fi 
	[ -d $(EXECDIR) ] || mkdir -p "$(EXECDIR)"; \
	[ -d $(OBJDIR) ] || mkdir -p "$(OBJDIR)"; \
	[ -d $(OUTDIR) ] || mkdir -p "$(OUTDIR)"; \
	[ -d $(NMLDIR) ] || mkdir -p "$(NMLDIR)"; \
	cd src; make TESTING=$(TESTING) COMPILER=$(COMPILER) \
			OBJDIR=../$(OBJDIR) EXECDIR=../$(EXECDIR) all

doc: force_look
	cd doc; make all

clean: force_look
	cd src; make OBJDIR=../$(OBJDIR) EXECDIR=../$(EXECDIR) clean
	rm -f src/case.used

doc_clean: force_look
	cd doc; make OBJDIR=$(OBJDIR) EXECDIR=$(EXECDIR) clean

really_clean: clean doc_clean
	rm -f output/*.nc \
	output/*.nml \
	run_data/*.data \
	output/figures/*

force_look:
	true
