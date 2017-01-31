Please see the documentation in doc/KiD.pdf for further details.

Updates for KiD version 2.1.1804+
====================================

1) Update to a 2D driver

Version 2.X of the KiD allows for 2D flow fields to be specified. The
driver can still be used as a single column as before, but now it has
the option for 2D (by setting the number of horizontal gridpoints,
nx). Since many of the coupled microphysics schemes require fixed
dimensions, KiD must be recompiled for each different configuration (see
preprocessing and build below)
Please see the documentation for more details on 2D KiD.

2) Preprocessing, build and running

(Note the following assumes gnu make)

The makefile now has a few more options and by default preprocesses all
fortran source files.  For a standard build of the original 1D model,
once you have unpacked the model simply cd into the top level directory
(NOT /src) and type 

make CASE=1D all
./bin/KiD_1D.exe

This will run the model with the default namelist. By default the ifort
compiler is used, to change this modify COMPILER in the top level
makefile (note COMPILERS is used for building on several compilers -
useful for developing new code).  You can also specify this at the
command line, i.e.

make COMPILER=gfortran CASE=1D all

The preprocessor directives that define the grid are passed through to
the model from src/defines.inc. It is recommended that you add in new
definitions if wishing to run at different resolutions.

There are currently 5 different case configurations defined;
1D, SC_2D, CU_2D, ISDAC_2D, SQUALL_2D
to build exectuables for them all, use

make build_cases

(again a different compiler can be specified as above).  To run with all
the default namelists, use

make run_cases

This should run the cases sequentially and  produce netcdf files (and
corresponding namelists) in output (Note that 2D cases take a lot longer
to run than previous 1D cases). If you wish to exploit multiple cpus on
your system (currently no OMP support in KiD), you can use the -j
option, i.e. 

make -j run_cases

(use the -l option to limit load on the system).

As with previous versions if your compiler allows it, custom namelists
and output can be specified at the command line when running, e.g.

./bin/KiD_SC_2D.exe namelists/SC_2D.nml output/my_output_file.nc

If no command line arguments are given, default namelists will be used.
If no output file is specified, the default naming convention for the
output will be used.


3) 2D Test cases

As noted above there are 4 new 2D test cases defined:

SC_2D:     An eddy-pair based on Morrison and Grabowski (2007)
CU_2D:     A time varying convective case based on Morrison and Grabowski
           (2007)
ISDAC_2D:  Similar to SC_2D, but modified for a cold environment
SQUALL_2D: A flow field designed to mimic a squall-line (similar to 
           Slawinska et al 2009, QJRMS)

The latter two are not well tested or finalized, so are subject to
change in later releases.


