# set the appropriate resolutions for building each test case.
# 
# The resolutions chosen here are as specified for any intercomparison

#
# Uncomment the appropriate choice below or specify in the call to make
# i.e. make CASE=CU_2D all
#
CASE=1D
#CASE=SC_2D
#CASE=CU_2D
#CASE=ISDAC_2D
#CASE=SQUALL_2D
#CASE=OROG_2D

ifeq ($(CASE), 1D)
D1=DEF_NX=1
D2=DEF_NZ=120
endif

ifeq ($(CASE), SC_2D)
D1=DEF_NX=100
D2=DEF_NZ=52
endif

ifeq ($(CASE), CU_2D)
D1=DEF_NX=180
D2=DEF_NZ=60
endif

ifeq ($(CASE), ISDAC_2D)
D1=DEF_NX=40
D2=DEF_NZ=26
endif

ifeq ($(CASE), SQUALL_2D)
D1=DEF_NX=320
D2=DEF_NZ=48
endif

ifeq ($(CASE), OROG_2D)
D1=DEF_NX=200
D2=DEF_NZ=120
endif

ifeq ($(CASE), ICMW_SC)
D1=DEF_NX=60
D2=DEF_NZ=61
endif

ifeq ($(CASE), CMW_BOX)
D1=DEF_NX=1
D2=DEF_NZ=2
endif

# Allow input at the command line
COMMANDLINE=1
# Use the following if it doesn't work with your compiler
#ifeq ($(COMPILER), ????)
#COMMANDLINE=0
#endif

# extra code components to include
DEF_UM=#Default is off
DEF_4A=#Default is off

ifdef DEF_UM
UM_MICRO=1
else
UM_MICRO=0
endif

ifdef DEF_4A
SHIPWAY_MICRO=1
else
SHIPWAY_MICRO=0
endif

# Use allocatable arrays within derived types instead of pointers
ALLOCATABLE_TYPE=0

# Defines to pass to compiler
DEFS=-D$(D1) -D$(D2) -DDEF_CASE=\'$(CASE)\' -DCOMMANDLINE=$(COMMANDLINE) -DMODEL_KiD=1 -DMODEL_LEM=2 -DMODEL_UM=3 -DDEF_MODEL=MODEL_KiD -DALLOCATABLE_TYPE=$(ALLOCATABLE_TYPE) -DSHIPWAY_MICRO=$(SHIPWAY_MICRO) -DUM_MICRO=$(UM_MICRO)


