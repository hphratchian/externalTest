#
# This is a simple makefile for building spin-squared calculation code.
#
MQCDir       = $(mqcinstall)
MQCMODS      = $(MQCDir)/PGI/mod
MQCLIB       = $(MQCDir)/PGI/lib
LIBS         = -llapack -lblas -L$(MQCLIB)
F03Flags     = 
RunF         = pgfortran -i8 -r8 -Mallocatable=03
#RunF         = pgfortran -i8 -r8
#
#
# The 'all' rule.
#
all: external.exe

#
# Generic rules for building module (*.mod) and object (*.o) files.
#
#%.mod: %.f90
#	$(RunF) -c $*.f90

%.o: %.f90
	$(RunF) -I$(MQCMODS) -c $*.f90

%.o: %.f03
	$(RunF) $(F03Flags) -I$(MQCMODS) -c $*.f03

#
# Generic rule for building general executable program (*.exe) from a standard
# f90 source (*.f90) file.
#

%.exe: %.f03 %_mod.f03 $(MQCLIB)/libmqc.a
	$(RunF) $(LIBS) $(Prof) -I$(MQCMODS) -o $*.exe $*.f03 $(MQCLIB)/libmqc.a

#
# Some clean rules.
#
cleanGTests:
	(cd GTests ; rm -f *.chk ; rm -f *.mat ; rm -f *.log)	
