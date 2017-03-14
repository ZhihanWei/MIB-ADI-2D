# ----- makefile to compile the program heat2dll ----- #

#++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#               compilation environment                #
#   HOST: pc,dmc,palmetto; COMPILER: gfortran,ifort;   #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++#
HOST     = pc
#PROGRAM  = mib2dll
#COMPILER = gfortran
COMPILER = ifort
VPATH    = ./src
LIBS     = libs/libarpack.a libs/libslatec.a

#++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#                    code setups                       #
#             EXAMPLE: 1-7; SOLVER: 1-4                #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++#
EXAMPLE  = 5
SOLVER   = 2
DT       = 1.0E-6
PROGRAM  = mib2d_e$(EXAMPLE)s$(SOLVER)

#++++++++++++ host and compiler selection +++++++++++++#
ifeq ($(HOST),pc)
   ifeq ($(COMPILER),gfortran)
      CFLAG = -g -Wall -fmessage-length=0 -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8
   endif

   ifeq ($(COMPILER),ifort)
      CFLAG = -g -i8 -r8
   endif
endif

ifeq ($(HOST),dmc)
   COMPILER  = ifort
   CFLAG     = -g -i8 -r8
endif

#+++++++++++ source code modules selection ++++++++++++#
DATA = data_comdata.F90
MIB  = mod_mib2d.F90
TIME = mod_time.F90
MAIN = program_main.F90

DATAOBJ = $(DATA:.F90=.o)
MIBOBJ  = $(MIB:.F90=.o)
TIMEOBJ = $(TIME:.F90=.o)
MAINOBJ = $(MAIN:.F90=.o)

#+++++++++++++++++ compilation options ++++++++++++++++#
default:heat2d

heat2d:$(DATAOBJ) $(MIBOBJ) $(TIMEOBJ) $(MAINOBJ) $(DATA) $(MIB) $(TIME) $(MAIN)
	#echo $(COMPILER)
	#echo $(CFLAG)
	#echo $(VPATH)
	#echo $(LIBS)
	$(COMPILER) $(CFLAG) -o $(PROGRAM) $(DATAOBJ) $(MIBOBJ) $(TIMEOBJ) $(MAINOBJ)
	@echo "     "
	@echo ">>> compiled on `hostname -s` with $(COMPILER) <<<"
	@echo "     "

dmc:$(DATAOBJ) $(MIBOBJ) $(TIMEOBJ) $(MAINOBJ) $(DATA) $(MIB) $(TIME) $(MAIN)
	module purge
	module load intel/11.1.072 lapack/3.4.2
	$(COMPILER) $(CFLAG) -o $(PROGRAM) $(DATAOBJ) $(MIBOBJ) $(TIMEOBJ) $(MAINOBJ) $(LIBS)

run:
	./$(PROGRAM)

clean:
	rm -r *.o *.mod
wipe:
	rm -i $(PROGRAM) test*

%.o:%.f90
	$(COMPILER) $(CFLAG) -c $(VPATH)/$*.f90 
	
%.o:%.F90
	$(COMPILER) $(CFLAG) -c -Dexample=$(EXAMPLE) -Dtime=$(SOLVER) $(VPATH)/$*.F90
