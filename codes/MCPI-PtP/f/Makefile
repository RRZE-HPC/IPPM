FC = mpiifort #Intel Classical Compilers: deprecated
CC = mpiicc #Intel Classical Compilers: deprecated
LD = mpiifort #Intel Classical Compilers: deprecated
#FC = mpiifx #Intel LLVM compilers for now and future
#CC = mpiicx #Intel LLVM compilers for now and future
#LD = mpiifx #Intel LLVM compilers for now and future

# flags with optimization
FFLAGS  = -O3 -traceback
LDFLAGS = -O3 -traceback
#FFLAGS  = -O0 -C -g  -traceback
#LDFLAGS = -O0 -C -g  -traceback

# # flags for debugging
# FFLAGS  = -O0 -g -debug extended -fstack-protector-all -debug inline-debug-info -traceback -warn all -C -check all,nooutput_conversion,noformat
# LDFLAGS = -O0 -g -debug extended -fstack-protector-all -debug inline-debug-info -traceback -warn all -C -check all,nooutput_conversion,noformat

DEPS = $(MAKEFILE_LIST)

OBJS = div.o mcpi-ptp.o

.PHONY: all
all: mcpi-ptp

div.c: ../utils/div.c div.h
	[ -L div.c ] || (ln -s ../utils/div.c div.c)

div.h: ../utils/div.h
	[ -L div.h ] || (ln -s ../utils/div.h div.h)

mcpi-ptp: $(OBJS) $(DEPS)
	$(LD) $(LDFLAGS) -o $@ $(filter-out $(DEPS),$^)

mcpi-ptp.o: div.o

%.o: %.F90 $(DEPS)
	$(FC) $(FFLAGS) -c -o $@ $<

%.o: %.c $(DEPS)
	$(CC) $(FFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	-rm -fv *.o *.mod mcpi-ptp
	[ ! -L div.c ] || rm div.c
	[ ! -L div.h ] || rm div.h
