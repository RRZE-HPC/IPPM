FC = mpiifort #Intel Classical Compilers: deprecated
LD = mpiifort #Intel Classical Compilers: deprecated
#FC = mpiifx #Intel LLVM compilers for now and future
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

OBJS = utils.o j1d_collectives_nonblocking.o

.PHONY: all
all: j1d_collectives_nonblocking


j1d_collectives_nonblocking: $(OBJS) $(DEPS)
	$(LD) $(LDFLAGS) -o $@ $(filter-out $(DEPS),$^)

j1d_collectives_nonblocking.o: utils.o

%.o: %.F90 $(DEPS)
	$(FC) $(FFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	-rm -fv *.o *.mod j1d_collectives_nonblocking
