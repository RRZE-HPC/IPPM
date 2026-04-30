FC = mpiifort
LD = mpiifort

# flags with optimization
FFLAGS  = -O3 -traceback
LDFLAGS = -O3 -traceback
#FFLAGS  = -O0 -C -g  -traceback
#LDFLAGS = -O0 -C -g  -traceback

# # flags for debugging
# FFLAGS  = -O0 -g -debug extended -fstack-protector-all -debug inline-debug-info -traceback -warn all -C -check all,nooutput_conversion,noformat
# LDFLAGS = -O0 -g -debug extended -fstack-protector-all -debug inline-debug-info -traceback -warn all -C -check all,nooutput_conversion,noformat

DEPS = $(MAKEFILE_LIST)

OBJS = utils.o Jacobi_1D.o

.PHONY: all
all: jacobi_1d

utils.F90:
	[ -L utils.F90 ] || (ln -s ../utils.F90 utils.F90)

jacobi_1d: $(OBJS) $(DEPS)
	$(LD) $(LDFLAGS) -o $@ $(filter-out $(DEPS),$^)

Jacobi_1D.o: utils.o

%.o: %.F90 $(DEPS)
	$(FC) $(FFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	-rm -fv *.o *.mod jacobi_1d
	[ ! -L utils.F90 ] || rm utils.F90
