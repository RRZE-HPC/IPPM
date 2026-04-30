CC = mpiicc
LD = mpiicc

# flags with optimization
FFLAGS  = -O3 -traceback -fno-alias
LDFLAGS = -O3 -traceback -fno-alias
#FFLAGS  = -O0 -g -traceback
#LDFLAGS = -O0 -g -traceback

# # flags for debugging
# FFLAGS  = -O0 -g -debug extended -fstack-protector-all -debug inline-debug-info -traceback -warn all -C -check all,nooutput_conversion,noformat
# LDFLAGS = -O0 -g -debug extended -fstack-protector-all -debug inline-debug-info -traceback -warn all -C -check all,nooutput_conversion,noformat

DEPS = $(MAKEFILE_LIST)

OBJS = utils.o Jacobi_1D.o

.PHONY: all
all: jacobi_1d

utils.cpp: utils.h
	[ -L utils.cpp ] || (ln -s ../utils.cpp utils.cpp)

utils.h:
	[ -L utils.h ] || (ln -s ../utils.h utils.h)

jacobi_1d: $(OBJS) $(DEPS)
	$(LD) $(LDFLAGS) -o $@ $(filter-out $(DEPS),$^)

Jacobi_1D.o: utils.o

%.o: %.cpp $(DEPS)
	$(CC) $(FFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	-rm -fv *.o *.mod jacobi_1d
	[ ! -L utils.cpp ] || rm utils.cpp
	[ ! -L utils.h ] || rm utils.h
