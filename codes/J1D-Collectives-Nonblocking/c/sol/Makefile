CC = mpiicpc #Intel Classical Compilers: deprecated
LD = mpiicpc #Intel Classical Compilers: deprecated
#CC = mpiicpx #Intel LLVM compilers for now and future
#LD = mpiicpx #Intel LLVM compilers for now and future

# flags with optimization
CFLAGS  = -O3 -traceback -fno-alias
LDFLAGS = -O3 -traceback -fno-alias
#CFLAGS  = -O0 -g -traceback
#LDFLAGS = -O0 -g -traceback

# # flags for debugging
# CFLAGS  = -O0 -g -debug extended -fstack-protector-all -debug inline-debug-info -traceback -warn all -C -check all,nooutput_conversion,noformat
# LDFLAGS = -O0 -g -debug extended -fstack-protector-all -debug inline-debug-info -traceback -warn all -C -check all,nooutput_conversion,noformat

DEPS = $(MAKEFILE_LIST)

OBJS = utils.o j1d_collectives_nonblocking.o

.PHONY: all
all: j1d_collectives_nonblocking

utils.cpp: utils.h
	[ -L utils.cpp ] || (ln -s ../utils.cpp utils.cpp)

utils.h:
	[ -L utils.h ] || (ln -s ../utils.h utils.h)

j1d_collectives_nonblocking: $(OBJS) $(DEPS)
	$(LD) $(LDFLAGS) -o $@ $(filter-out $(DEPS),$^)

j1d_collectives_nonblocking.o: utils.o

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	-rm -fv *.o *.mod j1d_collectives_nonblocking
	[ ! -L utils.cpp ] || rm utils.cpp
	[ ! -L utils.h ] || rm utils.h
