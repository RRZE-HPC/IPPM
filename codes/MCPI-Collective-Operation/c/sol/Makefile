CC = mpiicc #Intel Classical Compilers: deprecated
LD = mpiicc #Intel Classical Compilers: deprecated
#CC = mpiicx #Intel LLVM compilers for now and future
#LD = mpiicx #Intel LLVM compilers for now and future

# flags with optimization
CFLAGS  = -O3 -traceback -fno-alias
LDFLAGS = -O3 -traceback -fno-alias
#CFLAGS  = -O0 -C -g  -traceback
#LDFLAGS = -O0 -C -g  -traceback

# # flags for debugging
# CFLAGS  = -O0 -g -debug extended -fstack-protector-all -debug inline-debug-info -traceback -warn all -C -check all,nooutput_conversion,noformat
# LDFLAGS = -O0 -g -debug extended -fstack-protector-all -debug inline-debug-info -traceback -warn all -C -check all,nooutput_conversion,noformat

DEPS = $(MAKEFILE_LIST)

OBJS = div.o mcpi-collective-operation.o

.PHONY: all
all: mcpi-collective-operation

div.c: ../../utils/div.c div.h
	[ -L div.c ] || (ln -s ../../utils/div.c div.c)

div.h: ../../utils/div.h
	[ -L div.h ] || (ln -s ../../utils/div.h div.h)

mcpi-collective-operation: $(OBJS) $(DEPS)
	$(LD) $(LDFLAGS) -o $@ $(filter-out $(DEPS),$^)

mcpi-collective-operation.o: div.o

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	-rm -fv *.o *.mod mcpi-collective-operation
	[ ! -L div.c ] || rm div.c
	[ ! -L div.h ] || rm div.h
