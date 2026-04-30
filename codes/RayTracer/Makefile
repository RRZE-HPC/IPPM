CC = mpiicc #Intel Classical Compilers: deprecated
LD = mpiicc #Intel Classical Compilers: deprecated
#CC = mpiicx #Intel LLVM compilers for now and future
#LD = mpiicx #Intel LLVM compilers for now and future


CFLAGS	=  -O3 -xHOST -fno-alias
LDFLAGS = -O3 -xHOST -fno-alias 

all: ray_par

ray_par: ray_par.o
	$(LD) $(LDFLAGS) -o ray_par ray_par.o

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o ray_par

