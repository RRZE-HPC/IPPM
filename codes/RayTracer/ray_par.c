#include <stdio.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>

double objs[] = {
  0. ,0., -100.5, 10000., 0., 0., 0., 0.25,
  0.272166, 0.272166, 0.544331, .027777,
  0.643951, 0.172546, 0., .027777,
  0.172546, 0.643951, 0., .027777,
 -0.371785, 0.099620, 0.544331, .027777,
 -0.471405, 0.471405, 0., .027777,
 -0.643951,-0.172546, 0., .027777,
  0.099620,-0.371785, 0.544331, .027777,
 -0.172546,-0.643951, 0., .027777,
  0.471405,-0.471405, 0., .027777,
  4., 3., 2., 1., -4., 4., -3., 1., 5.
};

double *intersect(double x, double y, double z, double dx, double dy, double dz, double *maxp)
{
  int i;
  double *o = objs, *oo = 0;
  double max = *maxp;
  double xx, yy, zz, b, t;

  for (i = 0; i < 11; i++)
    {
      xx = *o++ - x; yy = *o++ - y; zz = *o++ - z;
      b = xx * dx + yy * dy + zz * dz;
      if ((t = b * b - xx * xx - yy * yy - zz * zz + *o++) < 0 || (t = b-sqrt(t)) < 1e-6 || t > max)
	continue;
      oo = o - 4;
      max = t;
    }
  *maxp = max;
  return oo;
}

double shade(double x, double y, double z, double dx, double dy, double dz, int de)
{
  double max = 1e6, c = 0, r, k, *o;
  double nx, ny, nz, ldx, ldy, ldz, rdx, rdy, rdz;
  int i;

  if (!(o = intersect(x, y, z, dx, dy, dz, &max)))
    return 0;
  x += max * dx; y += max * dy; z += max * dz;
  nx = x - *o++; ny = y - *o++; nz = z - *o++;
  r = sqrt(nx * nx + ny * ny + nz * nz);
  nx /= r; ny /= r; nz /= r;
  k = nx * dx + ny * dy + nz * dz;
  rdx = dx - 2 * k * nx; rdy = dy - 2 * k * ny; rdz = dz - 2 * k * nz;
  o = objs + 44;
  for (i = 0; i < 3; i++)
    {
      ldx = *o++ - x; ldy = *o++ - y; ldz = *o++ - z;
      r = sqrt(ldx * ldx + ldy * ldy + ldz * ldz);
      ldx /= r; ldy /= r; ldz /= r;
      if (intersect(x, y, z, ldx, ldy, ldz, &r))
	continue;
      if ((r = ldx * nx + ldy * ny + ldz * nz) < 0)
	continue;
      c += r;
      if ((r = rdx * ldx + rdy * ldy + rdz * ldz) > 0)
	c += 2 * pow(r, 15.);
    }
  if (de < 10)
    c += .5 * shade(x, y, z, rdx, rdy, rdz, de + 1);
  return c;
}


void calc_tile(int size, int xstart, int ystart, int tilesize, unsigned char* tile)
{
  double dx, dy, dz, c, r;
  int x, y, i;
  
  i=0;

  for (y = ystart; y < ystart+tilesize; y++)
    for (x = xstart; x < xstart+tilesize; x++)
      {
        double xx = x / (float)(size-1);
        double yy = 1. - y / (float)(size-1);
        dx = -0.847569 - xx * 1.30741 - yy * 1.19745;
        dy = -1.98535  + xx * 2.11197 - yy * 0.741279;
        dz = -2.72303                 + yy * 2.04606;
        r = sqrt(dx * dx + dy * dy + dz * dz);
	c = 100 * shade(2.1, 1.3, 1.7, dx / r, dy / r, dz / r, 0);
	if (c < 0)
	  c = 0;
	if (c > 255)
	  c = 255;
	tile[(y-ystart)*tilesize+(x-xstart)]=(unsigned char)c;
      }
}

int main(int argc, char** argv)
{
  int size; //10000
  int tilesize; //=100
  int xtiles,ytiles,xc,yc,count,tilebase,i;

  unsigned char *tile, *picture;
  char c;
  int my_rank, p, source, rcount, scount, ts1, ts2;
  MPI_Status status;
  int **tiledata;
  int rtiledata[4];
  FILE *fd;
  char string[100];

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  if(argc<3) {
      printf("ERROR: two arguments (size and tilesize) are required!\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  else {
      size=atoi(argv[1]);
      tilesize=atoi(argv[2]);
  }
  if(size%tilesize!=0) {
      printf("ERROR: size is not a multiple of tilesize!\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }

  if((tile=(unsigned char*)malloc(tilesize*tilesize*sizeof(unsigned char)))==NULL)
    {
      fprintf(stderr,"Could not allocate tile memory!\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }      

  if(0==my_rank)
    {
      if((picture=(unsigned char*)malloc(size*size*sizeof(unsigned char)))==NULL)
	{
	  fprintf(stderr,"Could not allocate picture memory!\n");
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}
      tiledata=(int**)malloc(p*sizeof(int*));
      for(i=0; i<p ; i++)
	{
	  tiledata[i]=(int*)malloc(4*sizeof(int));
	}
    }
  
//  fprintf(stderr, "Hello!\n");  
  
  double stime=MPI_Wtime();
  MPI_Request r;
  
  /* number of tiles in x and y direction */
  xtiles=size/tilesize;
  ytiles=size/tilesize;
  
  if(0==my_rank){
    rcount=0; scount=0;
    while(rcount<xtiles*ytiles)
      {
	xc=scount % xtiles;
	yc=scount / xtiles;

	/* wait for a client to request a tile */
	MPI_Recv(tile, tilesize*tilesize, 
		 MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	source=status.MPI_SOURCE;
	ts1=tiledata[source][1]; ts2=tiledata[source][2];
	/* give tile data to client */
	if(scount<xtiles*ytiles)
	  {
	    tiledata[source][0]=size; tiledata[source][1]=xc*tilesize; 
	    tiledata[source][2]=yc*tilesize; tiledata[source][3]=tilesize;
        // usleep(1);
	    MPI_Send(tiledata[source], 4, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
	    scount++;
	  }
	else			/* no tiles available; shut down client with tag=2 */
	  {
	    MPI_Send(tiledata, 4, MPI_INT, status.MPI_SOURCE, 2, MPI_COMM_WORLD);
	  }
	/* was there a tile sent with the request? */
	if(1==status.MPI_TAG)
	  {
	    /* fprintf(stderr, "Hello-tile!\n");   */
	    /* if so, copy to picture buffer */
	    for(i=0; i<tilesize; i++)
	      {
		tilebase=ts2*tilesize*xtiles+ts1;
		memcpy((void*)(picture+tilebase+i*tilesize*xtiles),
		       (void*)(tile+i*tilesize),
		       tilesize*sizeof(unsigned char));
	      }
	    rcount++; /* fprintf(stderr,"[%d]\n",rcount); */
	  }
      }

    stime=MPI_Wtime()-stime;
    printf("Procs: %d Time: %8.3lf   MP/s: %6.1lf   size: %6d    tilesize: %5d\n",p,stime,(double)size*size/stime/1000000.0,size,tilesize);
    fd=fopen("/tmp/ray.pnm","w");
    fprintf(fd,"P5\n%d %d\n255\n",size,size);
    fwrite(picture,sizeof(unsigned char),size*size,fd);
    fclose(fd);
    /* for(count=0; count<size*size; count++)
       putchar(picture[count]); */
    
    for(i=0; i<p ; i++)
      free(tiledata[i]);
    
    free(tiledata);
    free(picture);
  }
  else{
    /* we are a client; first send is a pure request (tag=0) */
    MPI_Send(tile, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    for(;;)
      {
	MPI_Recv(rtiledata, 4, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	/* no tiles available? quit! */
	if(2==status.MPI_TAG)
	  break;
	calc_tile(rtiledata[0],rtiledata[1],rtiledata[2],rtiledata[3],tile);
	/* send tile to master (tag=1) */
	MPI_Send(tile, rtiledata[3]*rtiledata[3], MPI_CHAR, 0, 1, MPI_COMM_WORLD);
      }
  }
  
  free(tile);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  exit(0);
}
