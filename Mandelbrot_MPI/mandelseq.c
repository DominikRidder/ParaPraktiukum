/*
 * mandelseq.c
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include <ppmwrite.h>
#include <mpi.h>

#define MASTER_BLOCKSIZE 16
#define ROOT 0
#define TAG 0

double esecond(void) {

  struct timeval tp;
  struct timezone tzp;
  
  gettimeofday(&tp, &tzp);
  return tp.tv_sec + (tp.tv_usec * 1e-6);
}

void usage(char *name) {
  fprintf(stderr, "Usage: %s options\n\nwith the following optional options (default values in parathesis):\n\n",name);

  fprintf(stderr, "  [-x <x0> <x1> <y0> <y1>]  coordinates of initial area (-1.5 0.5 -1.0 1.0)\n");
  fprintf(stderr, "  [-w <width>]              image width in pixels (256)\n");
  fprintf(stderr, "  [-h <height>]             image height in pixels (256)\n");
  fprintf(stderr, "  [-i <maxiter>]            max. number of iterations per pixel (256)\n");
  fprintf(stderr, "  [-t <type>]               0=stride, 1=stripe, 2=blockmaster, 3=serial\n");
  fprintf(stderr, "  [-v]                      verbose (off)\n\n");
  exit(1);
}

void calc(int *iterations, int width, int height, int myid, int numprocs,
	  double xmin, double xmax, double ymin, double ymax, int maxiter );
void calcRow(int *iterations, int *recvbuf, int width, int height, int myid, int numprocs,
	  double xmin, double xmax, double ymin, double ymax, int maxiter, double *mpitime, double *waittime );
void calcBlock(int *iterations, int *recvbuf, int width, int height, int myid, int numprocs,
	  double xmin, double xmax, double ymin, double ymax, int maxiter, double *mpitime, double *waittime );
void calcMaster(int *iterations, int *recvbuf, int width, int height, int myid, int numprocs,
	  double xmin, double xmax, double ymin, double ymax, int maxiter, double *mpitime, double *waittime );


int main(int argc, char *argv[]) {
  /* default values for command line parameters */

  double xmin = -1.5;  /* coordinates of rectangle */
  double xmax =  0.5;
  double ymin = -1.0;
  double ymax =  1.0;
  int width   = 256;   /* size of rectangle in pixels */
  int height  = 256;
  int maxiter = 256;   /* max. number of iterations */
  int verbose = 0;     /* per default only print error messages */
  int type = 0;        /* per default only print error messages */
  int *iterations,*recvbuffer;
  int ix,iy;

  int    numprocs,myid;
  int    i;
  double st,st_total,timeused;
  double calctime=0.0, waittime=0.0, iotime=0.0, mpitime=0.0, runtime=0.0;
  char   filename[1024];
  
  st_total = esecond();
  MPI_Init(&argc, &argv);
  
  ppminitsmooth(1);

  /* parse command line */
  i=1;
  while( i < argc ) {
    if( argv[i][0] == '-' ) {
      switch( argv[i][1] ) {
      case 'x':
        xmin = atof(argv[++i]);
        xmax = atof(argv[++i]);
        ymin = atof(argv[++i]);
        ymax = atof(argv[++i]);
        break;
      case 'i':
        maxiter = atoi(argv[++i]);
        break;
      case 'w':
        width = atoi(argv[++i]);
        break;
      case 'h':
        height = atoi(argv[++i]);
        break;
      case 't':
        type = atoi(argv[++i]);
        break;
      case 'v':
        verbose++;
        break;
      default:
        usage(argv[0]);
      }
    } else {
      usage(argv[0]);
    }
    i++;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  
  /* initialize recvbuffer */
  if (myid == ROOT) {
    recvbuffer = malloc(sizeof(int) * width * height);
    
    for (ix=0; ix<width; ++ix) {
      for (iy=0; iy<height; ++iy) {
	recvbuffer[ix*height+iy] = 0;
      }
    }
  }

  /* start calculation */
  if(verbose) {
    printf("start calculation (x=%8.5g ..%8.5g,y=%10.7g ..%10.7g) ... \n",
           xmin,xmax,ymin,ymax);
    fflush(stdout);
  }

  st = esecond();
  switch(type) {
    case 0:
      calcRow(iterations, recvbuffer, width, height, myid, numprocs, xmin, xmax, ymin, ymax, maxiter, &mpitime, &waittime ); break;
    case 1:
      calcBlock(iterations, recvbuffer, width, height, myid, numprocs, xmin, xmax, ymin, ymax, maxiter, &mpitime, &waittime ); break;
    case 2:
      calcMaster(iterations, recvbuffer, width, height, myid, numprocs, xmin, xmax, ymin, ymax, maxiter, &mpitime, &waittime ); break;
    case 3:
      if (myid == ROOT) { // single + serial
	calc(recvbuffer, width, height, myid, numprocs, xmin, xmax, ymin, ymax, maxiter );
      }
      break;
  default:
      printf("Please use a type from 0-3!\n");
      MPI_Abort(MPI_COMM_WORLD, type);
  }
  timeused = esecond()-st;
  calctime += timeused;

  timeused = esecond()-st_total;
  runtime += timeused;

  st = esecond();
  if (myid == ROOT) {
    ppmwrite(recvbuffer,width,height,0,maxiter,"mandelcol.ppm");
  }

  timeused = esecond()-st;
  iotime += timeused;
  if(verbose) printf("PE%02d: calc=%7.4f, mpi=%7.4f, wait=%7.4f, io=%7.4f, run=%7.4f\n",
                     myid,calctime,mpitime,waittime,iotime,runtime);

  MPI_Finalize();
  exit(0);
}

void calc(int *iterations, int width, int height, int myid, int numprocs,
         double xmin, double xmax, double ymin, double ymax, int maxiter ) {
  double dx,dy,x,y;
  int    ix,iy;

  dx = (xmax - xmin) / width;
  dy = (ymax - ymin) / height;

  y = ymin;
  for (iy=0; iy<height; ++iy) {
    x = xmin;
    for (ix=0; ix<width; ix++) {
      double zx=0.0,zy=0.0,zxnew;
      int count = 0;
      while ( zx*zx+zy*zy < 16*16 && count < maxiter ) {
        /* z = z*z + (x + i y) */
        zxnew = zx*zx-zy*zy + x;
        zy    = 2*zx*zy     + y;
        zx    = zxnew;
        ++count;
      }
      iterations[iy*width+ix] = count;
      x += dx;
    }
    y += dy;
  }
}

void calcRow(int *iterations, int *recvbuf, int width, int height, int myid, int numprocs,
	     double xmin, double xmax, double ymin, double ymax, int maxiter, double *mpitime, double *waittime ) {
  double dx,dy,x,y;
  int    ix,iy;
  double st;

  iterations = malloc(sizeof(int) * width * height);

  for (int i = 0; i < width * height; i++) {
    iterations[i] = 0;
  }

  dx = (xmax - xmin) / width;
  dy = (ymax - ymin) / height;

  y = ymin + dy*myid;
  for (iy=myid; iy<height; iy+=numprocs) {
    x = xmin;
    for (ix=0; ix<width; ix++) {
      double zx=0.0,zy=0.0,zxnew;
      int count = 0;
      while ( zx*zx+zy*zy < 16*16 && count < maxiter ) {
        /* z = z*z + (x + i y) */
        zxnew = zx*zx-zy*zy + x;
        zy    = 2*zx*zy     + y;
        zx    = zxnew;
        ++count;
      }
      iterations[iy*width+ix] = count;
      x += dx;
    }
    y += dy*numprocs;
  }

  st = esecond();
  MPI_Barrier(MPI_COMM_WORLD);
  *waittime = esecond()-st;

  st = esecond();
  MPI_Reduce(iterations, recvbuf, width * height, MPI_INT, MPI_SUM, ROOT, MPI_COMM_WORLD);
  *mpitime = esecond()-st;
}

void calcBlock(int *iterations, int *recvbuf, int width, int height, int myid, int numprocs,
	       double xmin, double xmax, double ymin, double ymax, int maxiter, double *mpitime, double *waittime ) {
  double dx,dy,x,y;
  int    ix,iy;
  double st;
  int *recvcounts, *displs;

  dx = (xmax - xmin) / width;
  dy = (ymax - ymin) / height;

  recvcounts = (int *) malloc(sizeof(int) * numprocs);
  displs = (int *) malloc(sizeof(int) * numprocs);
    
  for (int i = 0; i < numprocs; i++) {
    recvcounts[i] = height / numprocs * width;
    if (i < height % numprocs) {
      recvcounts[i] += width;
    }
    if (i > 0) {
      displs[i] = displs[i-1] + recvcounts[i-1];
    } else {
      displs[i] = 0;
    }
  }
  
  iterations = (int *) malloc(sizeof(int) * recvcounts[myid]);

  y = ymin + dy * displs[myid]/width;
  for (iy=0; iy<recvcounts[myid]/width; ++iy) {
    x = xmin;
    for (ix=0; ix<width; ix++) {
      double zx=0.0,zy=0.0,zxnew;
      int count = 0;
      while ( zx*zx+zy*zy < 16*16 && count < maxiter ) {
        /* z = z*z + (x + i y) */
        zxnew = zx*zx-zy*zy + x;
        zy    = 2*zx*zy     + y;
        zx    = zxnew;
        ++count;
      }
      iterations[iy*width+ix] = count;
      x += dx;
    }
    y += dy;
  }

  st = esecond();
  MPI_Barrier(MPI_COMM_WORLD);
  *waittime = esecond()-st;

  st = esecond();
  MPI_Gatherv(iterations, recvcounts[myid], MPI_INT, recvbuf, recvcounts, displs, MPI_INT, ROOT, MPI_COMM_WORLD);
  *mpitime = esecond()-st;
}

void calcMaster(int *iterations, int *recvbuf, int width, int height, int myid, int numprocs,
		double xmin, double xmax, double ymin, double ymax, int maxiter, double *mpitime, double *waittime ) {
  double dx,dy,x,y;
  double st;
  int    ix,iy,xoff,yoff;  
  int    blockID = 0;  
  int    nBlocks = (height * width) / (MASTER_BLOCKSIZE * MASTER_BLOCKSIZE);
  int    inwork = 0; // number of SLAVES calculating
  int    cancel = - 1; // SLAVE release msg
  MPI_Status status;

  if (numprocs == 1) {
    printf("There is no slave for master-slave type.\n");
    MPI_Abort(MPI_COMM_WORLD, myid);
  }

  iterations = (int *) malloc(sizeof(int) * (MASTER_BLOCKSIZE * MASTER_BLOCKSIZE + 1)); // block-id + block

  dx = (xmax - xmin) / width;
  dy = (ymax - ymin) / height;

  if (myid == ROOT) { // MASTER

    for (int i = 1; i < numprocs; i++) {
      //MASTER initial work spreading
      MPI_Send(&blockID, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      blockID++;
      inwork++;
    } 

    while(inwork != 0) {
      //MASTER waiting for more slaves
      st = esecond();
      MPI_Recv(iterations, MASTER_BLOCKSIZE * MASTER_BLOCKSIZE + 1, MPI_INT, MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, &status);
      *mpitime += esecond()-st;
      inwork--;

      xoff = (iterations[0] * MASTER_BLOCKSIZE % width);
      yoff = (iterations[0] * MASTER_BLOCKSIZE / width) * MASTER_BLOCKSIZE * width;
      for (iy=0; iy<MASTER_BLOCKSIZE; ++iy) {
	for (ix=0; ix<MASTER_BLOCKSIZE; ix++) {
	  recvbuf[yoff+xoff+iy*width+ix] = iterations[iy*MASTER_BLOCKSIZE+ix+1];
	}
      }


      if (blockID < nBlocks) {
	MPI_Send(&blockID, 1, MPI_INT, status.MPI_SOURCE, TAG, MPI_COMM_WORLD);
	blockID++;
	inwork++;
      } else {
	MPI_Send(&cancel, 1, MPI_INT, status.MPI_SOURCE, TAG, MPI_COMM_WORLD);
      }
    }
  } else { // SLAVE
    while(1) {
      MPI_Recv(&blockID, 1, MPI_INT, ROOT, 0, MPI_COMM_WORLD, &status);

      if (blockID == cancel) {
	//SLAVE canceld
	break;
      }

      y = (blockID * MASTER_BLOCKSIZE / width) * MASTER_BLOCKSIZE * dy + ymin;
       for (iy=0; iy<MASTER_BLOCKSIZE; ++iy) {
	 x = ((blockID * MASTER_BLOCKSIZE) % width) * dx + xmin;
	 for (ix=0; ix<MASTER_BLOCKSIZE; ix++) {
	   double zx=0.0,zy=0.0,zxnew;
	   int count = 0;
	   while ( zx*zx+zy*zy < 16*16 && count < maxiter ) {
	     /* z = z*z + (x + i y) */
	     zxnew = zx*zx-zy*zy + x;
	     zy    = 2*zx*zy     + y;
	     zx    = zxnew;
	     ++count;
	   }
	   iterations[iy*MASTER_BLOCKSIZE+ix+1] = count;
	   x += dx;
	 }
	 y += dy;
       }

       iterations[0] = blockID;

       st = esecond();
       MPI_Send(iterations, MASTER_BLOCKSIZE * MASTER_BLOCKSIZE + 1, MPI_INT, ROOT, TAG, MPI_COMM_WORLD);
       *mpitime += esecond()-st;
    }
  }

  st = esecond();
  MPI_Barrier(MPI_COMM_WORLD);
  *waittime = esecond()-st;
}
