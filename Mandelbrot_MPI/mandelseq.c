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
  fprintf(stderr, "  [-t <type>]               0=stride, 1=stripe, 2=blockmaster\n");
  fprintf(stderr, "  [-v]                      verbose (off)\n\n");
  exit(1);
}

void calc(int *iterations, int width, int height, int myid, int numprocs,
         double xmin, double xmax, double ymin, double ymax, int maxiter );
void calcRow(int *iterations, int *recvbuf, int width, int height, int myid, int numprocs,
	  double xmin, double xmax, double ymin, double ymax, int maxiter );
void calcBlock(int *iterations, int *recvbuf, int width, int height, int myid, int numprocs,
	  double xmin, double xmax, double ymin, double ymax, int maxiter );
void calcMaster(int *iterations, int *recvbuf, int width, int height, int myid, int numprocs,
	  double xmin, double xmax, double ymin, double ymax, int maxiter );


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
  int    ix,iy;

  int    numprocs,myid;
  int    i;
  double st,timeused,calctime=0.0, waittime=0.0, iotime=0.0;
  char   filename[1024];
  
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
  
  /* initialize arrays */
  if (myid == 0) {
    recvbuffer = malloc(width*height*sizeof(int));
    
    for (ix=0; ix<width; ++ix) {
      for (iy=0; iy<height; ++iy) {
	recvbuffer[ix*height+iy] = 0;
      }
    }
  }

  
  if (type == 0 || type == 3) {
    iterations = malloc(width*height*sizeof(int));
    
    for (ix=0; ix<width; ++ix) {
      for (iy=0; iy<height; ++iy) {
	iterations[ix*height+iy] = 0;
      }
    }
  } else if (type == 1) {
    iterations = malloc(height / numprocs * width * sizeof(int));
    
    for (ix=0; ix<width; ++ix) {
      for (iy=0; iy<height / numprocs; ++iy) {
	iterations[iy*height+ix] = 0;
      }
    }
  } else if (type == 2) {
    iterations = malloc(MASTER_BLOCKSIZE * MASTER_BLOCKSIZE * sizeof(int));
    
    for (ix=0; ix<MASTER_BLOCKSIZE; ++ix) {
      for (iy=0; iy<MASTER_BLOCKSIZE; ++iy) {
	iterations[iy*MASTER_BLOCKSIZE+ix] = 0;
      }
    }
  }

  //numprocs = 1;
  //myid     = 0;

  /* start calculation */
  if(verbose) {
    printf("start calculation (x=%8.5g ..%8.5g,y=%10.7g ..%10.7g) ... \n",
           xmin,xmax,ymin,ymax);
    fflush(stdout);
  }
  printf("Start (id = %d, #procs = %d)\n", myid, numprocs);

  st = esecond();
  switch(type) {
    case 0:
      calcRow(iterations, recvbuffer, width, height, myid, numprocs, xmin, xmax, ymin, ymax, maxiter ); break;
    case 1:
      calcBlock(iterations, recvbuffer, width, height, myid, numprocs, xmin, xmax, ymin, ymax, maxiter ); break;
    case 2:
      calcMaster(iterations, recvbuffer, width, height, myid, numprocs, xmin, xmax, ymin, ymax, maxiter ); break;
    case 3:
      if (myid == 0) { // single + serial
	calc(iterations, width, height, myid, numprocs, xmin, xmax, ymin, ymax, maxiter );
      }
      recvbuffer = iterations;
      break;
  default:
      printf("Please use a type from 0-3!\n");
      MPI_Abort(1, MPI_COMM_WORLD);
  }
  timeused = esecond()-st;
  calctime += timeused;

  printf("Returned from calculation\n");

  st = esecond();
  //ppmwrite(iterations,width,height,0,maxiter,"mandelcol.ppm");
  if (myid == 0) {
    ppmwrite(recvbuffer,width,height,0,maxiter,"mandelcol.ppm");
  }

  timeused = esecond()-st;
  iotime += timeused;
  if(verbose) printf("PE%02d: calc=%7.4f,wait=%7.4f, io=%7.4f\n",
                     myid,calctime,waittime,iotime);

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
         double xmin, double xmax, double ymin, double ymax, int maxiter ) {
  double dx,dy,x,y;
  int    ix,iy;
  //MPI_type type;
  
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
  
  //MPI_Type_vector(height / numprocs, width, width * numprocs, MPI_INT, &type);
  //MPI_Type_commit(&type);
  //MPI_Reduce(&iterations[myid], recvbuf, 1, type, MPI_SUM, 0, MPI_COMM_WORLD);
  //MPI_Type_free(&type);
  MPI_Reduce(iterations, recvbuf, width * height, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
}

void calcBlock(int *iterations, int *recvbuf, int width, int height, int myid, int numprocs,
         double xmin, double xmax, double ymin, double ymax, int maxiter ) {
  double dx,dy,x,y;
  int    ix,iy;

  dx = (xmax - xmin) / width;
  dy = (ymax - ymin) / height;

  y = ymin + dy * myid * height/numprocs;
  for (iy=0; iy<height/numprocs; ++iy) {
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
  
  MPI_Gather(iterations, height/numprocs * width, MPI_INT, recvbuf, height/numprocs * width, MPI_INT, 0, MPI_COMM_WORLD);
}

void calcMaster(int *iterations, int *recvbuf, int width, int height, int myid, int numprocs,
         double xmin, double xmax, double ymin, double ymax, int maxiter ) {
  double dx,dy,x,y;
  int    ix,iy;
  int    inwork = 0;
  int    *buffer;
  int    curr;
  int    cancel = - 1;
  int    todo = (height * width) / (MASTER_BLOCKSIZE * MASTER_BLOCKSIZE);
  MPI_Status status;

  printf("TEST (procs = %d)\n", numprocs);
  buffer = malloc(sizeof(int) * (MASTER_BLOCKSIZE * MASTER_BLOCKSIZE + 1));

  dx = (xmax - xmin) / width;
  dy = (ymax - ymin) / height;

  y = ymin;
  curr = 0;

  if (myid == 0) {

    for (int i = 1; i < numprocs; i++) {
      printf("MASTER: Initial work spreading\n");
      MPI_Send(&curr, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      curr++;
      inwork++;
    } 

    while(inwork != 0) {
      //printf("MASTER: waiting for more jobs. (cur=%d)\n", curr);
      MPI_Recv(buffer, MASTER_BLOCKSIZE * MASTER_BLOCKSIZE + 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
      inwork--;

      y = (buffer[0] * MASTER_BLOCKSIZE/width) * dy + ymin;
      for (iy=0; iy<MASTER_BLOCKSIZE; ++iy) {
	x = ((buffer[0] * height) % MASTER_BLOCKSIZE) * dy + xmin;
	for (ix=0; ix<MASTER_BLOCKSIZE; ix++) {
	  recvbuf[iy*width+ix] = buffer[iy*MASTER_BLOCKSIZE+ix+1];
	  x += dx;
	}
	y += dy;
      }


      if (curr < todo) {
	MPI_Send(&curr, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
	inwork++;
	curr++;
      } else {
	MPI_Send(&cancel, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
      }
    }
  } else {
    while(1) {
      MPI_Recv(&curr, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

      if (curr == cancel) {
	printf("Worker canceld\n");
	break;
      }
      
      printf("Next Block: x = %d, y = %d\n", ((curr * height) % MASTER_BLOCKSIZE) * dx + xmin, (curr * width / MASTER_BLOCKSIZE) * dy + ymin);

      y = (curr * width / MASTER_BLOCKSIZE) * dy + ymin; 
       for (iy=0; iy<MASTER_BLOCKSIZE; ++iy) {
	 x = ((curr * height) % MASTER_BLOCKSIZE) * dx + xmin ;
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
	   buffer[iy*MASTER_BLOCKSIZE+ix + 1] = count;
	   x += dx;
	 }
	 y += dy;
       }

       buffer[0] = curr;

       MPI_Send(buffer, MASTER_BLOCKSIZE*MASTER_BLOCKSIZE + 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
  }
}
