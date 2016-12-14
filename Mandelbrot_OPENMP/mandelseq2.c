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
  fprintf(stderr, "  [-t <type>]               0=stride, 1=stripe, 2=block, 3=blockmaster\n");
  fprintf(stderr, "  [-v]                      verbose (off)\n\n");
  exit(1);
}

void calc(int *iterations, int width, int height, int myid, int numprocs,
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
  int *iterations;
  int    ix,iy;

  int    numprocs,myid;
  int    i;
  double st,timeused,calctime=0.0, waittime=0.0, iotime=0.0;
  char   filename[1024];

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

  /* initialize arrays */
  iterations = malloc(width*height*sizeof(int));

  for (ix=0; ix<width; ++ix) {
    for (iy=0; iy<height; ++iy) {
      iterations[ix*height+iy] = 0;
    }
  }

  numprocs = 1;
  myid     = 0;


  /* start calculation */
  if(verbose) {
    printf("start calculation (x=%8.5g ..%8.5g,y=%10.7g ..%10.7g) ... \n",
           xmin,xmax,ymin,ymax);
    fflush(stdout);
  }

  st = esecond();
  calc(iterations, width, height, myid, numprocs, xmin, xmax, ymin, ymax, maxiter );
  timeused = esecond()-st;
  calctime += timeused;


  st = esecond();
  ppmwrite(iterations,width,height,0,maxiter,"mandelcol.ppm");

  timeused = esecond()-st;
  iotime += timeused;
  if(verbose) printf("PE%02d: calc=%7.4f,wait=%7.4f, io=%7.4f\n",
                     myid,calctime,waittime,iotime);

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
