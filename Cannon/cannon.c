/******************************************************************************/
/*  MPI-template for the parallel calculation of C = C + A*B based on         */
/*  Cannon's algorithm                                                        */
/******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define N 16 /* matrix dimension (global) */

#define TAG_A 100
#define TAG_B 101
#define TAG_C 102


/******************************************************************************/
/*  Function MatrixMultiply()                                                 */
/*   - calculates C = C + A*B                                                 */
/******************************************************************************/
void MatrixMultiply(int n, double *A, double *B, double *C){

  int i, j, k;

  for (i=0; i<n; i++){
    for(j=0; j<n; j++)
      for(k=0; k<n; k++)
	C[i*n+j] += A[i*n+k] * B[k*n+j];
  }
}


/******************************************************************************/
/*  Function CannonMatrixMultiply                                             */
/*   - multiplies two quadratic matrices based on Cannon's algorithm          */
/******************************************************************************/
void CannonMatrixMultiply(int n, double *A, double *B, double *C, int blocksize,
                          int num_blocks, int mycoords[2], MPI_Comm comm_2d){

  int i;
  int right, left, up, down;
  int shiftsource, shiftdest;
	int shiftcoords[2];
	int ndim[2];
  MPI_Status status;

	MPI_Cartdim_get(comm_2d, ndim);
  /* TODO: compute ranks of all four neighbors */
	//up
	shiftcoords[0] = mycoords[0];
	shiftcoords[1] = (mycoords[1]-1+ndim[0])%ndim[0];
	MPI_Cart_rank(comm_2d, shiftcoords, &up);
	//down
	shiftcoords[0] = mycoords[0];
	shiftcoords[1] = (mycoords[1]+1)%ndim[0];
	MPI_Cart_rank(comm_2d, shiftcoords, &down);
	//left
	shiftcoords[0] = (mycoords[0]-1+ndim[0])%ndim[0];
	shiftcoords[1] = mycoords[1];
	MPI_Cart_rank(comm_2d, shiftcoords, &left);
	//right
	shiftcoords[0] = (mycoords[0]+1)%ndim[0];
	shiftcoords[1] = mycoords[1];
	MPI_Cart_rank(comm_2d, shiftcoords, &left);

  /* TODO: perform the initial matrix alignment for A and B */
	for (i = 0; i < mycoords[1]; i++) {
		MPI_Sendrecv_replace(A, blocksize, MPI_DOUBLE, left, 45, right, 45, comm_2d, &status);
	}
	for (i = 0; i < mycoords[0]; i++) {
		MPI_Sendrecv_replace(B, blocksize, MPI_DOUBLE, up, 45, down, 45, comm_2d, &status);
	}

  /* get into the main computation loop */
  for (i=0; i<num_blocks; ++i){

    /* TODO: compute C = C + A*B */
	  MatrixMultiply(n, A, B, C);

    /* TODO: shift matrix A left by one */
	  MPI_Sendrecv_replace(A, blocksize, MPI_DOUBLE, left, 45, right, 45, comm_2d, &status);

    /* TODO: shift matrix B up by one */
	  MPI_Sendrecv_replace(B, blocksize, MPI_DOUBLE, up, 45, down, 45, comm_2d, &status);
  }

  /* TODO: restore the original distribution of A and B */
	for (i = 0; i < mycoords[1]; i++) {
		MPI_Sendrecv_replace(A, blocksize, MPI_DOUBLE, right, 45, left, 45, comm_2d, &status);
	}
	for (i = 0; i < mycoords[0]; i++) {
		MPI_Sendrecv_replace(B, blocksize, MPI_DOUBLE, down, 45, up, 45, comm_2d, &status);
	}
}


/******************************************************************************/
/*                              M A I N                                       */
/******************************************************************************/

int main(int argc, char *argv[]) 
{

  int    i, j, k;
  double matA[N * N], matB[N * N], matC[N * N];  /* total matrices */
  double *locA, *locB, *locC;  /* total matrices */
  double errmax, vergl;

  int num_test_prints;
  int num_procs, my_rank;
  int blocksize, num_blocks, blockstart;
  int proc, count;
  int dims[2], periods[2];
  int my_rank_2d, mycoords[2];
	int blocklen;

  MPI_Status status;
  MPI_Datatype blockmat;
  MPI_Comm comm_2d;
  
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&num_procs);

  num_blocks = sqrt(num_procs);

  /* N divisible by sqrt(num_procs) ? */ 
  if ( (fabs(sqrt(num_procs) - num_blocks) > 0.01) || ((N % num_blocks) != 0) ){
    printf("N = %d must be divisible by sqrt(num_procs)!\n",N);
    exit(1);
  }

  blocksize = N / num_blocks;
	blocklen = sqrt(blocksize);

  /* TODO: ALL PROCESSES: allocate memory for local part of the matrices */
	locA = malloc(blocksize*sizeof(double));
	locB = malloc(blocksize*sizeof(double));
	locC = malloc(blocksize*sizeof(double));

  /* TODO: ALL PROCESSES: create the Cartesian topology, without rank reordering */
  dims[0] = blocksize;
  dims[1] = blocksize;
  periods[0] = 1;
  periods[1] = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_2d);
 
  /* TODO: get rank and coordinates with respect to the new communicator */
	MPI_Comm_rank(comm_2d, &my_rank_2d);
	MPI_Cart_coords(comm_2d, my_rank_2d, 2, mycoords);

  /* TODO: ALL PROCESSES: data type for the blockwise distribution of the matrices */
	MPI_Type_vector(blocklen, blocklen, num_blocks*blocklen, MPI_DOUBLE, &blockmat);
    MPI_Type_commit(&blockmat);

  /* MASTER: initialize matrices */
  if (my_rank == 0){
    k=1;
    for (i=0; i<N; i++) {
      for (j=0;j<N; j++) { 
	matA[i*N+j]= k++; 
	matB[i*N+j]=(i<=j)?1.0:0.0;
	matC[i*N+j]=0;
      }
    }
    /* test prints */
    num_test_prints = N < 10 ? N : 10;
    for (j=0;j<num_test_prints; j++) { /* j := Zeile */
      printf("A(2,%d)= %14.8f\n",j,matA[2*N+j]);
    }
    for (j=0;j<num_test_prints ; j++) { /* j := Zeile */
      printf("B(2,%d)= %14.8f\n",j,matB[2*N+j]);
    }

    /* initialize local matrix blocks */
    count = 0;
    for(i=0; i<blocksize; ++i){
      for (j=0; j<blocksize; ++j){
	locA[count] = matA[i*N + j];
	locB[count] = matB[i*N + j];
	count++;
      }
    }
	
	MPI_Request request;//TODO ?
    /* TODO: distribute matrices blockwise among processes*/
	int start;
	for (i = 0; i < num_procs; i++) {
		start = i * blocklen + (i/num_blocks) * blocklen*blocklen*blocklen;
		MPI_Isend(&matA[start], 1, blockmat, i, 42, comm_2d, &request);
	}
	for (i = 0; i < num_procs; i++) {
		start = i * blocklen + (i/num_blocks) * blocklen*blocklen*blocklen;
		MPI_Isend(&matB[start], 1, blockmat, i, 43, comm_2d, &request);
	}
  }

  /* TODO: WORKER: recieve matrix blocks */
  /*else{
	
  }*/
	//MPI_Status status;//???
	MPI_Recv(locA, 1, blockmat, 0, 42, comm_2d, &status);
	MPI_Recv(locB, 1, blockmat, 0, 43, comm_2d, &status);

  /* ALL PROCESSES: initialize matric C */
  for (i=0; i<blocksize*blocksize; ++i)
    locC[i] = 0;


  /*********************************************************/
  /* ALL PROCESSES:   call funktion CannonMatrixMultiply() */
  /*********************************************************/
  CannonMatrixMultiply(blocksize, locA, locB, locC, blocksize, num_blocks, mycoords, comm_2d);


  /* TODO: WORKER: send result to master */
  if(my_rank != 0){
	MPI_Send(locC, blocksize, MPI_DOUBLE, 0, 44, comm_2d);
  }

  /* MASTER: collect results and check correctness */
  else{

    /* TODO: collect results */
	  for (i = 0; i < blocklen; i++) {
		  for (j = 0; j < blocklen; j++) {
			 matC[j+num_blocks*i] = locC[i*N+j];
		  }
	  }
	  double* buffer = malloc(blocksize*sizeof(double));
	  for (i = 1; i < num_procs; i++) {
		  MPI_Recv(buffer, blocksize, MPI_DOUBLE, i, 44, comm_2d, &status);
		  for (j = 0; j < blocklen; j++) {
			  for (k = 0; k < blocklen; k++) {
				 matC[k+num_blocks*j] = buffer[blocklen*i+j*N+k];
			  }
		  }
	  }

    /* copy own results */
    count = 0;
    for(i=0; i<blocksize; ++i){
      for (j=0; j<blocksize; ++j){
	matC[i*N + j] = locC[count];
	count++;
      }
    }      
     
    /* check results */
    errmax=0.0;
    for (i=0; i<N; i++) {
      for (j=0;j<N; j++) { 
	vergl=(double)(i*(j+1.0)*N+(j+2.0)*(j+1.0)/2.0);
	if( fabs(matC[i*N+j]-vergl) > errmax)
	  {
	    errmax=fabs(matC[i*N+j]-vergl);
	    if(errmax>1.0)
	      {
		printf("C(%d,%d)=%14.8f, vergl=%14.8f \n",i,j,matC[i*N+j],vergl);
	      }
	  }
      }
    }
    for (j=0;j<num_test_prints; j++) { 
      printf("C(2,%d)= %14.8f, Vergleichswert=%d \n",j,matC[2*N+j],2*(j+1)*N+(j+2)*(j+1)/2);
    }
    printf("maximum error: %14.8f\n",errmax);
  }

  /* TODO: free communicator */

  MPI_Finalize();

  return 0;
}
