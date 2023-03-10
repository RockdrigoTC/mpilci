/**********************************************************************************
 BCsinc.c
 prueba de sincronicidad del BC
**********************************************************************************/
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#define N 100
int main(int argc, char** argv)
{
 int pid, A[N];
 MPI_Init(&argc, &argv);
 MPI_Comm_rank(MPI_COMM_WORLD, &pid);

 sleep(pid*2);
 printf("\n >> P%d: comienzo del BC \n", pid);
 MPI_Bcast(A, N, MPI_INT, 0, MPI_COMM_WORLD);
 printf("\n fin del BC en P%d \n", pid);
 MPI_Finalize();
 return 0;
} 

