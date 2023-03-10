/**********************************************************************************
 scatterv.c ejemplo del uso de la funcion scatterv
 Reparto de un vector de 8 elementos entre 4 proces.: 2, 4, 1 y 1 elementos
 Solo para el caso particular de 4 procesos
********************************************************************************/
#include <stdio.h>
#include <mpi.h>
int main(int argc, char** argv)
{
 int pid, npr, i;
 int A[8], B[8], tam[4], desp[4];
 MPI_Init(&argc, &argv);
 MPI_Comm_rank(MPI_COMM_WORLD, &pid);
 MPI_Comm_size(MPI_COMM_WORLD, &npr);
/* Inicializaciones en todos los procesos */
 for(i=0; i<8; i++) { A[i] = -1; B[i] = -1; }
 if (pid == 0) {
 for(i=0; i<8; i++) A[i] = i;
 printf("\n A en P0 = %2d %2d %2d %2d %2d %2d %2d %2d\n",A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7]);
 }
/* Reparto de A (8 elementos) desde P0 en 4 trozos: 2, 4, 1, 1 elementos */
/* vector de tamainos de los trozos */
 tam[0] = 2; tam[1] = 4; tam[2] = tam[3] = 1;
/* vector de desplazamientos desde el origen a cada trozo */
 desp[0] = 0; desp[1] = 2; desp[2] = 6; desp[3] = 7;
 MPI_Scatterv(&A[0], tam, desp, MPI_INT, &B[0], tam[pid], MPI_INT, 0, MPI_COMM_WORLD);
 printf("\n SCATTERV: P%d B =%2d %2d %2d %2d %2d %2d %2d %2d\n",pid,B[0],B[1], B[2],B[3],B[4],B[5],B[6],B[7]);
 MPI_Finalize();
 return 0;
}
