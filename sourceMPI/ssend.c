/**********************************************************************************
 ssend.c
 Ping-pong entre dos procesadores
 Comunicacion sincrona ssend
**********************************************************************************/
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <unistd.h>
#define VUELTAS 4
double calculo()
{
 double aux;

 sleep(1);
 aux = rand() % 100;

 return(aux);
}
int main(int argc, char** argv)
{
 double t0, t1, dat= 0.0, dat1, dat_rec;
 int pid, i;
 MPI_Status status;
 MPI_Init(&argc, &argv);
 MPI_Comm_rank(MPI_COMM_WORLD, &pid);
 if (pid == 0) t0 = MPI_Wtime();
 for(i=0; i<VUELTAS; i++)
 {
 if (pid == 0)
 {
 MPI_Ssend(&dat, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
 dat1 = calculo();
 MPI_Recv(&dat_rec, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status);
 dat = dat1 + dat_rec;
 }
 else
 {
 dat1 = calculo();
 MPI_Recv(&dat_rec, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
 dat = dat1 + dat_rec;
 MPI_Ssend(&dat, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
 }
 }
 if (pid == 0)
 {
 t1 = MPI_Wtime();
 printf("\n Tiempo de ejecucion = %f s \n", t1-t0);
 printf("\n Dat = %1.3f \n\n", dat);
 }
 MPI_Finalize();
 return 0;
}
