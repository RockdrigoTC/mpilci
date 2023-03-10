#include <mpi.h>
#include <stdio.h>
#define N 10
int main (int argc, char **argv)
{
 int pid, npr, origen, destino, ndat, tag;
 int VA[N], i;
 MPI_Status info;
 MPI_Init(&argc, &argv);
 MPI_Comm_rank(MPI_COMM_WORLD, &pid);
 MPI_Comm_size(MPI_COMM_WORLD, &npr);
 for (i=0; i<N; i++) VA[i] = 0;
 if (pid == 0)
 {
 for (i=0; i<N; i++) VA[i] = i;
 destino = 1; tag = 0;
 MPI_Send(&VA[0], N, MPI_INT, destino, tag, MPI_COMM_WORLD);
 }
 else if (pid == 1)
 {
 printf("\nvalor de VA en P1 antes de recibir datos\n\n");
 for (i=0; i<N; i++) printf("%4d", VA[i]);
 printf("\n\n");
 origen = 0; tag = 0;
 MPI_Recv(&VA[0], N, MPI_INT, origen, tag, MPI_COMM_WORLD, &info);
 MPI_Get_count (&info, MPI_INT, &ndat);
 printf("Pr 1 recibe VA de Pr %d: tag %d, ndat %d \n\n",
 info.MPI_SOURCE, info.MPI_TAG, ndat);
 for (i=0; i<ndat; i++) printf("%4d", VA[i]);
 printf("\n\n");
 }
 MPI_Finalize();
 return 0;
} 
