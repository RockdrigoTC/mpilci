#include <mpi.h>
#include <stdio.h>
#define DECBIN(n,i) ((n&(1<<i))?1:0)
void test (int pid, int z)
{
 int v[16], i;
 for (i=0; i<16; i++) v[i] = DECBIN(z,i);
 if ((v[0] || v[1]) && (!v[1] || !v[3]) && (v[2] || v[3])
 && (!v[3] || !v[4]) && (v[4] || !v[5])
 && (v[5] || !v[6]) && (v[5] || v[6])
 && (v[6] || !v[15]) && (v[7] || !v[8])
 && (!v[7] || !v[13]) && (v[8] || v[9])
 && (v[8] || !v[9]) && (!v[9] || !v[10])
 && (v[9] || v[11]) && (v[10] || v[11])
 && (v[12] || v[13]) && (v[13] || !v[14])
 && (v[14] || v[15]))
 {
 printf(" %d) %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d (%d)\n", pid, v[15],v[14],v[13],
 v[12],v[11],v[10],v[9],v[8],v[7],v[6],v[5],v[4],v[3],v[2],v[1],v[0], z);
 fflush(stdout);
 }
}
int main (int argc, char *argv[])
{
 int i, pid, npr;
 MPI_Init(&argc, &argv);
 MPI_Comm_size(MPI_COMM_WORLD, &npr);
 MPI_Comm_rank(MPI_COMM_WORLD, &pid);
 for (i=pid; i<65536; i += npr) test(pid, i);
 MPI_Finalize();
 return 0;
} 
