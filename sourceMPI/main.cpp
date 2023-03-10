#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

using namespace std;
#define NE 100000

int* arrayDesorden;
int* arrayOrdS;
int* arrayOrdP;
int* arrayBloque = NULL;

void leeDatos(const char* nombre, int* v) 
{
    FILE* archivo;
    fopen_s(&archivo, nombre, "r");
    int* buf = (int*)malloc(NE * sizeof(int));

    for (int i = 0; i < NE; i++) {
        fscanf_s(archivo, "%ld", &buf[i]);
        v[i] = buf[i];
    }
    fclose(archivo);
}

void imprimeDatos(const char* nombre, int* a) 
{
    FILE* archivo;
    fopen_s(&archivo, nombre, "w");
    for (int i = 0; i < NE; i++) {
        fprintf(archivo, "%7ld", a[i]);
        fprintf(archivo, "\n");
    }
    fclose(archivo);
}

void bubbleSort(int* array, int n) 
{
	int i, j;                
	int temp = 0;             

	for (i = 1;i < n;i++)
	{
		for (j = 0; j < n - i;j++)
		{
			if (array[j] > array[j + 1])
			{
				temp = array[j];
				array[j] = array[j + 1];
				array[j + 1] = temp;
			}
		}
	}
}

int buscaVM(int* array)
{
    int maximo = array[0];
    for (int i = 0; i < NE; i++)
    {
        if (maximo < array[i]) {
            maximo = array[i];
        }
    }
    return maximo;
}

int comparaInt(const void* a, const void* b) 
{

    int int_a = *((int*)a);
    int int_b = *((int*)b);

    if (int_a == int_b) return 0;
    else if (int_a < int_b) return -1;
    else return 1;

}

void ensamblaBloque(int* origen, int* destino, int inicio, int fin) 
{
    int k = 0;
    for (int i = inicio; i < fin; i++)
    {
        destino[i] = origen[k];
        k++;
    }
}

double coeficienteCorrelacion(int* a, int* b, int n)
{
	double sumaA = 0, sumaB = 0, sumaA2 = 0, sumaB2 = 0, sumaAB = 0;

	for (int i = 0; i < n; i++)
	{
		sumaA += a[i];
		sumaB += b[i];
		sumaAB += a[i] * b[i];
		sumaA2 += a[i] * a[i];
		sumaB2 += b[i] * b[i];

	}

	for (int i = 0; i < n; i++)
	{
		if (abs(a[i] - b[i]) > 10)
		{
			printf("Son diferentes %i %i %i\n", i, a[i],b[i]);
			system("pause");
		}
	}

	return (n * sumaAB - sumaA * sumaB) / ((sqrt(n * sumaA2 - sumaA * sumaA)) * (sqrt(n * sumaB2 - sumaB * sumaB)));
}

int main(int argc, char** argv) 
{
    const char* archivo1 = "numerosDesordenados.dat"; 
    const char* archivo2 = "numerosOrdenados.dat"; 
    int numEB, rango, i, nprocs, id, posInicial, posFin;
    double TS, TP, start, end;
    MPI_Status status;
    MPI_Request request;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    int NB = nprocs - 1; 
    vector<int>* bloques = new vector<int>[NB];

    if (id == 0) 
    {
        arrayDesorden = (int*)malloc(NE * sizeof(int));
        arrayOrdP = (int*)malloc(NE * sizeof(int));
        arrayOrdS = (int*)malloc(NE * sizeof(int));
       
        leeDatos(archivo1,arrayDesorden);
        arrayOrdP = arrayDesorden;
        arrayOrdS = arrayDesorden;

        start = MPI_Wtime(); 
        printf(" \n -Inicia secuencial- \n"); fflush(stdout);
        //qsort(arrayOrdS, NE, sizeof(int), comparaInt); 
        bubbleSort(arrayOrdS, NE);                         
        end = MPI_Wtime(); 
        printf(" -Termina secuencial- \n"); fflush(stdout);
        TS = end - start;

        rango = buscaVM(arrayDesorden) / NB;
        if (rango == 0) rango = 1; 
        for (i = 0; i < NE; i++)
        {
            int dato = arrayDesorden[i];
            int idx = dato / rango;
            if (idx > NB - 1)
            {
                idx -= idx - (NB - 1);
            }
            bloques[idx].push_back(dato);
        }
        
        printf(" \n -Inicia MPI- \n"); fflush(stdout);
        start = MPI_Wtime();
        int destino = 0;
        for (i = 0; i < NB; i++)
        {
            numEB = bloques[i].size();
            arrayBloque = (int*)malloc(numEB * sizeof(int));
            copy(bloques[i].begin(), bloques[i].end(), arrayBloque);
            destino++;
            MPI_Isend(&numEB, 1, MPI_INT, destino, 1, MPI_COMM_WORLD, &request);
            printf("Nodo 0 envia bloque desordenado al nodo %d \n", destino); fflush(stdout);
            MPI_Isend(arrayBloque, numEB, MPI_INT, destino, 1, MPI_COMM_WORLD, &request);
        }

        posInicial = 0, posFin = 0, destino = 0;
        for (i = 0; i < NB; i++)
        {
            destino++;
            MPI_Recv(&numEB, 1, MPI_INT, destino, 1, MPI_COMM_WORLD, &status);
            arrayBloque = (int*)malloc(numEB * sizeof(int));
            printf("Nodo 0 recibe bloque ordenado del nodo %d\n", destino); fflush(stdout);
            MPI_Recv(arrayBloque, numEB, MPI_INT, destino, 1, MPI_COMM_WORLD, &status);

            posInicial = posFin;
            posFin += numEB;
            ensamblaBloque(arrayBloque, arrayOrdP, posInicial, posFin);
        }
        printf(" -Termina MPI- \n"); fflush(stdout);
        end = MPI_Wtime();
        TP = end - start;
        printf(" \n Tiempo Secuencial = %.4lf Tiempo Paralelo = %.4lf Speedup = %.4lf\n", TS, TP, (TS / TP));
        printf(" Coeficiente de correlacion %lf\n", coeficienteCorrelacion(arrayOrdP, arrayOrdS, NE));
        imprimeDatos(archivo2, arrayOrdP);
    }
    else {
        MPI_Recv(&numEB, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        arrayBloque = (int*)malloc(numEB * sizeof(int));
        printf("Nodo %d recibe bloque desordenado del nodo 0\n", id); fflush(stdout);
        MPI_Recv(arrayBloque, numEB, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        //qsort(arrayBloque, numEB, sizeof(int), comparaInt); 
        bubbleSort(arrayBloque, numEB);                        
        MPI_Isend(&numEB, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &request);
        printf("Nodo %d envia bloque ordenado al nodo 0 \n", id); fflush(stdout);
        MPI_Isend(arrayBloque, numEB, MPI_INT, 0, 1, MPI_COMM_WORLD, &request);
    }
    MPI_Finalize();
}