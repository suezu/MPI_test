#include "mpi.h"
#include <string>
#include <iostream>
using std::string;
using std::cout;

string printArray(int* array,int N){
    string text = "";
    text += "[ ";
    for(int i=0;i<N;i++){
        text += std::to_string(array[i]) + " ,";
    }
    text += "]\n";
    return text;
}

string printArray(int** array,int N,int M){
    string text = "";
    text += "[ ";
    for(int i=0;i<N;i++){
        text += "[ ";
        for(int j=0;j<M;j++){
            text += std::to_string(array[i][j]) + " ,";
        }
        text += "]\n";
    }
    text +=  "]\n";
    return text;
}

int main(int argc, char* argv[]) {
    const int N=4;
    const int M = 4;

    int  A[M][N] = {{1, 2, 3, 4},
                    {5, 6, 7, 8},
                    {9, 10, 11 ,12},
                    {13, 14, 15, 16}};
    int local_A[N][M];
    int  X[N] = {1, 2, 3, 4}, local_X[N];
    int  Y[N], local_Y[N];
    int local_m, local_n;
    int i,j;

    int p;//ランク数
    int myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    local_m = M/p;
    local_n = N/p;

    //行列Aをlocal_Aに分散
    MPI_Scatter(A, local_m*N, MPI_INT, local_A,local_m*N, MPI_INT, 0, MPI_COMM_WORLD);

    //ベクトルXをloacl_Xに分散
    MPI_Scatter(X, local_n, MPI_INT, local_X, local_n, MPI_INT,0, MPI_COMM_WORLD);
    

    MPI_Allgather(local_X, local_n, MPI_INT, X, local_n, MPI_INT,MPI_COMM_WORLD);
    for (i = 0; i < local_m; i++) {
        local_Y[i] = 0;
        for (j = 0; j <  N; j++)
        local_Y[i] = local_Y[i] + local_A[i][j]*X[j];
    }

    MPI_Gather(local_Y, local_m, MPI_INT, Y, local_m, MPI_INT, 0,MPI_COMM_WORLD);

    if (myrank == 0){
        for (i=0; i < M; i++) printf("%d ", Y[i]);
        printf("\n");
    }

    MPI_Finalize();

    return 0;
}