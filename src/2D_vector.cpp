#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
using std::cout;
using std::endl;
using std::string;

void conv2to1(int dim0,int dim1,int **array2D,int *array1D){
    //2次元はいれつを[[dim1][dim1]....]に変換
    for(int i=0;i<dim0;i++){
        for(int j=0;j<dim1;j++){
            array1D[i*dim1+j] = array2D[i][j];
        }
    }
}

void conv1to2(int dim0,int dim1,int *array1D,int **array2D){
    //1次元配列[[dim1][dim1]]を[[dim1],[dim1],...に変換]
    for(int i=0;i<dim0;i++){
        for(int j=0;j<dim1;j++){
            array2D[i][j] = array1D[i*dim1+j];
        }
    }
}

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

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    int N = 12;
    int M = 5;
    int **A,**local_A;
    int *A_1D,*local_A_1D;
    int *B;
    int taskSize = N/world_size;
    
    /*
    **配列の確保
    */
    A = new int*[N];
    B = new int[N];
    local_A = new int*[taskSize];
    for(int i=0; i<N; i++) {
        A[i] = new int[M];
    }
    for(int i=0;i<taskSize;i++){
        local_A[i] = new int[M];
    }

    local_A_1D = new int[taskSize*M];
    A_1D = new int[N*M];

    if(myrank==0){
        if(N%world_size!=0){
            exit(0);
        }

        // 以下のコードでN行m列の配列を生成します
        // 配列をランダムに初期化します
        for(int i=0; i<N; i++) {
            B[i] = i;
            for(int j=0; j<M; j++) {
                //A[i][j] = rand() % 100;
                A[i][j] = 10*i + j;
            }
        }

        conv2to1(N,M,A,A_1D);

        cout << "init array A : " << endl;
        cout << printArray(A,N,M);

        //cout << "init array B : " << endl;
        //cout << printArray(B,N);
    }

    //--ここですべてのプロセスが終わるのを待つ
    MPI_Barrier(MPI_COMM_WORLD);

    //Bcastで配列を送る
    MPI_Bcast(&B[0],N,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&A_1D[0],N*M,MPI_INT,0,MPI_COMM_WORLD);
    conv1to2(N,M,A_1D,A);

    for(int i=0;i<taskSize;i++){
        for(int j=0;j<M;j++){
            int ind = myrank*taskSize + i;
            local_A[i][j] = A[ind][j] + 1000*myrank;
        }
    }
    conv2to1(taskSize,M,local_A,local_A_1D);

    //--ここですべてのプロセスが終わるのを待つ
    MPI_Barrier(MPI_COMM_WORLD);

    //allgatherで送信
    MPI_Allgather(&local_A_1D[0],taskSize*M,MPI_INT,&A_1D[0],taskSize*M,MPI_INT,MPI_COMM_WORLD);
    conv1to2(N,M,A_1D,A);


    cout << "rank : " + std::to_string(myrank) + " : \n" + printArray(A,N,M);

    // メモリを解放する
    for(int i=0; i<N; i++) {
        delete[] A[i];
    }
    delete[] A;
    delete[] B;
    MPI_Finalize();
    return 0;
}
