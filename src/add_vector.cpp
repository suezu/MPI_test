#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <omp.h>
#include "mpi.h"
using namespace std;


//実行ではmpirunですること
int main(int argc, char **argv)
{   
    bool ISdebug = true;
    int n = 16; // ベクトルの要素数
    int *a, *b; // ベクトルaとb
    int partial_sum = 0.0; // 各プロセスの部分和
    int global_sum = 0.0; // 全プロセスの和


    int np_size, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&np_size);//npに、MPIで並列化するサイズが入る
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);//myrankに、MPIで並列化したときの自身のスレッドが入る
    //omp_set_num_threads(2);

    //ベクトルの初期化
    if(myrank==0){//ここでは、0番目のスレッドのみ実行される
        cout << "proogram run" << endl;
        a = new int[n];
        b = new int[n];

        for(int i=0;i<n;i++){
            a[i] = i;
            b[i] = 2*i;
        }

        //初期値の出力
        if(ISdebug){
            cout << "a : [ ";
            for(int i=0;i<n;i++){
                cout << a[i] << ", ";
            }
            cout << "]" << endl << "b : [ ";
            for(int i=0;i<n;i++){
                cout << b[i] << ", ";
            }
            cout << "]" << endl;
        }
    }

    //--ここですべてのプロセスが終わるのを待つ
    MPI_Barrier(MPI_COMM_WORLD);

    //各プロセスにベクトルの割り当て
    int partialSize = n/np_size;
    int *local_a = new int[partialSize];
    int *local_b = new int[partialSize];
    
    MPI_Scatter(a, partialSize, MPI_INT, local_a, partialSize, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(b, partialSize, MPI_INT, local_b, partialSize, MPI_INT, 0, MPI_COMM_WORLD);

    //各プロセスに割り当てられた配列
    if(ISdebug==true){
        string text = "myrannk : " + to_string(myrank) + ", array [ ";
        for(int i=0;i<partialSize;i++){
            text += to_string(local_a[i]) + ", ";
        }
        text += "]\n";
        cout << text;
    }

    //各プロセスに割り当てられた内積を計算
    for(int i=0;i<partialSize;i++){
        partial_sum += local_a[i] * local_b[i];
    }
    
    //プロセス0に集約
    MPI_Reduce(&partial_sum, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    //プロセス0が和を出力する
    if(myrank==0){
        cout << "Dot = " << global_sum << endl;
    }
    
    MPI_Finalize();

    return 0;
}
