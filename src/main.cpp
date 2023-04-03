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
    int np, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&np);//npに、MPIで並列化するサイズが入る
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);//myrankに、MPIで並列化したときの自身のスレッドが入る
    omp_set_num_threads(2);

    if(myrank==0){//ここでは、0番目のスレッドのみ実行される
        cout << "this is main waiting 3 sec" << endl;
        double t1 = MPI_Wtime();
        double t2 = t1;
        while(t2-t1 < 3.0){
            t2 = MPI_Wtime();
        }
        cout << "after count" << endl;
    }

    //--ここですべてのプロセスが終わるのを待つ
    MPI_Barrier(MPI_COMM_WORLD);

    int i,thread_num;
    #pragma omp parallel for private(i)//,thread_num)//ここで指定した変数は、ループ内で独立してる
    for(i=0;i<3;i++){
        thread_num = omp_get_thread_num();
        //cout << "this myrank is : " << myrank << " , size = : " << np << " , thread = : " << thread_num << endl;

        string text = "this myrank is : " + to_string(myrank) + " , size = : " + to_string(np) + " , thread = : " + to_string(thread_num) + "\n";
        cout << text;
    }
    MPI_Finalize();

    return 0;
}
