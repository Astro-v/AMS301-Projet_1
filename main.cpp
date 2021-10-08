/*
PROJET 1 - AMS301
MICHEL Valentin
Probleme sur grille structure

mpicxx *.cpp

mpirun -np 4 ./a.out 'a' 'b' 'U0' 'alpha' 'Nx' 'Ny'
mpirun -np 4 ./a.out 1 1 1 0.5 11 11

*/

// LIBRARY
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>


// MY FILE
#include "constant.hpp"
#include "grid.hpp"
#include "jacobiSequentiel.hpp"
#include "jacobiParallelise.hpp"
#include "gaussSeidelSequentiel.hpp"

using namespace std;

int main(int argc, char* argv[])
{
    int myRank, nbTasks;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nbTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    /*JacobiSequentiel seq(argc,argv);
    double timeSeq = seq.resolve();
    seq.saveData();
    if (myRank==0)
    {
        cout << "Temps séquentiel : " << timeSeq << endl;
    }*/

    MPI_Barrier(MPI_COMM_WORLD);
            
    JacobiParallelise par(argc,argv);
    double timePar = par.resolve();
    if (myRank==0)
    {
        cout << "Temps pour jacobi parallelise pour " << nbTasks << " coeur : " << timePar << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
            
    GaussSeidelSequentiel seq(argc,argv);
    double timeSeq = seq.resolve();
    if (myRank==0)
    {
        cout << "Temps pour Gauss-Seidel sequentiel pour " << nbTasks << " coeur : " << timeSeq << endl;
    }
    

    // Finalize MPI
    MPI_Finalize();

    return 0;
}

