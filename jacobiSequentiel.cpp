#include "jacobiSequentiel.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>

using namespace std;

// -------- CREATOR -------- //

JacobiSequentiel::JacobiSequentiel(int argc, char* argv[]):
_a(atof(argv[1])), _b(atof(argv[2])),
_u0(atof(argv[3])), _alpha(atof(argv[4])),
_nx(atoi(argv[5])+1), _ny(atoi(argv[6])+1)
{
    MPI_Comm_size(MPI_COMM_WORLD, &_nbTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &_myRank);

    _grid = new Grid(0,0,_a,_b,_nx-2,_ny-2,_u0);
    _dx = _a/(_nx-1);
    _dy = _b/(_ny-1);

    // Limit condition
    if (MODE == 1)
    {
        initLimitCondition1();
    }
    else if (MODE == 2)
    {
        initLimitCondition2();
    }
}

// -------- DESTRUCTOR -------- //

JacobiSequentiel::~JacobiSequentiel()
{
    delete _grid;
    delete _left;
    delete _right;
    delete _up;
    delete _down;
}

// -------- ACCESSOR -------- //

int JacobiSequentiel::getRank() const
{
    return _myRank;
}

// -------- OTHER -------- // 

double JacobiSequentiel::resolve()
{
    if (_myRank == 0)
    {
        double timeInit = MPI_Wtime();
        int k = 0;
        double diff;
        do
        {
            k += 1;

            // Compute next step
            diff = jacobi();

        }while (k<=MAX_STEP && diff>=MAX_DIFF);
        double timeEnd = MPI_Wtime();
        return timeEnd-timeInit;
    }
    return 0;
}

void JacobiSequentiel::saveData() const
{
    if (_myRank == 0)
    {
        ofstream file;
        file.open("jacobiSeq.txt", ios::out);
        file << *_grid;
        file.close();
    }
}

double JacobiSequentiel::jacobi()
{
    double diff = 0; // Difference between the new and the old grid
    Grid *newGrid;
    newGrid = new Grid(0,0,_a,_b,_nx-2,_ny-2);
    double constant = _dx*_dx*_dy*_dy/(2.*_dx*_dx+2.*_dy*_dy);
    for (int i = 0;i<_nx-2;++i)
    {
        for (int j = 0;j<_ny-2;++j)
        {
            if (i>0)
            {
                newGrid->get(i,j) += constant*_grid->get(i-1,j)/(_dx*_dx);
            }
            else
            {
                newGrid->get(i,j) += constant*_left[j]/(_dx*_dx);
            }
            if (j>0)
            {
                newGrid->get(i,j) += constant*_grid->get(i,j-1)/(_dy*_dy);
            }
            else
            {
                newGrid->get(i,j) += constant*_down[i]/(_dy*_dy);
            }
            if (i+1<_nx-2)
            {
                newGrid->get(i,j) += constant*_grid->get(i+1,j)/(_dx*_dx);
            }
            else
            {
                newGrid->get(i,j) += constant*_right[j]/(_dx*_dx);
            }
            if (j+1<_ny-2)
            {
                newGrid->get(i,j) += constant*_grid->get(i,j+1)/(_dy*_dy);
            }
            else
            {
                newGrid->get(i,j) += constant*_up[i]/(_dy*_dy);
            }
            if (MODE == 1)
            {
                newGrid->get(i,j) -= constant*f1(i*_dx+_dx,j*_dy+_dy);
            }
            else if (MODE == 2)
            {
                newGrid->get(i,j) -= constant*f2(i*_dx+_dx,j*_dy+_dy);
            }
            diff += (newGrid->get(i,j)-_grid->get(i,j))*(newGrid->get(i,j)-_grid->get(i,j));
        }
    }
    delete _grid;
    _grid = newGrid;
    return sqrt(diff)/((_nx-2)*(_ny-2));
}

void JacobiSequentiel::initLimitCondition1()
{
    // Limit condition
    _left = new double[_ny];
    _right = new double[_ny];
    for (int i=0;i<_ny;++i)
    {
        _left[i]=_u0*(1+_alpha*(1+cos(2*M_PI*((i+1)*_dy-_b/2.)/_b)));
        _right[i]=_u0;
    }
    _up = new double[_nx];
    _down = new double[_nx];
    for (int i=0;i<_nx;++i)
    {
        _up[i]=_u0;
        _down[i]=_u0;
    }
}

void JacobiSequentiel::initLimitCondition2()
{
    // Limit condition
    _left = new double[_ny];
    _right = new double[_ny];
    for (int i=0;i<_ny;++i)
    {
        _left[i]=sin(M_PI*(1-(i+1)*_dy/_b));
        _right[i]=sin(M_PI*(1-(i+1)*_dy/_b));
    }
    _up = new double[_nx];
    _down = new double[_nx];
    for (int i=0;i<_nx;++i)
    {
        _up[i]=0;
        _down[i]=0;
    }
}

double JacobiSequentiel::f1(const double &x, const double &y) const
{
    return 0;
}

double JacobiSequentiel::f2(const double &x, const double &y) const
{
    return -M_PI*M_PI/(_b*_b)*sin(M_PI*(1-y/_b));
}

double JacobiSequentiel::error() const
{
    double err(0);
    if (MODE==2)
    {
        for (int i=0;i<_nx-2;++i)
        {
            for (int j=0;j<_ny-2;++j)
            {
                err+=(_grid->get(i,j)-sin(M_PI*(1-(j+1)*_dy/_b)))*(_grid->get(i,j)-sin(M_PI*(1-(j+1)*_dy/_b)));
            }
        }
    }
    return err/((_nx-2)*(_ny-2));
}