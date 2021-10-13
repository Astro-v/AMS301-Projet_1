#include "gaussSeidelSequentiel.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>

using namespace std;

// -------- CREATOR -------- //

GaussSeidelSequentiel::GaussSeidelSequentiel(int argc, char* argv[]):
_a(atof(argv[1])), _b(atof(argv[2])),
_u0(atof(argv[3])), _alpha(atof(argv[4])),
_nx(atoi(argv[5])+1), _ny(atoi(argv[6])+1)
{
    MPI_Comm_size(MPI_COMM_WORLD, &_nbTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &_myRank);

    _grid = new Grid(0,0,_a,_b,_nx-2,_ny-2);
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

GaussSeidelSequentiel::~GaussSeidelSequentiel()
{
    delete _grid;
    delete _left;
    delete _right;
    delete _up;
    delete _down;
}

// -------- ACCESSOR -------- //

int GaussSeidelSequentiel::getRank() const
{
    return _myRank;
}

// -------- OTHER -------- // 

double GaussSeidelSequentiel::resolve()
{
    if (_myRank == 0)
    {
        double timeInit = MPI_Wtime();
        int k = 0;
        do
        {
            k += 1;

            // Compute next step for red
            gaussSeidel(0);

            // Compute next step for black
            gaussSeidel(1);

        }while (k<=MAX_STEP);
        double timeEnd = MPI_Wtime();
        return timeEnd-timeInit;
    }
    return 0;
}

void GaussSeidelSequentiel::saveData() const
{
    if (_myRank == 0)
    {
        ofstream file;
        file.open("gaussSeidelSeq.txt", ios::out);
        file << *_grid;
        file.close();
    }
}

double GaussSeidelSequentiel::gaussSeidel(int step)
{
    double constant = _dx*_dx*_dy*_dy/(2.*_dx*_dx+2.*_dy*_dy);
    for (int i = 0;i<_nx-2;++i)
    {
        for (int j = 0;j<_ny-2;++j)
        {
            if ((i+j)%2==step) // if step = 0 we compute on even and odd if step = 1
            {
                _grid->get(i,j) = 0;
                if (i>0)
                {
                    _grid->get(i,j) += constant*_grid->get(i-1,j)/(_dx*_dx);
                }
                else
                {
                    _grid->get(i,j) += constant*_left[j]/(_dx*_dx);
                }
                if (j>0)
                {
                    _grid->get(i,j) += constant*_grid->get(i,j-1)/(_dy*_dy);
                }
                else
                {
                    _grid->get(i,j) += constant*_down[i]/(_dy*_dy);
                }
                if (i+1<_nx-2)
                {
                    _grid->get(i,j) += constant*_grid->get(i+1,j)/(_dx*_dx);
                }
                else
                {
                    _grid->get(i,j) += constant*_right[j]/(_dx*_dx);
                }
                if (j+1<_ny-2)
                {
                    _grid->get(i,j) += constant*_grid->get(i,j+1)/(_dy*_dy);
                }
                else
                {
                    _grid->get(i,j) += constant*_up[i]/(_dy*_dy);
                }
                if (MODE == 1)
                {
                    _grid->get(i,j) -= constant*f1(i*_dx+_dx,j*_dy+_dy);
                }
                else if (MODE == 2)
                {
                    _grid->get(i,j) -= constant*f2(i*_dx+_dx,j*_dy+_dy);
                }
            }
        }
    }
}

void GaussSeidelSequentiel::initLimitCondition1()
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

void GaussSeidelSequentiel::initLimitCondition2()
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

double GaussSeidelSequentiel::f1(const double &x, const double &y) const
{
    return 0;
}

double GaussSeidelSequentiel::f2(const double &x, const double &y) const
{
    return -M_PI*M_PI/(_b*_b)*sin(M_PI*(1-y/_b));
}