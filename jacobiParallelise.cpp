#include "jacobiParallelise.hpp"
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <cmath>

using namespace std;

// -------- CREATOR -------- //

JacobiParallelise::JacobiParallelise(int argc, char* argv[]):
_a(atof(argv[1])), _b(atof(argv[2])),
_u0(atof(argv[3])), _alpha(atof(argv[4])),
_nx(atoi(argv[5])+1), _ny(atoi(argv[6])+1),
_dx(_a/(_nx-1)), _dy(_b/(_ny-1))
{
    // Initialize MPI
    MPI_Comm_size(MPI_COMM_WORLD, &_nbTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &_myRank);

    // New position for each process
    _x0 = (_a-_dx)*ceil((_nx-2)*double(_myRank)/_nbTasks+1)/double(_nx-2);
    _y0 = _dy;
    _xf = (_a-_dx)*ceil((_nx-2)*double(_myRank+1)/_nbTasks)/double(_nx-2);
    _yf = _b-_dy;
    _nxf = round(1+(_xf-_x0)*double(_nx-2)/(_a-_dx));
    _nyf = _ny-2;


    // Grid matrice
    _grid = new Grid(_x0,_y0,_xf,_yf,_nxf,_nyf,_u0);

    // Limit condition
    _left = new double[_nyf];
    _right = new double[_nyf];
    _up = new double[_nxf];
    _down = new double[_nxf];

    if (MODE == 1){init1();}
    else if (MODE == 2){init2();}
}

// -------- DESTRUCTOR -------- //

JacobiParallelise::~JacobiParallelise()
{
    delete _left;
    delete _right;
    delete _up;
    delete _down;
    delete _grid;
}

// -------- ACCESSOR -------- //

int JacobiParallelise::getRank() const
{
    return _myRank;
}

// -------- OTHER -------- // 

void JacobiParallelise::init1()
{
    for (int i=0;i<_nyf;++i)
    {
        if (_myRank+1 == _nbTasks)
        {
            _right[i] = _u0;
        }
        else 
        {
            _right[i] = 0;
        }
        if (_myRank == 0)
        {
            _left[i] = _u0*(1+_alpha*(1+cos(2*M_PI*((_y0+i*_dy)-_b/2.)/_b)));
        }
        else 
        {
            _left[i] = 0;
        }
    }
    for (int i=0;i<_nxf;++i)
    {
        _up[i] = _u0;
        _down[i] = _u0;
    }
}

void JacobiParallelise::init2()
{
    for (int i=0;i<_nyf;++i)
    {
        if (_myRank+1 == _nbTasks)
        {
            _right[i] = sin(M_PI*(1-(_y0+i*_dy)/_b));
        }
        else 
        {
            _right[i] = 0;
        }
        if (_myRank == 0)
        {
            _left[i] = sin(M_PI*(1-(_y0+i*_dy)/_b));
        }
        else 
        {
            _left[i] = 0;
        }
    }
    for (int i=0;i<_nxf;++i)
    {
        _up[i] = 0;
        _down[i] = 0;
    }
}

double JacobiParallelise::resolve()
{
    // Check time
    MPI_Barrier(MPI_COMM_WORLD);
    double timeInit;
    if(_myRank == 0){
        timeInit = MPI_Wtime();
    }
    int k = 0;
    do
    {
        _error = 0; 
        k += 1;

        // Exchanging data
        exchangeData();

        // Compute next step
        jacobi();

        // Exchanging error
        exchangeError();

    } while (k<=MAX_STEP && _error>=MAX_DIFF);
    // Check time
    MPI_Barrier(MPI_COMM_WORLD);
    double timeEnd;
    if(_myRank == 0){
        timeEnd = MPI_Wtime();
    }
    return timeEnd-timeInit;
}

void JacobiParallelise::saveData() const
{
    for (int myRankPrint=0;myRankPrint<_nbTasks;++myRankPrint)
    {
        if (_myRank==myRankPrint)
        {
            ofstream file;
            if(_myRank == 0){
                file.open("jacobiPara.txt", ios::out);
            } else {
                file.open("jacobiPara.txt", ios::app);
            }
            file << *_grid;
            file.close();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void JacobiParallelise::jacobi()
{
    Grid *newGrid;
    newGrid = new Grid(_x0,_y0,_xf,_yf,_nxf,_nyf);
    double constant = _dx*_dx*_dy*_dy/(2.*_dx*_dx+2.*_dy*_dy);
    for (int i = 0;i<_nxf;++i)
    {
        for (int j = 0;j<_nyf;++j)
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
            if (i+1<_nxf)
            {
                newGrid->get(i,j) += constant*_grid->get(i+1,j)/(_dx*_dx);
            }
            else
            {
                newGrid->get(i,j) += constant*_right[j]/(_dx*_dx);
            }
            if (j+1<_nyf)
            {
                newGrid->get(i,j) += constant*_grid->get(i,j+1)/(_dy*_dy);
            }
            else
            {
                newGrid->get(i,j) += constant*_up[i]/(_dy*_dy);
            }
            if (MODE == 1)
            {
                newGrid->get(i,j) -= constant*f1(i*_dx+_x0,j*_dy+_y0);
            }
            else if (MODE == 2)
            {
                newGrid->get(i,j) -= constant*f2(i*_dx+_x0,j*_dy+_y0);
            }
            _error += (newGrid->get(i,j)-_grid->get(i,j))*(newGrid->get(i,j)-_grid->get(i,j));
        }
    }
    delete _grid;
    _grid = newGrid;
}

void JacobiParallelise::exchangeError()
{
    double error;
    MPI_Allreduce(&_error, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    _error = error;
    _error = sqrt(_error)/((_nx-2)*(_ny-2));
}

void JacobiParallelise::exchangeData()
{
    // MPI Exchanges (with left side)
    MPI_Request reqSendLeft, reqRecvLeft;
    if(_myRank > 0)
    {
        double toSend[_nyf];
        _grid->getSide(toSend, LEFT);
        MPI_Isend(&toSend, _nyf, MPI_DOUBLE, _myRank-1, 0, MPI_COMM_WORLD, &reqSendLeft);
        MPI_Irecv(_left, _nyf, MPI_DOUBLE, _myRank-1, 0, MPI_COMM_WORLD, &reqRecvLeft);
    }


    // MPI Exchanges (with right side)
    MPI_Request reqSendRight, reqRecvRight;
    if(_myRank < _nbTasks-1)
    {
        double toSend[_nyf];
        _grid->getSide(toSend, RIGHT);
        MPI_Isend(&toSend, _nyf, MPI_DOUBLE, _myRank+1, 0, MPI_COMM_WORLD, &reqSendRight);
        MPI_Irecv(_right, _nyf, MPI_DOUBLE, _myRank+1, 0, MPI_COMM_WORLD, &reqRecvRight);    
    }


    
    // MPI Exchanges (check everything is send/recv)
    if(_myRank > 0){
        MPI_Wait(&reqSendLeft, MPI_STATUS_IGNORE);
        MPI_Wait(&reqRecvLeft, MPI_STATUS_IGNORE);
    }
    if(_myRank < _nbTasks-1){
        MPI_Wait(&reqSendRight, MPI_STATUS_IGNORE);
        MPI_Wait(&reqRecvRight, MPI_STATUS_IGNORE);
    }
}

/*Grid JacobiParallelise::gather()
{
    if (myRank == 0)
    {
        Grid grid(0,0,_a,_b,_nx,_ny);
        for (int i = 0;i<_nxf;++i)
        {
            for (int j = 0;j<_nyf;++j)
            {
                
            }
        }
        for (int rank = 1;rank<_nbTasks;++rank)
        {

        }
    }
}*/

double JacobiParallelise::f1(const double &x, const double &y) const
{
    return 0;
}

double JacobiParallelise::f2(const double &x, const double &y) const
{
    return -M_PI*M_PI/(_b*_b)*sin(M_PI*(1-y/_b));
}