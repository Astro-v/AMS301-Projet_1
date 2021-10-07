#ifndef __JACOBIPARALLELISE_HPP__
#define __JACOBIPARALLELISE_HPP__

#include "grid.hpp"
#include "constant.hpp"

class JacobiParallelise
{
private:
    double _a, _b;
    double _x0, _y0, _xf, _yf;
    double _u0, _alpha;
    int _nx, _ny, _nxf, _nyf;
    double _dx, _dy;
    Grid *_grid;
    int _myRank, _nbTasks;
    double *_left, *_right, *_up, *_down;

    // -------- OTHER -------- // 
    void init();
    void jacobi(); // Jacobi method
    void exchangeData(); // Exchange data with the others process

public:
    // -------- CREATOR -------- //
    JacobiParallelise(int argc, char* argv[]);

    // -------- DESTRUCTOR -------- //
    ~JacobiParallelise();

    // -------- ACCESSOR -------- //
    int getRank() const;

    // -------- OTHER -------- // 
    double resolve();
    void saveData() const;
    /*Grid gather();*/
};

#endif // __JACOBIPARALLELISE_HPP__