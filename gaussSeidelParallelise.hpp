#ifndef __GAUSSSEIDELPARALLELISE_HPP__
#define __GAUSSSEIDELPARALLELISE_HPP__

#include "grid.hpp"
#include "constant.hpp"

class GaussSeidelParallelise
{
private:
    double _a, _b;
    double _x0, _y0, _xf, _yf;
    double _u0, _alpha;
    int _nx, _ny, _nxf, _nyf;
    double _dx, _dy;
    Grid *_grid;
    int _myRank, _nbTasks;
    double *_left, *_right;
    double *_up, *_down;
    int _parity; // Parity of the pos (0,0) of _grid
    double _error;

    // -------- OTHER -------- // 
    void init1();
    void init2();
    void gaussSeidel(int step); // Gauss Seidel method
    void exchangeData(); // Exchange data with the others process
    void exchangeError(); // Exchange error with the others process
    double f1(const double &x, const double &y) const;
    double f2(const double &x, const double &y) const;

public:
    // -------- CREATOR -------- //
    GaussSeidelParallelise(int argc, char* argv[]);

    // -------- DESTRUCTOR -------- //
    ~GaussSeidelParallelise();

    // -------- ACCESSOR -------- //
    int getRank() const;

    // -------- OTHER -------- // 
    double resolve();
    void saveData() const;
    /*Grid gather();*/
};

#endif // __GAUSSSEIDELPARALLELISE_HPP__