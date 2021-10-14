#ifndef __JACOBISEQUENTIEL_HPP__
#define __JACOBISEQUENTIEL_HPP__

#include "grid.hpp"

class JacobiSequentiel
{
private:
    double _a, _b;
    double _u0, _alpha;
    int _nx, _ny;
    double _dx, _dy;
    Grid *_grid;
    double *_left, *_right, *_up, *_down;
    int _myRank, _nbTasks;

    // -------- OTHER -------- // 
    double jacobi();
    void initLimitCondition1();
    void initLimitCondition2();
    double f1(const double &x, const double &y) const;
    double f2(const double &x, const double &y) const;

public:
    // -------- CREATOR -------- //
    JacobiSequentiel(int argc, char* argv[]);

    // -------- DESTRUCTOR -------- //
    ~JacobiSequentiel();

    // -------- ACCESSOR -------- //
    int getRank() const;

    // -------- OTHER -------- // 
    double resolve();
    void saveData() const;
    double error() const;
    
};




#endif // __JACOBISEQUENTIEL_HPP__