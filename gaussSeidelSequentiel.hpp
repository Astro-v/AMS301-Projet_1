#ifndef __GAUSSSEIDELSEQUENTIEL_HPP__
#define __GAUSSSEIDELSEQUENTIEL_HPP__

#include "grid.hpp"

class GaussSeidelSequentiel
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
    double gaussSeidel();
    void initLimitCondition1();
    void initLimitCondition2();
    double f1(const double &x, const double &y) const;
    double f2(const double &x, const double &y) const;

public:
    // -------- CREATOR -------- //
    GaussSeidelSequentiel(int argc, char* argv[]);

    // -------- DESTRUCTOR -------- //
    ~GaussSeidelSequentiel();

    // -------- ACCESSOR -------- //
    int getRank() const;

    // -------- OTHER -------- // 
    double resolve();
    void saveData() const;
    
};


#endif // __GAUSSSEIDELSEQUENTIEL_HPP__