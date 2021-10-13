#include <iostream>
#include "grid.hpp"
#include "constant.hpp"

using namespace std;

// -------- CREATOR -------- //

Grid::Grid():
_x0(0), _y0(0),
_xf(1), _yf(1),
_nx(NX_INIT), _ny(NY_INIT)
{
    _dx2 = (_xf-_x0)/_nx;
    _dx2 = _dx2*_dx2;
    _dy2 = (_yf-_y0)/_ny;
    _dy2 = _dy2*_dy2;

    _grid = new double*[_nx];
    for (int i(0);i<_nx;++i)
    {
        _grid[i] = new double[_ny];
        for (int j(0);j<_ny;++j)
        {
            _grid[i][j] = 0;
        }
    }
}

Grid::Grid(double x0, double y0, double xf, double yf, int nx, int ny):
_x0(x0), _y0(y0),
_xf(xf), _yf(yf),
_nx(nx), _ny(ny)
{
    _dx2 = (_xf-_x0)/(_nx-1);
    _dx2 = _dx2*_dx2;
    _dy2 = (_yf-_y0)/(_ny-1);
    _dy2 = _dy2*_dy2;

    _grid = new double*[_nx];
    for (int i(0);i<_nx;++i)
    {
        _grid[i] = new double[_ny];
        for (int j(0);j<_ny;++j)
        {
            _grid[i][j] = 0;
        }
    }
}

Grid::Grid(double x0, double y0, double xf, double yf, int nx, int ny, double vec[]):
_x0(x0), _y0(y0),
_xf(xf), _yf(yf),
_nx(nx), _ny(ny)
{
    _dx2 = (_xf-_x0)/(_nx-1);
    _dx2 = _dx2*_dx2;
    _dy2 = (_yf-_y0)/(_ny-1);
    _dy2 = _dy2*_dy2;

    _grid = new double*[_nx];
    for (int i(0);i<_nx;++i)
    {
        _grid[i] = new double[_ny];
        for (int j(0);j<_ny;++j)
        {
            _grid[i][j] = vec[j+_ny*i];
        }
    }
}

// -------- DESTRUCTOR -------- //

Grid::~Grid()
{
    for (int i(0);i<_nx;++i)
    {
        delete _grid[i];
    }
    delete _grid;
}

// -------- ACCESSOR -------- //

int Grid::dim(int const& i) const
{
    if (i == 1)
    {
        return _nx;
    }
    return _ny;
}

void Grid::getSide(double vec[], const Side side)
{
    switch (side)
    {
    case DOWN:
        for (int i=0;i<_nx;++i)
        {
            vec[i]=_grid[i][0];
        }
        break;
    case UP:
        for (int i=0;i<_nx;++i)
        {
            vec[i]=_grid[i][_ny-1];
        }
        break;
    case LEFT:
        for (int i=0;i<_ny;++i)
        {
            vec[i]=_grid[0][i];
        }
        break;
    case RIGHT:
        for (int i=0;i<_ny;++i)
        {
            vec[i]=_grid[_nx-1][i];
        }
        break;
    }
}

// -------- OPERATOR -------- //

double& Grid::operator ()(int const& i, int const& j)
{
    return _grid[i][j];
}

double Grid::operator ()(int const& i, int const& j) const
{
    return _grid[i][j];
}

double& Grid::get(int const& i, int const& j)
{
    return _grid[i][j];
}

// -------- OTHER -------- // 

void Grid::print()
{
    for (int i=0;i<_nx;++i)
    {
        for (int j=0;j<_ny;++j)
        {
            cout << _grid[i][j] << " ";
        }
        cout << endl;
    }
}

void Grid::toVector(double* vec) const
{
    for (int i=0;i<_nx;++i)
    {
        for (int j=0;j<_ny;++j)
        {
            vec[j+_ny*i] = _grid[i][j];
        }
    }
}

// -------- ACCESSOR -------- //

// FLUX OPERATOR

ostream& operator <<(ostream& stream, Grid const& grid)
{
    for (int i=0;i<grid.dim(1);++i)
    {
        for (int j=0;j<grid.dim(2);++j)
        {
            stream << grid(i,j) << ' ';
        }
        stream << endl;
    }
    return stream;
}