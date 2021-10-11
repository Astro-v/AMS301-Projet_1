
#ifndef __GRID_HPP__
#define __GRID_HPP__

    #include <iostream>
    #include "constant.hpp"

    class Grid
    {
    private:
        int _nx, _ny;
        double _x0, _y0, _xf, _yf;
        double _dx2, _dy2;
        double** _grid;

    public:
        // -------- CREATOR -------- //
        Grid();
        Grid(double x0, double y0, double xf, double yf, int nx, int ny);
        Grid(double x0, double y0, double xf, double yf, int nx, int ny, double vec[]);
        
        // -------- DESTRUCTOR -------- //
        ~Grid();

        // -------- ACCESSOR -------- //
        int dim(int const&) const; // READ only (_nx, _ny)
        void getSide(double vec[], const Side side); // READ only (_grid)
        void getSide(double vec[], const Side side, const int step);

        // -------- OPERATOR -------- //
        double& operator ()(int const&, int const&); // READ & WRITE (_grid)
        double operator ()(int const&, int const&) const; // READ only (_grid)
        double& get(int const&, int const&); // READ & WRITE (_grid)

        // -------- OTHER -------- //   
        // DISPLAY    
        void print();

        void toVector(double* vec) const;
        
    };

    // -------- ACCESSOR -------- //
    // FLUX OPERATOR
    std::ostream& operator <<(std::ostream&, Grid const&);


#endif // __GRID_HPP__