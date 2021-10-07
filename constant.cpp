#include <cmath>

double f(double x, double y, double b)
{
    return -M_PI*M_PI/(b*b)*sin(M_PI*(1-y/b));
}