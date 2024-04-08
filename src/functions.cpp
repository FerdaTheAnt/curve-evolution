#include "functions.h"

#include <cmath>

const double PI = 3.1415926535;

/**
* 2D curves
*/
void fold(double *x, double p_u)
{
    x[0] =  (1+0.6*cos(7*2*PI*p_u))*cos(2*PI*p_u);
    x[1] = (1+0.6*cos(7*2*PI*p_u))*sin(2*PI*p_u);
}


void circle(double *x, double p_u)
{
    x[0] = cos(2*PI*p_u);
    x[1] = sin(2*PI*p_u);
}


void lissajous(double *x, double p_u)
{
    x[0] = sin(5*2*PI*p_u + PI/2);
    x[1] = sin(4*2*PI*p_u);
}

/**
* 3D curves
*/
void spherical(double *x, double p_u)
{
    double r = 1/sqrt(1+16*cos(5*2*PI*p_u)*cos(5*2*PI*p_u));
    x[0] = r * cos(2*PI*p_u);
    x[1] = r * sin(2*PI*p_u);
    x[2] = r * 4* cos(5*2*PI*p_u);
}

void wavy_circle(double *x, double p_u)
{
    x[0] = cos(2*PI*p_u);
    x[1] = sin(2*PI*p_u);
    x[2] = cos(12*PI*p_u);
}

void inclined_circle(double *x, double p_u)
{
    x[0] = 0.5*sqrt(2)*cos(2*PI*p_u);
    x[1] = sin(2*PI*p_u);
    x[2] = 0.5*sqrt(2)*cos(2*PI*p_u);
}

double phi(double p_u)
{
    return 0.7*p_u*p_u + 0.25*p_u;
}

void example(double *x, double p_u)
{
    x[0] = cos(phi(2*PI*p_u));
    x[1] = sin(phi(2*PI*p_u));
    x[2] = sin(2*PI*p_u) - 0.5*sin(4*PI*p_u);
}

void torus_helix(double *x, double p_u)
{
    x[0] = (1 + 0.3*cos(2*PI*10*p_u))*cos(2*PI*p_u);
    x[1] = (1 + 0.3*cos(2*PI*10*p_u))*sin(2*PI*p_u);
    x[2] = 0.3*sin(2*PI*10*p_u);
}

void spherical2(double *x, double p_u)
{
    double ni = 2.0*sin(2*2*PI*p_u) - 1.99*pow(sin(2*2*PI*p_u), 3);
    x[0] = cos(ni)*cos(2*PI*p_u);
    x[1] = sin(ni)*cos(2*PI*p_u);
    x[2] = sin(2*PI*p_u);
}

void cyllindrical(double *x, double p_u)
{
    p_u = 2*PI*p_u;
    double r = 1 + 0.1*sin(p_u);
    x[0] = r*cos(4*p_u);
    x[1] = r*cos(p_u);
    x[2] = cos(p_u) - sin(2*p_u);
}
