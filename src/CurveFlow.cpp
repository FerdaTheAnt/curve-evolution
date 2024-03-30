#include "CurveFlow.h"

#include <cmath>
#include <iostream>
#include <fstream>

CurveFlow::CurveFlow(int dimension, int n)
{
    this->dimension = dimension;
    this->n = n;
    this->h = 1.0/n;
    this->bd = new double[dimension*(n+2)]();
    this->g = new double[n+2]();
}

int CurveFlow::getDegreesOfFreedom()
{
    return dimension*(n+2);
}

void CurveFlow::setInitialCondition(double* u, void (*func)(double*, double))
{
    /**
    * watch out for memory leaks and array size
    */
    double *temp = new double[this->dimension]();
    for (int i = 0; i<n; i++) 
    {
        func(temp, i*h);
        for(int j = 0; j<dimension; j++)
        {
            u[i*dimension + j] = temp[j];
        }
    }
    /**
    * period boundary conditions
    */
    for(int j = 0; j<dimension; j++)
    {
        u[n*dimension + j] = u[j];
        u[(n+1)*dimension + j] = u[dimension + j];
    }

    delete[] temp;
}

void CurveFlow::getRightHandSide(const double& t, double* u, double* fu)
{
    getBackDifference(u);
    getLocalLength();

    for(int j = 0; j<dimension; j++)
    {
        for(int i = 0; i<n; i++)
        {
            fu[dimension*i + j] = (bd[(i+1)*dimension + j]/g[i+1]-bd[i*dimension + j]/g[i])/(g[i]*h);
        }
        /**
        * period boundary condition
        */
        fu[n*dimension + j] = fu[j];
        fu[(n+1)*dimension + j] = fu[dimension + j];
    }
}

void CurveFlow::getBackDifference(const double* u)
{
    for(int j = 0; j<dimension; j++)
    {
        for(int i = 1; i<n+1; i++)
        {
            bd[i*dimension + j] = (u[dimension*i+j] - u[dimension*(i-1)+j])/h;
        }
        bd[j] = bd[n*dimension + j];
        bd[(n+1)*dimension + j] = bd[dimension + j];
    }
}

void CurveFlow::getLocalLength()
{
    for(int i = 1; i<n+1; i++)
    {
        g[i] = 0;
        for(int j = 0; j<dimension; j++)
        {
            g[i] += bd[i*dimension + j] * bd[i*dimension + j];
        }
        g[i] = sqrt(g[i]);
    }

    g[0] = g[n];
    g[n+1] = g[1];
}

bool CurveFlow::writeSolution(const double &t, int step, const double *u)
{
    std::fstream file;
    if(step == 0)
        file.open("curve.dat", std::ios::out | std::ios::trunc);
    else
        file.open("curve.dat", std::ios::out | std::ios::app);
    if(!file)
    {
        std::cerr << "Unable to open file curve.dat" <<std::endl;
        return false;
    }
    file << "#t = " << t << "; step = " << step << "\n";
    for(int i = 0; i<n+1; i++)
    {
        for(int j = 0; j<dimension; j++)
        {
            file << u[i*dimension + j] << " "; 
        }
        file << "\n";
    }
    file<<"\n\n\n";
    return true;
}

CurveFlow::~CurveFlow()
{
    if(this->bd) delete[] bd;
    if(this->g) delete[] g;
}
