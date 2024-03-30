#include "CurveFlow.h"

#include <cmath>
#include <iostream>
#include <fstream>

CurveFlow::CurveFlow(int dimension, int n)
{
    this->dimension = dimension;
    this->n = n;
    this->h = 1.0/n;
    this->d = new double[dimension*(n+2)]();
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
    getPartialLength(u);

    for(int j = 0; j<dimension; j++)
    {
        for(int i = 1; i<n+1; i++)
        {
            fu[dimension*i + j] = 2.0 / (d[i]+d[i+1]) * 
                            ((u[(i+1)*dimension+j] - u[i*dimension+j])/d[i+1] -
                            (u[i*dimension+j] - u[(i-1)*dimension+j])/d[i]);

        }
        /**
        * period boundary condition
        */
        fu[j] = fu[n*dimension + j];
        fu[(n+1)*dimension + j] = fu[dimension + j];
    }
}

void CurveFlow::getPartialLength(const double* u)
{
    for(int i = 1; i<n+1; i++)
    {
        d[i] = 0;
        for(int j = 0; j<dimension; j++)
        {
            d[i] += (u[i*dimension + j] - u[(i-1)*dimension+j])*
                    (u[i*dimension + j] - u[(i-1)*dimension+j]);
        }
        d[i] = std::sqrt(d[i]);
    }

    /**
    * boundary condition
    */
    d[0] = d[n];
    d[n+1] = d[1];

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
    if(this->d) delete[] d;
}
