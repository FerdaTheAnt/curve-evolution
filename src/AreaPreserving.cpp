#include "AreaPreserving.h"
#include "CurveFlow.h"
#include <cmath>

const int DIM = 2; //dimension for area preserving flow is just 2
const double PI = 3.1415926535;

AreaPreserving::AreaPreserving(int dimension, int n)
:CurveFlow(dimension, n)
{
}

void AreaPreserving::getRightHandSide(const double& t, double* u, double* fu)
{
    double norm(0);
    getPartialLength(u);

    getAreaForce();

    double *N = new double[dimension]();
    double normN = 0;

    for(int i = 1; i<n+1; i++)
    {
        
        normN = 0;
        for(int j = 0; j<dimension; j++)
        {
            N[j] = ((u[(i+1)*dimension + j] - u[i*dimension + j])/d[i+1] -
                (u[i*dimension + j] - u[(i-1)*dimension + j])/d[i]);
            normN += N[j]*N[j];
        }
        normN = sqrt(normN);

        for(int j = 0; j<dimension; j++)
        {
            norm = 2.0/(d[i]*d[i] + d[i+1]*d[i+1]);
            fu[dimension*i + j] = norm*(u[(i+1)*dimension + j]-2.0*u[i*dimension + j]+u[(i-1)*dimension+j]);
            fu[i*dimension + j] += force*(1.0/normN)*N[j];
        }
        ///**
        //* norm now means normal vector
        //*/
        //norm = (u[i*dimension+1]-u[(i-1)*dimension+1])/d[i];
        //fu[dimension*i] += force*norm;

        //norm = -(u[i*dimension]-u[(i-1)*dimension])/d[i];
        //fu[dimension*i+1] += force*norm;
    }

    /**
    * period boundary condition
    */
    for(int j = 0; j<dimension; j++)
    {
        fu[j] = fu[n*dimension+j];
        fu[(n+1)*dimension + j] = fu[dimension + j];
    }
}

void AreaPreserving::getAreaForce()
{
    double length = getLength();
    force = -2*PI/length;
}

double AreaPreserving::getLength()
{
    double length = 0;
    for(int i = 1; i<n+1; i++){
        length += d[i];
    }
    return length;
}
