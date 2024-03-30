#include "AreaPreserving.h"
#include "CurveFlow.h"

const int DIM = 2; //dimension for area preserving flow is just 2
const double PI = 3.1415926535;

AreaPreserving::AreaPreserving(int dimension, int n)
:CurveFlow(dimension, n)
{
}

void AreaPreserving::getRightHandSide(const double& t, double* u, double* fu)
{
    double norm(0);
    getBackDifference(u);
    getLocalLength();

    getAreaForce();
    for(int i = 0; i<n; i++)
    {
        for(int j = 0; j<dimension; j++)
        {
            norm = 0.5*(g[i]*g[i] + g[i+1]*g[i+1]);
            fu[dimension*i + j] = (bd[(i+1)*dimension + j]-bd[i*dimension + j])/(norm*h);
        }
        /**
        * norm now means normal vector
        */
        norm = bd[i*dimension+1]/g[i];
        fu[dimension*i] += force*norm;

        norm = -bd[i*dimension]/g[i];
        fu[dimension*i+1] += force*norm;
    }

    /**
    * period boundary condition
    */
    for(int j = 0; j<dimension; j++)
    {
        fu[n*dimension + j] = fu[j];
        fu[(n+1)*dimension + j] = fu[dimension + j];
    }
}

void AreaPreserving::getAreaForce()
{
    double length = getLength();
    force = 2*PI/length;
}

double AreaPreserving::getLength()
{
    double length = 0;
    for(int i = 0; i<n; i++){
        length += h*g[i];
    }
    return length;
}
