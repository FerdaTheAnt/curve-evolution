#include "DeTurck.h"
#include "CurveFlow.h"

DeTurck::DeTurck(int dimension, int n) : CurveFlow(dimension, n){}

void DeTurck::getRightHandSide(const double& t, double* u, double* fu)
{
    double norm(0);
    getBackDifference(u);
    getLocalLength();

    for(int j = 0; j<dimension; j++)
    {
        for(int i = 0; i<n; i++)
        {
            norm = 0.5*(g[i]*g[i] + g[i+1]*g[i+1]);
            fu[dimension*i + j] = (bd[(i+1)*dimension + j]-bd[i*dimension + j])/(norm*h);
        }
        /**
        * period boundary condition
        */
        fu[n*dimension + j] = fu[j];
        fu[(n+1)*dimension + j] = fu[dimension + j];
    }
}
