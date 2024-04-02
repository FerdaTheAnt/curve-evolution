#include "DeTurck.h"
#include "CurveFlow.h"

DeTurck::DeTurck(int dimension, int n) : CurveFlow(dimension, n){}

void DeTurck::getRightHandSide(const double& t, double* u, double* fu)
{
    double norm(0);
    getPartialLength(u);

    for(int j = 0; j<dimension; j++)
    {
        for(int i = 1; i<n+1; i++)
        {
            norm = 2.0/(d[i]*d[i] + d[i+1]*d[i+1]);
            fu[dimension*i + j] = norm*(u[(i+1)*dimension+j] - 2*u[i*dimension+j] + u[(i-1)*dimension+j]);
        }
        /**
        * period boundary condition
        */
        fu[j] = fu[n*dimension + j];
        fu[(n+1)*dimension + j] = fu[dimension + j];
    }
}
