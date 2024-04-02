#include "CurveFlowForce.h"
#include "CurveFlow.h"

CurveFlowForce::CurveFlowForce(int dimension, int n)
:CurveFlow(dimension, n)
{
}

void CurveFlowForce::getRightHandSide(const double& t, double* u, double* fu)
{
    getPartialLength(u);
    
    for(int j = 0; j<dimension; j++)
    {
        for(int i = 1; i<n+1; i++)
        {
            fu[dimension*i + j] = 2.0/(d[i+1]+d[i])*
                                    ((u[(i+1)*dimension + j] - u[i*dimension + j])/d[i+1] -
                                     (u[i*dimension + j] - u[(i-1)*dimension + j])/d[i]);
        }
        /**
        * boundary condition
        */
        fu[j] = fu[n*dimension + j];
        fu[(n+1)*dimension + j] = fu[dimension + j];
    }
}

CurveFlowForce::~CurveFlowForce()
{
}
