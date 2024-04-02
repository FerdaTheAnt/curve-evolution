#include "CurveFlowForce.h"
#include "CurveFlow.h"
#include <cmath>

CurveFlowForce::CurveFlowForce(int dimension, int n)
:CurveFlow(dimension, n)
{
    F = new double[(n+2)*dimension]();
}

void CurveFlowForce::getRightHandSide(const double& t, double* u, double* fu)
{
    getPartialLength(u);
    getForce();
    
    double* N = new double[dimension]();
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
            fu[i*dimension + j] = 2.0/(d[i+1]+d[i])*N[j];
            fu[i*dimension +j] += F[i]*(1.0/normN)*N[j];
        }
    }

    /**
        * boundary condition
        */
    for(int j = 0; j<dimension; j++)
    {
        fu[j] = fu[n*dimension + j];
        fu[(n+1)*dimension + j] = fu[dimension + j];
    }

    delete[] N;
}

void CurveFlowForce::getForce()
{
    for(int i = 1; i<n+1; i++)
    {
        F[i] = -10.0;
    }

    /*
    * Boundary condition
    */
    F[0] = F[n];
    F[n+1] = F[1];
}

CurveFlowForce::~CurveFlowForce()
{
    if(this->F) delete[] F;
}
