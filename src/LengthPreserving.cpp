#include "LengthPreserving.h"
#include "CurveFlow.h"

#include <cmath>

LengthPreserving::LengthPreserving(int dimension, int n)
:CurveFlow(dimension, n)
{
    k = new double[(n+2)*dimension]();
}

void LengthPreserving::getRightHandSide(const double& t, double* u, double* fu)
{
    getPartialLength(u);
    getCurvature(u);
    getForce();
    
    double norm(0);
    double* N = new double[dimension]();
    double normN = 0;

    for(int i = 1; i<n+1; i++)
    {
        // normN = 0;
        // for(int j = 0; j<dimension; j++)
        // {
        //     N[j] = ((u[(i+1)*dimension + j] - u[i*dimension + j])/d[i+1] -
        //         (u[i*dimension + j] - u[(i-1)*dimension + j])/d[i]);
        //     normN += N[j]*N[j];
        // }
        // normN = sqrt(normN);

        for(int j = 0; j<dimension; j++)
        {
            fu[i*dimension + j] = 2.0/(d[i+1]*d[i+1]+d[i]*d[i])*
                                (u[(i+1)*dimension+j] - 2*u[i*dimension+j] + u[(i-1)*dimension+j]);
            //fu[i*dimension +j] += F*(1.0/normN)*N[j];
        }

        /**
        * norm now means normal vector
        */
        F = 5;
        norm = (u[i*dimension+1]-u[(i-1)*dimension+1])/d[i];
        fu[dimension*i] += F*norm;

        norm = -(u[i*dimension]-u[(i-1)*dimension])/d[i];
        fu[dimension*i+1] += F*norm;
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

void LengthPreserving::getCurvature(const double* u)
{
    double aux(0.0);
    /**
     * 3d form of curvature
     */
    for(int i = 1; i<n+1; i++)
    {
        k[i] = 0;
        for(int j = 0; j < dimension; j++)
        {
            aux = (
                (u[(i+1)*dimension+j] - u[i*dimension+j])/d[i+1] - 
                (u[i*dimension+j] - u[(i-1)*dimension+j])/d[i]
            );
            k[i] += aux*aux;
        }
        k[i] = std::sqrt(k[i]);
        k[i] *= (2.0/(d[i+1] + d[i]));
    }

    /**
     * 2d form of curvature
     */
    //double norm[] = {0, 0};
    //for(int i = 1; i<n+1; i++)
    //{
    //    k[i] = 0;

    //    norm[0] = (u[i*dimension+1]-u[(i-1)*dimension+1])/d[i];

    //    norm[1] = -(u[i*dimension]-u[(i-1)*dimension])/d[i];

    //    for(int j = 0; j < dimension; j++)
    //    {
    //        aux = (
    //            (u[(i+1)*dimension+j] - u[i*dimension+j])/d[i+1] - 
    //            (u[i*dimension+j] - u[(i-1)*dimension+j])/d[i]
    //        );
    //        k[i] += aux*norm[j];
    //    }
    //    k[i] *= (2.0/(d[i+1] + d[i]));
    //}

    /**
    * boundary condition 
    */
    k[0] = k[n];
    k[n+1] = k[1];
}

void LengthPreserving::getForce()
{
    double elastic = 0;
    double curvature = 0;
    for(int i = 1; i<n+1; i++)
    {
        elastic += k[i]*k[i]*d[i];
        curvature += k[i]*d[i];
    }
    F = elastic/(2*3.1415926);
}

LengthPreserving::~LengthPreserving()
{
    if(this->k) delete[] k;
}
