#include "RedistributionInSpace.h"
#include "CurveFlow.h"

#include <cmath>

RedistributionSpace::RedistributionSpace(int dimension, int n)
:CurveFlow(dimension, n)
{
    length = 0.0;
    omega = 0.0;
    k = new double[n+2]();
    alpha = new double[n+2]();
}


void RedistributionSpace::getRightHandSide(const double& t, double* u, double* fu)
{
    getPartialLength(u);
    getLength();
    getCurvature(u);
    getOmega();
    getSeries();
    getSum();
    getAlpha();

    for(int j = 0; j<dimension; j++)
    {
        for(int i = 1; i<n+1; i++)
        {
            fu[dimension*i + j] = alpha[i]*(u[(i+1)*dimension+j] - u[(i-1)*dimension+j])/(d[i+1] + d[i]);
            fu[dimension*i + j] += 2.0/(d[i+1]+d[i])*
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

void RedistributionSpace::getCurvature(const double* u)
{
    double aux(0.0);
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
    * boundary condition 
    */
    k[0] = k[n];
    k[n+1] = k[1];
}

void RedistributionSpace::getOmega()
{
    const double kappa2 = 5.0;
    omega = 0.0;
    for(int i = 1; i < n+1; i++)
    {
        omega += k[i]*k[i]*d[i];
    }
    omega *= (kappa2/length);

    //TODO: consider removing this, it was used only for debugging
    //if(STEPS < 6)
    //    std::cout<< "omega is " << omega << "\n";
    //STEPS++;
}

void RedistributionSpace::getSeries()
{
    double averageElastic (0.0);
    for(int i = 1; i<n+1; i++)
    {
        averageElastic += k[i]*k[i]*d[i];
    }
    averageElastic *= (1.0/length);
    for(int i = 1; i < n+1; i++)
    {
        alpha[i] = k[i]*k[i]*d[i] - averageElastic*d[i] + omega*(length/n - d[i]);
        //alpha[i] *= 0.5*(d[i+1]+d[i]);
    }

    /* 
    * boundary condition 
    */
    alpha[0] = alpha[n];
    alpha[n+1] = alpha[1];

    //TODO: consider removing this, it was used only for debugging
    //if(STEPS < 6)
    //    std::cout<< "a series is " << alpha[0] << " and " << alpha[2] << "\n";
    //STEPS++;
}

void RedistributionSpace::getSum()
{
    double sum(0.0);
    //double lastAlph = alpha[1];
    for(int i = 2; i < n+1; i++)
    {
        sum += alpha[i];
        //lastAlph = alpha[i];
        alpha[i] = sum;
    }

    //TODO: consider removing this, it was used only for debugging
    //if(STEPS < 6)
    //    std::cout<< "sum series is " << alpha[0] << " and " << alpha[2] << "\n";
    //STEPS++;
}

void RedistributionSpace::getAlpha()
{
    /**
    * compute first alpha
    */
    double aux(0.0);
    double auxSumLength(0.0);
    auxSumLength += 0.5*(d[2] + d[1]);
    for(int i = 2; i<n+1; i++)
    {
        aux += 0.5*alpha[i]*(d[i+1]+d[i]);
        auxSumLength += 0.5*(d[i+1]+d[i]);
    }
    alpha[1] = -aux/auxSumLength;
 
    for(int i = 2; i<n+1; i++)
    {
        alpha[i] += alpha[1];
    }
    
    /**
    * boundary condition
    */
    alpha[0] = alpha[n];
    alpha[n+1] = alpha[1];
}

RedistributionSpace::~RedistributionSpace()
{
    if(k != nullptr) delete k;
    if(alpha != nullptr) delete alpha;
}
