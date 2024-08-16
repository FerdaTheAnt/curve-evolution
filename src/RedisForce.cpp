#include "RedisForce.h"
#include "CurveFlow.h"

#include <cmath>

RedisForce::RedisForce(int dimension, int n)
:CurveFlow(dimension, n)
{
    length = 0.0;
    omega = 0.0;
    alpha = new double[n+2]();
    d = new double[n+2]();
    k = new double[n+2]();

    F = new double[n+2]();    
}

void RedisForce::getRightHandSide(const double& t, double* u, double* fu)
{
    getPartialLength(u);
    getLength();
    getCurvature(u);
    getForce(u);
    getOmega();
    getSeries();
    getSum();
    getAlpha();

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
            fu[dimension*i + j] = alpha[i]*(u[(i+1)*dimension+j] - u[(i-1)*dimension+j])/(d[i+1]+d[i]);
            fu[dimension*i + j] += 2.0/(d[i]+d[i+1])*N[j];
            fu[dimension*i + j] += F[i]*(1.0/normN)*N[j];
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

void RedisForce::getPartialLength(const double* u)
{
    for(int i = 1; i < n+1; i++)
    {
        d[i] = 0;
        for(int j = 0; j<dimension; j++)
        {
            d[i] += (u[i*dimension + j]-u[(i-1)*dimension+j])*
                (u[i*dimension + j]-u[(i-1)*dimension+j]);
        }
        d[i] = std::sqrt(d[i]);
    }

    /**
    * boundary condition
    */
    d[0] = d[n];
    d[n+1] = d[1];
}

void RedisForce::getLength()
{
    length = 0;
    for(int i = 1; i < n+1; i++)
        length += d[i];
}

void RedisForce::getCurvature(const double* u)
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

void RedisForce::getForce(const double* u)
{
    for(int i = 1; i < n+1; i++)
    {
        F[i] = -5.0; //shape of force function in here
    }
    F[0] = F[n];
    F[n+1] = F[1];
}

void RedisForce::getOmega()
{
    const double kappa2 = 5.0;
    omega = 0.0;
    for(int i = 1; i < n+1; i++)
    {
        omega += (k[i] + F[i])*k[i]*d[i];
    }
    omega *= (kappa2/length);
}

void RedisForce::getSeries()
{
    /*
    * Careful choice of derivative norm aproximation is necessary
    */
    double averageElastic (0.0);
    for(int i = 1; i<n+1; i++)
    {
        averageElastic += (k[i] + F[i])*k[i]*d[i];
    }
    averageElastic *= (1.0/length);
    for(int i = 1; i < n+1; i++)
    {
        alpha[i] = (k[i]+F[i])*k[i]*d[i] - averageElastic*d[i] + omega*(length/n - d[i]);
    }

    /* 
    * boundary condition 
    */
    alpha[0] = alpha[n];
    alpha[n+1] = alpha[1];
}

void RedisForce::getSum()
{
    double sum(0.0);
    //double lastAlph = alpha[1];
    for(int i = 2; i < n+1; i++)
    {
        sum += alpha[i];
        //lastAlph = alpha[i];
        alpha[i] = sum;
    }
}

void RedisForce::getAlpha()
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

RedisForce::~RedisForce()
{
    if(d != nullptr) delete d;
    if(k != nullptr) delete k;
    if(alpha != nullptr) delete alpha;
    if(F != nullptr) delete F;
}
