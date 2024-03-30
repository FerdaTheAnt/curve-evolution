#include "Merson.h"

#include "ODEProblem.h"

#include <iostream>
#include <cmath>

Merson::Merson()
{
    k1 = 0;
    k2 = 0;
    k3 = 0;
    k4 = 0;
    k5 = 0;
    aux = 0;

    adaptivity = 10e-6;
}

bool Merson::setup(const int dof)
{
    k1 = new double[dof];
    k2 = new double[dof];
    k3 = new double[dof];
    k4 = new double[dof];
    k5 = new double[dof];
    aux = new double[dof];
    if( ! k1 || ! k2 || ! k3 || ! k4 || ! k5 || ! aux )
      return false;
   return true;
}

void Merson::setAdaptivity(const double& adaptivity)
{
    this->adaptivity = adaptivity;
}

bool Merson::solve(
    const double integrationTimeStep,
    const double stopTime,
    double *time,
    ODEProblem *problem,
    double *u
)
{
    const int dof = problem->getDegreesOfFreedom();
    double tau = std::min(integrationTimeStep, stopTime - *time);
    long int iteration(0);
    while(*time < stopTime)
    {
        /**
        * Compute k1
        */
        problem->getRightHandSide(*time, u, k1);

        /**
        * Compute k2
        */
        for(int i = 0; i<dof; i++)
            aux[i] = u[i] + tau * 1.0/3.0 * k1[i];
        problem->getRightHandSide(*time + 1.0/3.0 * tau, aux, k2);

        /**
        * Compute k3
        */
        for(int i = 0; i < dof; i++)
            aux[i] = u[i] + tau * 1.0/6.0 * (k1[i] + k2[i]);
        problem->getRightHandSide(*time + 1.0/3.0 * tau, aux, k3);

        /**
        * Compute k4
        */
        for(int i = 0; i < dof; i++)
            aux[i] = u[i] + tau * (0.125 * k1[i] + 0.375 * k3[i]);
        problem->getRightHandSide(*time + 1.0/2.0 * tau, aux, k4);

        /**
        * Compute k5
        */
        for(int i = 0; i < dof; i++)
            aux[i] = u[i] + tau * (0.5 * k1[i] - 1.5 * k3[i] + 2.0 * k4[i]);
        problem->getRightHandSide(*time + tau, aux, k5);

        /**
        * Compute error
        */
        double eps(0.0);
        for(int i = 0; i < dof; i++)
        {
            double err = tau/3.0 * std::fabs(0.2*k1[i] - 0.9*k3[i] + 0.8*k4[i] - 0.1*k5[i]);
            eps = std::max(eps, err);
        }

        if(this->adaptivity == 0.0 || eps < this->adaptivity)
        {
            for(int i = 0; i<dof; i++)
                u[i] += tau/6.0 * (k1[i] + 4.0*k4[i] + k5[i]);
            *time += tau;
            iteration++;
            if(iteration > 100000)
            {
                std::cerr << "The solver has reached maximum number of iterations. " << std::endl;
                return false;
            }
        }
        if(this->adaptivity)
            tau *= 0.8 * std::pow(this->adaptivity/eps, 0.2);
        tau = std::min(tau, stopTime - *time);
        //std::cout << "ITER: " << iteration << " \t tau = " << tau << " \t time = " << *time << "         \r " << std::flush;
    }
    //std::cout << std::endl;
    return true;
}

Merson::~Merson()
{
    if(k1) delete k1;
    if(k2) delete k2;
    if(k3) delete k3;
    if(k4) delete k4;
    if(k5) delete k5;
    if(aux) delete aux;
}