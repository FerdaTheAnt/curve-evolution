#include "Euler.h"

#include <algorithm>

Euler::Euler()
{
    this->k = 0;
}

bool Euler::setup( const int dof )
{
   this->k = new double[dof];
   if(!this->k)
      return false;
   return true;
}

bool Euler::solve( 
    const double integrationTimeStep,
    const double stopTime,
    double* time,
    ODEProblem* problem,
    double* u 
)
{
   const int dof = problem->getDegreesOfFreedom();
   long int iteration(0);
   while(*time < stopTime)
   {
     const double tau = std::min(integrationTimeStep, stopTime - *time);
     problem->getRightHandSide(*time, u, this->k);
     for(int i = 0; i < dof; i++)
        u[i] +=  tau * k[i];
      *time += tau;
      iteration++;
      //std::cout << "ITER: " << iteration << " \t tau = " << tau << " \t time= " << *time << "         \r " << std::flush;
   }
   //std::cout << std::endl;
   return true;
}

Euler::~Euler()
{
   if(k) delete[] k;
}
