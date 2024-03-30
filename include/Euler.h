#pragma once

#include "ODEProblem.h"
#include "ODESolver.h"

class Euler : public ODESolver
{
    public:
        Euler();
        bool setup(const int dof);
        bool solve(
            const double integrationTimeStep,
            const double stopTime,
            double *time,
            ODEProblem *problem,
            double* u
        );
        ~Euler();
    protected:
        double *k;
};