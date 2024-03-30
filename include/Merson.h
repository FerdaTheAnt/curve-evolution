#pragma once

#include "ODEProblem.h"
#include "ODESolver.h"

class Merson : public ODESolver
{
    public:
        Merson();
        bool setup(int dof);
        void setAdaptivity(const double& adaptivity);
        bool solve(
            const double integrationTimeStep,
            const double stopTime,
            double *time,
            ODEProblem *problem,
            double* u
        );
        ~Merson();
    protected:
        double *k1, *k2, *k3, *k4, *k5, *aux;
        double adaptivity;
};