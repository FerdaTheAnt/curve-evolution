#pragma once

#include "ODEProblem.h"

class ODESolver
{
    public:
        virtual bool setup(const int degreesOfFreedom) = 0;

        virtual bool solve(
            const double integrationTimeStep,
            const double stopTime,
            double *time,
            ODEProblem *problem,
            double *u
        ) = 0;
};