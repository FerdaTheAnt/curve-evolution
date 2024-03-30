#pragma once

#include "ODEProblem.h"

class CurveFlow : public ODEProblem
{
    public:
        CurveFlow(int dimension, int n);
        int getDegreesOfFreedom();
        void setInitialCondition(double* u, void (*func)(double*, double));
        virtual void getRightHandSide(const double& t, double* u, double* fu);
        bool writeSolution(const double &t, int step, const double *u);
        virtual ~CurveFlow();
    protected:
        int dimension;
        int n;
        double h;
        double* d;
        void getPartialLength(const double* u);
};
