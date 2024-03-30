#pragma once

#include "ODEProblem.h"

class CurveFlowForce : public ODEProblem
{
    public:
        CurveFlowForce(int dimension, int n);
        int getDegreesOfFreedom();
        void setInitialCondition(double* u, void (*func)(double*, double));
        virtual void getRightHandSide(const double& t, double* u, double* fu);
        bool writeSolution(const double &t, int step, const double *u);
        virtual ~CurveFlowForce();
    protected:
        int dimension;
        int n;
        double h;
        double* d;
        
        //void getBackDifference(const double* u);
        void getPartialLength(const double* u);
        //void getLocalLength();
    private:
        //void (*force)(double*, double, double);
        double* F;
        void getForce();
};
