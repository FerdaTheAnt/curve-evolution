#pragma once

#include "CurveFlow.h"

class CurveFlowForce : public CurveFlow
{
    public:
        CurveFlowForce(int dimension, int n);
        virtual void getRightHandSide(const double& t, double* u, double* fu) override;
        virtual ~CurveFlowForce();
    private:
        //void (*force)(double*, double, double);
        double* F;
        void getForce();
};
