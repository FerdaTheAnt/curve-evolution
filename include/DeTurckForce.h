#pragma once

#include "CurveFlow.h"

class DeTurckForce : public CurveFlow
{
    public:
        DeTurckForce(int dimension, int n);
        void getRightHandSide(const double& t, double* u, double* fu) override;
        virtual ~DeTurckForce();
    private:
        //void (*force)(double*, double, double);
        double* F;
        void getForce();
};
