#pragma once

#include "CurveFlow.h"

class AreaPreserving : public CurveFlow
{
    public:
        AreaPreserving(int dimension, int n);
        void getRightHandSide(const double& t, double* u, double* fu) override;
    protected:
        double force;
        void getAreaForce();
        double getLength();
};