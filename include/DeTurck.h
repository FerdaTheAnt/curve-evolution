#pragma once

#include "CurveFlow.h"

class DeTurck : public CurveFlow
{
    public:
        DeTurck(int dimension, int n);
        void getRightHandSide(const double& t, double* u, double* fu) override;
};