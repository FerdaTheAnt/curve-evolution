#pragma once

#include "CurveFlow.h"

class LengthPreserving : public CurveFlow
{
    public:
        LengthPreserving(int dimension, int n);
        virtual void getRightHandSide(const double& t, double* u, double* fu) override;
        virtual ~LengthPreserving();
    private:
        double* k;
        double F;
        void getForce();
        void getCurvature(const double *u);
};
