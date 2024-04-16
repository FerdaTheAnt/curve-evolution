#pragma once

#include "CurveFlow.h"

class RedistributionSpace : public CurveFlow
{
public:
    RedistributionSpace(int dimension, int n);
    ~RedistributionSpace();
    virtual void getRightHandSide(const double& t, double* u, double* fu) override;
protected:
    void getCurvature(const double* u);
    void getOmega();
    void getSeries();
    void getSum();
    void getAlpha();
private:
    double omega;
    double* alpha;
    double* k;
};
