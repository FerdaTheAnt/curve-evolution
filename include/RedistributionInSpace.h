#pragma once

#include "CurveFlow.h"

class RedistributionSpace : public CurveFlow
{
public:
    RedistributionSpace(int dimension, int n);
    ~RedistributionSpace();
    void getRightHandSide(const double& t, double* u, double* fu) override;
protected:
    void getLength();
    void getPartialLength(const double* u);
    void getCurvature(const double* u);
    void getOmega();
    void getSeries();
    void getSum();
    void getAlpha();
private:
    int STEPS = 1;
    double length;
    double omega;
    double* alpha;
    double* d;
    double* k;
};
