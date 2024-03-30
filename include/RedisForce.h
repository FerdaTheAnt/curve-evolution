#pragma once

#include "CurveFlow.h"

class RedisForce : public CurveFlow
{
public:
    RedisForce(int dimension, int n);
    ~RedisForce();
    void getRightHandSide(const double& t, double* u, double* fu) override;
protected:
    void getLength();
    void getPartialLength(const double* u);
    void getCurvature(const double* u);
    void getForce(const double* u);
    void getOmega();
    void getSeries();
    void getSum();
    void getAlpha();
private:
    double length;
    double omega;
    double* alpha;
    double* d;
    double* k;

    double* f;
};
