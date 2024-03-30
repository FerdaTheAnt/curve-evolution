#pragma once

class ODEProblem
{
    public:
        virtual int getDegreesOfFreedom() = 0;
        virtual void getRightHandSide(const double &t, double *u, double *fu) = 0;
        virtual bool writeSolution(const double &t, int step, const double *u) = 0;
};