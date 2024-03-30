#pragma once

#include "ODEProblem.h"
#include "ODESolver.h"

#include <cmath>
#include <algorithm>
#include <iostream>

bool solve(
    const double initialTime,
    const double finalTime,
    const double timeStep,
    const double integrationTimeStep,
    ODEProblem *problem,
    ODESolver * solver,
    double* u
)
{
    solver->setup(problem->getDegreesOfFreedom());
    const int timeStepsCount = std::ceil(std::max( 0.0, finalTime - initialTime ) / timeStep);
    double time(initialTime);
    problem->writeSolution(time, 0, u);
    for(int k = 1; k<timeStepsCount; k++)
    {
        std::cout << "Solving time step " << k << " / " << timeStepsCount << " => " <<
        (double) k / (double) timeStepsCount * 100.0 << "% " << std::endl;
        double currentTimeStep = std::min(timeStep, finalTime - time);
        if(!solver->solve(
                integrationTimeStep,
                time + currentTimeStep,  // stopTime
                &time,
                problem,
                u
            )
        )
            return false;
        problem->writeSolution( time, k, u );
   }
   std::cout << "Done." << std::endl;
   return true;
}
