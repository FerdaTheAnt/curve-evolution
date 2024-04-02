#include <chrono>

#include "CurveFlow.h"
#include "CurveFlowForce.h"
#include "DeTurck.h"
#include "AreaPreserving.h"
#include "RedistributionInSpace.h"
#include "RedisForce.h"
#include "Euler.h"
#include "Merson.h"
#include "curve-solve.h"
#include "functions.h"

const int dim = 2;
const int n = 800;
const double initialTime = 0.0;
const double finalTime = 0.3;
const double timeStep = 0.01;
const double integrationTimeStep = 0.45/(n*n);

int main(int argc, char **argv)
{
    auto start = std::chrono::steady_clock::now();

    AreaPreserving fdm(dim, n);
    //Euler solver;
    Merson solver;
    solver.setAdaptivity(10e-6);
    double *u = new double[(n+2)*dim]();

    fdm.setInitialCondition(u, fold);

    if( ! solve(
        initialTime,
        finalTime,
        timeStep,
        integrationTimeStep,
        &fdm,
        &solver,
        u))
    {
        delete[] u;

        auto end = std::chrono::steady_clock::now();
        auto diff = end - start;
        std::cout << "Computation time: " << std::chrono::duration <double> (diff).count() << " s" << std::endl;

        return EXIT_FAILURE;
    }
    delete[] u;

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << "Computation time: " << std::chrono::duration <double> (diff).count() << " s" << std::endl;

    return EXIT_SUCCESS;
}
