# curve-evolution

This is a curve flow solver for curves in plane and space used in my research thesis at the university.

There is no UI currently, it is necessary to recompile the code with new setup each time.

The program solves numerical represantation of equation for curves evolving by the curvature.

* For time integration, it is possible to use either Euler or Runge-Kutta-Merson method
* The equations are solved using **parametric** or also called **direct** approach, so it is possible to
use redistribution of points for better stability
* Currently it is possible to use algorithm without resitribution, with natural redistribution (also called DeTurck trick)
or asymptotically uniform redistribution
* Evolution with additional force term might be added in the future
