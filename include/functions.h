#pragma once

/**
* 2D curves
*/
void fold(double* x, double p_u);
void circle(double* x, double p_u);
void lissajous(double *x, double p_u);

/**
* 3D curves
*/
void spherical(double* x, double p_u);
void inclined_circle(double *x, double p_u);
void example(double *x, double p_u);
void torus_helix(double *x, double p_u);