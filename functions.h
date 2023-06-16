#define _USE_MATH_DEFINES

#include <iostream>
#include <armadillo>
#include <cmath>

#include "system.h"

vec getOmegaJs(System);

vec getCJs(System);

mat getBPrimeG(System);

cx_mat getBPrimeE(System, cx_mat, double);

cx_rowvec getKPrime(System, int);

void getDiagFuncs(cx_mat, cx_mat*, cx_mat*, cx_mat*, cx_mat*, cx_mat*, System);

double trapezoidInt(vec, vec, int);

double getMarcus(System, double);

