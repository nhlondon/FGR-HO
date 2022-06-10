#define _USE_MATH_DEFINES

#include <iostream>
#include <armadillo>
#include <cmath>

#include "system.h"

vec getOmegaJs(System);

vec getCJs(System);

mat getAPrime(System);

cx_mat getBPrime(System, cx_mat, double);

cx_rowvec getKPrime(System, int);

void getDiagFuncs(cx_mat, cx_mat*, cx_mat*, cx_mat*, cx_mat*, cx_mat*, System);

double trapezoidInt(vec, vec, int);

