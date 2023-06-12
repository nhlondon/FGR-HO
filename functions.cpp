
#define _USE_MATH_DEFINES

#include <iostream>
#include <armadillo>
#include <cmath>

#include "functions.h"

using namespace std;
using namespace arma;

vec getOmegaJs (System sys)
{
	vec omegaj(sys.f);
	for (int j = 0; j < sys.f; j++)
	{
		omegaj[j] = -sys.omegac * log((j+1 - 0.5)/sys.f);
	}
	return omegaj;
}

vec getCJs (System sys)
{
	vec cj(sys.f);
	for (int j = 0; j < sys.f; j++)
	{
		cj[j] = sys.omegaj[j] * sqrt((2.0 * sys.eta * sys.m * sys.omegac)/(sys.f * M_PI));
	}
	return cj;
}

mat getBPrimeG(System sys)
{
	//Generate Hessian Matrix of H_g
	mat hes(sys.n,sys.n,fill::zeros);
	for(int j=0; j < sys.f; j++)
	{
		hes(0,0) += (sys.cj[j] * sys.cj[j])/(sys.m * sys.omegaj[j] * sys.omegaj[j]);
		hes(0,j+1) = -sys.cj[j];
		hes(j+1,0) = hes(0,j+1);
		hes(j+1,j+1) = sys.m * sys.omegaj[j] * sys.omegaj[j];
	}
	hes(0,0) += sys.m * sys.omegas * sys.omegas;
	return hes;
}
cx_mat getBPrimeE(System sys, cx_mat bPrimeG, double epsilon)
{
	//Generate Hessian Matrix of H_e - identential to H_g
	cx_mat bPrimeE = bPrimeG; 
	return bPrimeE;
}

cx_rowvec getKPrime(System sys, int state)
{
	//Generate normal mode linear displacements for chosen state
	cx_rowvec kPrime(sys.n);
	kPrime.zeros();
	if (state == 1) kPrime(0) = cx_double(sys.bVal);
	else kPrime(0) = cx_double(-sys.bVal );

	return kPrime; 
}

void getDiagFuncs(cx_mat dS, cx_mat* chi, cx_mat* chiDag, cx_mat* phi, cx_mat* phiDag, cx_mat* psi, System sys)
{
	//Generate diagonal matrices for Bogoliubov transformation
	cx_mat identity(2*(sys.n), 2*(sys.n), fill::eye);
	cx_mat invDS = inv(dS);
	cx_mat expDS = expmat(dS);
	cx_mat expNegDS = expmat(-dS);
	cx_mat sinhDS = (expDS - expNegDS) / 2.0;
	*chi = (expDS - identity) * invDS;
	*chiDag = dS * inv(expDS - identity);
	*phi = (identity - expNegDS) * invDS;
	*phiDag = dS * inv(identity - expNegDS);
	*psi = (sinhDS - dS) * inv(powmat(dS,2));
	
	return;
}

double trapezoidInt(vec x,vec func,int numPoints)
{
	//Integration using trapezoid rule
	double integrand = 0.0;
	for (int i = 1; i < numPoints; i++)
	{
		integrand += (func(i-1) + func(i)) * abs(x(i) - x(i-1)) / 2.0;
	}
	return integrand;
}

