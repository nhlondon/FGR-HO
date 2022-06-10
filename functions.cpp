
#define _USE_MATH_DEFINES

#include <iostream>
#include <armadillo>
#include <cmath>

#include "functions.h"

using namespace std;
using namespace arma;
//vec getOmegaJs (double omegac, int f)
vec getOmegaJs (System sys)
{
	vec omegaj(sys.f);
	for (int j = 0; j < sys.f; j++)
	{
		omegaj[j] = -sys.omegac * log((j+1 - 0.5)/sys.f);
	}
	return omegaj;
}

//vec getCJs (vec omegaj, double omegac, int f, double eta, double m)
vec getCJs (System sys)
{
	vec cj(sys.f);
	for (int j = 0; j < sys.f; j++)
	{
		cj[j] = sys.omegaj[j] * sqrt((2.0 * sys.eta * sys.m * sys.omegac)/(sys.f * M_PI));
	}
	return cj;
}
/*
cx_mat getAPrime(System sys)
{
	/* Generate A' matrix of the form:
	B A
	A B 
	*/
	/*
	cx_mat aPrimeA(sys.n ,sys.n), aPrimeB(sys.n,sys.n);
	cx_mat aPrimeAB, aPrimeBA, aPrime;

	aPrimeA.zeros();
	aPrimeB.zeros();
	for(int j=0; j < sys.f; j++)
	{
		aPrimeA(0,0) += cx_double(1.0/(4.0 * sys.m * sys.m * sys.omegas) * 
			((sys.cj[j] * sys.cj[j]) / (sys.omegaj[j] * sys.omegaj[j])), 0.0);
		aPrimeA(0,j+1) = cx_double(-sys.cj[j] / (4.0 * sys.m * sqrt( sys.omegas * sys.omegaj[j])), 0.0);
		aPrimeA(j+1,0) = aPrimeA(0,j+1);
		aPrimeB(0,j+1) = aPrimeA(0,j+1);
		aPrimeB(j+1,0) = aPrimeA(0,j+1);
		aPrimeB(j+1,j+1) = cx_double(sys.omegaj[j] / 2.0, 0.0);
	}
	aPrimeB(0,0) = aPrimeA(0,0) + cx_double(sys.omegas / 2.0, 0.0);
	aPrimeAB = join_rows(aPrimeA,aPrimeB);
	aPrimeBA = join_rows(aPrimeB,aPrimeA);

	cx_mat eig1, eig2;
	vec eigVal1, eigVal2;

	eig_sym(eigVal1,eig1,aPrimeA+aPrimeB);
	cout << eigVal1 << endl;
	eig_sym(eigVal2,eig2,aPrimeB-aPrimeA);
	cout << eigVal2 << endl;
	aPrime = join_cols(aPrimeBA,aPrimeAB);
	return aPrime;
}
*/

mat getAPrime(System sys)
{
	//Generate Hessian Matrix
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
cx_mat getBPrime(System sys, cx_mat aPrime, double epsilon)
{
//	cx_mat identity(2*sys.n,2*sys.n,fill::eye);
	cx_mat bPrime = aPrime; //+(cx_double((2.0 * sys.bVal * sys.bVal) / (sys.m * sys.omegas * sys.omegas) - epsilon)
	//	* identity);
//	for(int j = 0; j < sys.n; j++)
//	{
//	bPrime(j, sys.n+j) -= cx_double(((2.0 * sys.bVal * sys.bVal) / ( sys.m * sys.omegas * sys.omegas ) - epsilon) /
//		(2.0 * sys.n));
//	bPrime(sys.n+j, j) += cx_double(((2.0 * sys.bVal * sys.bVal) / ( sys.m * sys.omegas * sys.omegas ) - epsilon) /
//		(2.0 * sys.n));
	//bPrime(j, j) -= cx_double((2.0 * sys.bVal * sys.bVal) / ( sys.m * sys.omegas * sys.omegas ) - epsilon);
	//bPrime(j+sys.n,j+sys.n) += cx_double((2.0 * sys.bVal * sys.bVal) / ( sys.m * sys.omegas * sys.omegas ) - epsilon);
	//}
	//cout << cx_double((2.0 * sys.bVal * sys.bVal) / ( sys.m * sys.omegas * sys.omegas ) - epsilon) << endl << endl;
	return bPrime;
}

cx_rowvec getKPrime(System sys, int state)
{
	cx_rowvec kPrime(sys.n);
	kPrime.zeros();
	//if (state == 1) kPrime(0) = cx_double(sys.bVal / sqrt(sys.m * sys.omegas));
	//else kPrime(0) = cx_double(-sys.bVal / sqrt(sys.m * sys.omegas));
	if (state == 1) kPrime(0) = cx_double(sys.bVal);
	else kPrime(0) = cx_double(-sys.bVal );

	return kPrime; 
}

void getDiagFuncs(cx_mat dS, cx_mat* chi, cx_mat* chiDag, cx_mat* phi, cx_mat* phiDag, cx_mat* psi, System sys)
{
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
	//cout << psi << endl;
	
	return;
}

double trapezoidInt(vec x,vec func,int numPoints)
{
	double integrand = 0.0;
	for (int i = 1; i < numPoints; i++)
	{
		integrand += (func(i-1) + func(i)) * abs(x(i) - x(i-1)) / 2.0;
	}
	return integrand;
}


