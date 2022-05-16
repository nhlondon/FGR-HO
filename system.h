#define _USE_MATH_DEFINES

#include <iostream>
#include <armadillo>
#include <cmath>
#include <string>

using namespace std;
using namespace arma;

class System 
{
	public:
		System(string);
		double aVal;
		double bVal;
		double delta;
		double m;
		int f;
		int n;
		double omegac;
		double omegas;
		double etaRed;
		double eta;
		double temp;
		double beta;

		vec epsilon;
		vec time;
	
		bool integrateC;
		bool drivingForceSweep;
		vec omegaj;
		vec cj;
};
