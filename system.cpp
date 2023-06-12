#define _USE_MATH_DEFINES

#include <iostream>
#include <armadillo>
#include <cmath>
#include <string>
#include <iomanip>

#include "system.h"

using namespace std;
using namespace arma;

System::System (string file)
{
	ifstream ip;
	ip.open(file);
	if(ip.fail()) {cout << "Can't open file: "<< file << endl; exit(0); }

	vec epsTemp(2);
	vec timeTemp(2);
	string intC;
	string dForce;
	
	ip >> aVal;
	ip >> bVal;
	ip >> delta;
	ip >> m;
	ip >> f;
	ip >> omegac;
	ip >> etaRed;
	ip >> temp;
	ip >> epsTemp[0] >> epsTemp[1];
	ip >> timeTemp[0] >> timeTemp[1];
	ip >> intC;
	ip >> dForce;

	ip.close();
	eta = etaRed * m * omegac;
	omegas = sqrt((2 * aVal)/m);
	beta = 1.0/(3.16681520371153e-6 * temp);
	epsilon = epsTemp;
	time = timeTemp;
	n = f + 1;
	
	if (intC == "true") integrateC = true;
	else integrateC = false;

	if (dForce == "true") drivingForceSweep = true;
	else drivingForceSweep = false;
	ip.close();
}
