
#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <string>
#include <time.h>
#include "functions.h"

using namespace std;
using namespace arma;

int main()
{
	wall_clock timer;
	timer.tic();
	
	System sys("system.txt");

	//Get bath frequencies and coupling	
	sys.omegaj = getOmegaJs(sys);
	sys.cj = getCJs(sys);
	
	//Generate tau rotation matrix
	cx_mat tauMat;
	cx_mat tauZ(sys.n, sys.n, fill::zeros);
	cx_mat tauI(sys.n, sys.n, fill::eye);
	tauMat = (join_cols(join_rows(tauZ, tauI), join_rows(-tauI, tauZ)));

	//H_g hessian
	mat bPrimeG;
	bPrimeG = getBPrimeG(sys);

	//Diagonalize hessian
	vec eigVal,zVec(sys.n,fill::zeros);
	mat eigVecsG;
	cx_mat bMatGSingle, bMatG;

	eig_sym(eigVal, eigVecsG, bPrimeG);
	cout << eigVal << endl;

	for(int i=0; i < sys.n; i++)
	{
		eigVal(i) = sqrt(eigVal(i));
	}
	cx_vec eigValsG((sqrt(1.0/sys.m)) * eigVal,zVec); //changed from root(2/m)
	bMatGSingle = diagmat(eigValsG);
	//Full B matrix for H_g
	bMatG = join_cols(join_rows(tauZ, bMatGSingle),join_rows(bMatGSingle,tauZ));
	
	//B matrix for H_e
	cx_mat bMatE;
	bMatE = getBPrimeE(sys, bMatG, sys.epsilon(0));
	
	//k matrix for both states in standard coordinates
	cx_rowvec kPrimeG, kPrimeE, kVec1E, kVec1G, kVecG, kVecE;
	kPrimeG = getKPrime(sys, 1);
	kPrimeE = getKPrime(sys, 2);

	//Convert to normal mode basis
	kVec1G = cx_double(1.0/sqrt(2.0),0.0)* kPrimeG * eigVecsG;
	for( int j = 0; j < sys.n; j++)
	{ 
		kVec1G(j) = kVec1G(j)/ sqrt(sys.m * eigValsG(j));
	}
	kVecG = join_rows(kVec1G,kVec1G);
	kVec1E = cx_double(1.0/sqrt(2.0),0.0)* kPrimeE * eigVecsG;
	for( int j = 0; j < sys.n; j++)
	{ 
		kVec1E(j) = kVec1E(j)/ sqrt(sys.m * eigValsG(j));
	}
	kVecE = join_rows(kVec1E,kVec1E);



	cx_mat sG, sE;
	cx_vec lambdaG, lambdaE, lambda;
	cx_double nu;
	int nTime = 4000;
	double tstep = (sys.time[1] - sys.time[0])/ nTime;
	vec time(nTime+1);
	cx_vec cTau(nTime+1);
	vec cTauReal(nTime+1);
	stringstream nTimeStr;
	nTimeStr << nTime;
	string fileName;
	fileName = "correlation.xvg";
	ofstream fp;
	fp.open(fileName);
	if(fp.fail()) {cout << "Can't open file" << endl; exit(0); }
	ofstream fp2;
	fp2.open("drivingForce.xvg");
	if(fp2.fail()) {cout << "Can't open file" << endl; exit(0); }

	int nDriving;
	double dFStep = (sys.epsilon[1] - sys.epsilon[0])/nDriving;

	if (sys.drivingForceSweep)
	{
		 nDriving = 50;
		ofstream fp2;
		fp2.open("drivingForce.xvg");
		if(fp2.fail()) {cout << "Can't open file" << endl; exit(0); }
		dFStep = (sys.epsilon[1] - sys.epsilon[0])/nDriving;
	}
	else 
	{
		nDriving = 0;
		dFStep  = 0;
	}
	vec epsilon(nDriving+1);

	fp2 << "epsilon  FGR rate  logFGR  Marcus rate   logMarcus" << endl;
	for (int j =0; j <= nDriving; j++)
	{
		epsilon(j) = (j * dFStep) +sys.epsilon[0];

		//double epsilon = sys.epsilon(0);
		fp << "#time " << "cTau " << "Re[cTau]" << endl;
		for (int i = 0; i <= nTime; i++)
		{
				time(i) = (i * tstep) + sys.time[0];
				double tauVal = 0.1;
				
				//Convert hamiltonians to S, lambda, nu form	
				sG = cx_double(-(sys.beta-tauVal), time(i)) * bMatG;
				sE = cx_double(-tauVal, -time(i)) * bMatE;
				lambdaG = strans( cx_double(-(sys.beta - tauVal), time(i)) * kVecG * inv(tauMat));
				lambdaE = strans( cx_double(-tauVal, -time(i)) * kVecE * inv(tauMat));
				cx_double constants = cx_double(-(sys.beta-tauVal), time(i)) * 
					cx_double(epsilon(j));

				cx_mat tauSG, tauSE, tauS, expTauS, sMat;
				cx_mat dSG, dSE, dS;
				cx_mat mSE, mSG, mS, invMSG, invMSE, invMS;
				cx_mat chiSE, chiDagSE, phiSE, phiDagSE, psiSE;
				cx_mat chiSG; 
				cx_mat chiDagSG, phiSG, phiDagSG, psiSG;
				cx_mat chiS, chiDagS, phiS, phiDagS, psiS;

				cx_vec eigValSE, eigValSG, eigValS;

				//Boboliubov transformation
				tauSG = tauMat * sG;
				tauSE = tauMat * sE;

				eig_gen(eigValSG, mSG, tauSG);
				eig_gen(eigValSE, mSE, tauSE);

				dSG = diagmat(eigValSG);
				dSE = diagmat(eigValSE);
				
				invMSG = inv(mSG);
				invMSE = inv(mSE);
			
				expTauS = mSG * expmat(dSG) * invMSG * mSE * expmat(dSE) * invMSE;
				
				tauS = logmat(expTauS);
				eig_gen(eigValS, mS, tauS);
				invMS = inv(mS);
				sMat = inv(tauMat) * tauS;
				dS = diagmat(eigValS);
				
				getDiagFuncs(dSG, &chiSG, &chiDagSG, &phiSG, &phiDagSG, &psiSG, sys);
				getDiagFuncs(dSE, &chiSE, &chiDagSE, &phiSE, &phiDagSE, &psiSE, sys);
				getDiagFuncs(dS, &chiS, &chiDagS, &phiS, &phiDagS, &psiS, sys);

				cx_mat nuMat;

				lambda = mS * chiDagS * invMS * mSG * chiSG * invMSG * lambdaG +
					mS * phiDagS * invMS * mSE * phiSE * invMSE * lambdaE;		
				nuMat = cx_double(0.5,0.0) * strans(lambdaG) * tauMat * mSG * psiSG * invMSG * lambdaG + 
					cx_double(0.5,0.0) * strans(lambdaE) * tauMat * mSE * psiSE * invMSE * lambdaE - 
					(cx_double(0.5,0.0) * strans(lambda) * tauMat * mS * psiS * invMS * lambda) +
					(cx_double(0.5,0.0) * strans(lambdaG) * tauMat * mSG * chiSG * invMSG * mSE * chiSE *
						invMSE * lambdaE);

				nu = nuMat(0,0);
				
				cx_mat identity(2* (sys.n), 2*(sys.n), fill::eye);

				cx_mat detMat, expMat;
				cx_double detVal, expVal;
				detVal = det(expTauS - identity);
				expMat = strans(lambda) * tauMat * inv(tauS) * lambda;
				expVal = expMat(0,0);
			
				//Correlation function
				cTau(i) = (cx_double(1.0,0.0) / sqrt( cx_double(pow(-1.0, sys.n),0.0) * detVal)) * 
					exp(nu - cx_double(0.5,0.0)*expVal)
					* exp(constants);
				cTauReal(i) = cTau(i).real();

				fp << time(i) << " " << cTau(i) << " " << cTauReal(i) << " " << cTau(i).imag() << endl;
			}
		
		fp.close();
	
		//Reaction partition function
		cx_mat identity(2* (sys.n), 2*(sys.n), fill::eye);
		sG = cx_double(-sys.beta, 0.0) * bMatG;
			lambdaG = strans( cx_double(-(sys.beta), 0.0) * kVecG * inv(tauMat));
		cx_mat tauSG = tauMat * sG;
		cx_double detVal = det(expmat(tauSG) - identity);
		cx_mat	expMat = strans(lambdaG) * tauMat * inv(tauSG) * lambdaG;
		cx_double	expVal = expMat(0,0);
		double z0 = (cx_double(1.0,0.0) / sqrt( cx_double(pow(-1.0, sys.n),0.0) * detVal) * 
				exp((-sys.beta*epsilon(j))-cx_double(0.5,0.0)*expVal)).real();
		cout << z0 << endl;

		//Integrate correlation function
		double integral, rate, logRate, marcus, logMarcus;
		integral = trapezoidInt(time, cTauReal, nTime);
		rate = sys.delta * sys.delta * integral / z0;
		logRate = log10(rate);
		marcus = getMarcus(sys, epsilon(j));
		logMarcus = log10(marcus);	

		cout << epsilon(j) << endl;
		cout << integral  << " " << integral/z0 << " " << rate << " " << logRate << 
			" " << marcus << " " << logMarcus << endl << endl;
		fp2 << epsilon(j) << " " << rate << " " << logRate << " " << marcus << " " << logMarcus << endl;
	}

  double n = timer.toc();
	cout << "Time taken: " << timer.toc() <<"s" << endl;
	
	return 0;
}
