
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
	
	sys.omegaj = getOmegaJs(sys);
	sys.cj = getCJs(sys);
	//cout << sys.omegas << endl << endl;
	//cout << sys.omegaj << endl << endl;
	//cout << sys.cj << endl << endl;
	
	cx_mat tauMat;
	cx_mat tauZ(sys.n, sys.n, fill::zeros);
	cx_mat tauI(sys.n, sys.n, fill::eye);
	tauMat = (join_cols(join_rows(tauZ, tauI), join_rows(-tauI, tauZ)));

	mat aPrime, atest;

	cx_double x;
	x = cx_double(1.0);
	//cout << x << endl;
	aPrime = getAPrime(sys);
	atest = getAtest(sys);

//	cout << aPrime << endl << endl;	
	vec eigVal,zVec(sys.n,fill::zeros);
	mat eigVecsA;
	cx_mat aMat1, aMat;

	eig_sym(eigVal, eigVecsA, aPrime);
	cout << eigVal << endl;
//	cout << eigVecsA << endl;

	//cout << 2.0 * trans(eigVecsA) * atest * eigVecsA;
	double z0 = 0.0;
	double z02 = 0.0;

	for(int i=0; i < sys.n; i++)
	{
		eigVal(i) = sqrt(eigVal(i));
	}
	cx_vec eigValsA((sqrt(2.0/sys.m)) * eigVal,zVec);
	aMat1 = diagmat(eigValsA);
//	aMat = join_cols(join_rows(aMat1, tauZ),join_rows(tauZ,aMat1));
	aMat = join_cols(join_rows(tauZ, aMat1),join_rows(aMat1,tauZ));
/*
	z0 = trace(expmat_sym(-cx_double(sys.beta) * aMat/ cx_double(2.0,0))).real();
	for(int i=0; i < (sys.n); i++)
		{	
			z02 += 2.0*exp(-sys.beta * (real(eigValsA(i)))/2.0 );
		}
	cout << z0 << endl;
	cout << z02 << endl;
*/
//	cout << trace(expmat_sym(-sys.beta * real(aPrime) / 2.0)) << endl;
//	cout << real(aMat) << endl << endl;;
//	cout << real(eigVecsA * aPrime / cx_double(2.0,0) * strans(eigVecsA)) << endl;
//	cout << real(eigVecsA * strans(eigVecsA)) << endl;

	cx_rowvec kPrime, kVec1, kVec;
	cx_mat bPrime, bMat;

	bMat = getBPrime(sys, aMat, sys.epsilon(0));
	kPrime = getKPrime(sys);

//	bMat = strans(eigVecsA) * (bPrime / cx_double(2.0,0.0)) * eigVecsA;
	kVec1 = cx_double(1.0/sqrt(2.0),0.0)* kPrime * eigVecsA;
	kVec = join_rows(kVec1,kVec1);

	//cout << bMat << endl << endl;
	//cout << kPrime << endl << endl;
	//cout << kVec << endl;


	cx_mat sG, sE;
	cx_vec lambdaE;
	int nTime = 2000;
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

	double epsilon = sys.epsilon(0);
	fp << "#time " << "cTau " << "Re[cTau]" << endl;
	for (int i = 0; i <= nTime; i++)
	{
		time(i) = (i * tstep) + sys.time[0];
	//	sG = cx_double(0.0,1.0) * (cx_double(time,0.0) + 
	//		cx_double(0.0,1.0) * ( cx_double(sys.beta,0.0) - cx_double(sys.beta / 2.0 , 0.0))) * aMat;
	//	sE = -cx_double(0.0,1.0) * (cx_double(time,0.0) - cx_double(sys.beta / 2.0,0.0)) * bMat;
	//	lambdaE = strans((-cx_double(0.0,1.0) * (cx_double(time,0.0) - cx_double(0.0,1.0)*cx_double(sys.beta / 2.0,0.0)) 
	//		* strans(kVec)) * inv(tauMat));
			double tauVal = 0.0;
			sG = cx_double(-(sys.beta-tauVal), time(i)) * aMat;
			sE = cx_double(-tauVal, -time(i)) * bMat;
			lambdaE = trans( cx_double(-tauVal, -time(i)) * kVec * inv(tauMat));
			cx_double constants = cx_double(-tauVal, -time(i)) * 
				cx_double((sys.bVal * sys.bVal)/( sys.m * sys.omegas * sys.omegas) - epsilon);
		
		if(time(i) != 0.0){	
		//cout << aMat << endl << endl;
		//cout << bMat << endl << endl;
		cx_mat tauSG, tauSE, tauS, expTauS, sMat;
		//cx_mat dSG(2*(sys.n),2*(sys.n), fill::zeros);  
		//cx_mat dSE(2*(sys.n),2*(sys.n), fill::zeros);  
		//cx_mat dS(2*(sys.n),2*(sys.n), fill::zeros);  
		cx_mat dSG, dSE, dS;
		cx_mat mSE, mSG, mS, invMSG, invMSE, invMS;
		cx_mat chiSE, chiDagSE, phiSE, phiDagSE, psiSE;
		cx_mat chiSG; 
		cx_mat chiDagSG, phiSG, phiDagSG, psiSG;
		cx_mat chiS, chiDagS, phiS, phiDagS, psiS;

		cx_vec eigValSE, eigValSG, eigValS;

		tauSG = tauMat * sG;
		tauSE = tauMat * sE;
//		cout << tauSE << endl;

		eig_gen(eigValSG, mSG, tauSG);
		eig_gen(eigValSE, mSE, tauSE);

		dSG = diagmat(eigValSG);
		dSE = diagmat(eigValSE);
		//cout << dSE << endl;
	//	cout << mSE << endl;
		invMSG = inv(mSG);
		invMSE = inv(mSE);
		
		expTauS = mSG * expmat(dSG) * invMSG * mSE * expmat(dSE) * invMSE;
		tauS = logmat(expTauS);
		eig_gen(eigValS, mS, tauS);
		invMS = inv(mS);
		sMat = inv(tauMat) * tauS;
		dS = diagmat(eigValS);

		//cout <<expmat(dSE)<< endl;

		getDiagFuncs(dSG, &chiSG, &chiDagSG, &phiSG, &phiDagSG, &psiSG, sys);
		getDiagFuncs(dSE, &chiSE, &chiDagSE, &phiSE, &phiDagSE, &psiSE, sys);
		getDiagFuncs(dS, &chiS, &chiDagS, &phiS, &phiDagS, &psiS, sys);

		cx_vec lambda;
		cx_double nu, nu2;
		cx_mat nuMat;

//		lambda = mS * chiDagS * invMS * mSE * chiSE * invMSE * lambdaE;
		lambda = mS * phiDagS * invMS * mSE * phiSE * invMSE * lambdaE;		
		nuMat = cx_double(0.5,0.0) * trans(lambdaE) * tauMat * mSE * psiSE * invMSE * lambdaE - 
			(cx_double(0.5,0.0) * trans(lambda) * tauMat * mS * psiS * invMS * lambda);

		nu = nuMat(0,0);

	//	cx_double cTau;
		cx_mat identity(2* (sys.n), 2*(sys.n), fill::eye);

		cx_mat detMat, expMat;
		cx_double detVal, expVal;
		//cout << expTauS << endl;
		detVal = det(expTauS - identity);
	//	detVal = detMat(0,0);
		expMat = trans(lambda) * tauMat * inv(tauS) * lambda;
		expVal = expMat(0,0);
		cTau(i) = (cx_double(1.0,0.0) / sqrt( cx_double(pow(-1.0, sys.n),0.0) * detVal)) * 
	//	cTau(i) =  sqrt( cx_double(pow(-1.0, sys.n),0.0) * detVal) * 
			exp(nu - cx_double(0.5,0.0)*expVal)
			* exp(constants);
		cTauReal(i) = cTau(i).real();
		}

		else{
		cx_mat tauSG, tauSE, tauS, expTauS, sMat;
		cx_mat identity(2* (sys.n), 2*(sys.n), fill::eye);

		cx_mat detMat, expMat;
		cx_double detVal, expVal;
		tauSG = tauMat * sG;
		expTauS = expmat(tauSG);
		detVal = det(expTauS - identity);
		cTau(i) = (cx_double(1.0,0.0) / sqrt( cx_double(pow(-1.0, sys.n),0.0) * detVal)); 
		cTauReal(i) = cTau(i).real();
		}
		fp << time(i) << " " << cTau(i) << " " << cTauReal(i) << endl;
		//cout << detVal << endl;
		//cout << expMat << endl;
		//cout << nuMat << endl;
	//	cout << cTau << endl << endl;;
}
	fp.close();
	cx_mat identity(2* (sys.n), 2*(sys.n), fill::eye);
	sG = cx_double(-sys.beta, 0.0) * aMat;
	cx_mat tauSG = tauMat * sG;
	cx_double detVal = det(expmat(tauSG) - identity);
	cout << detVal << endl;
	z0 = (cx_double(1.0,0.0) / sqrt( cx_double(pow(-1.0, sys.n),0.0) * detVal)).real(); 
	cout << z0 << endl;
	double integral, rate, logRate;
	integral = trapezoidInt(time, cTauReal, nTime);
	rate = sys.delta * sys.delta * integral / z0;
	logRate = log10(rate);
	
	cout << integral << " " << rate << " " << logRate << endl;
  double n = timer.toc();
	cout << "Time taken: " << timer.toc() <<"s" << endl;
	
return 0;
  }
