// Authors: Spencer Behling, Benjamin Fenker, Melissa Anholm - 2013

/*-----------------------------------------------------------------------------
  This is the re-write of my original version of the fermi-function. Originally
  I was using the gsl version of the complex gamma function. I was able to
  implement my own version that now appears in K37ComplexGammaFunction.cc. My
  version agrees with the Mathematica version and is much easier to use than
  The gsl version. Also using my own version has allowed me to eliminate a
  dependency that has caused problems in the past with linking.
  -----------------------------------------------------------------------------*/

#include "K37FermiFunction.hh"
#include <cmath>
#include "K37ComplexGammaFunction.hh"

K37FermiFunction::K37FermiFunction(double zChoice, double aChoice)
  :CGF(0) 
{
	CGF = new K37ComplexGammaFunction();
	A = aChoice;
	Z = zChoice;
	R = 2.5896E-3 *(1.2*pow(A, (1./3.)));
	alpha = 1./137.036;
	gamma = sqrt(1.-(alpha*alpha*Z*Z));
	twoGammaPlusTwo = 2.*(gamma + 1.);
	twoGammaMinusTwo = 2.*(gamma - 1.);
//	massOfElectron = 0.510998;  // in MeV.  
	massOfElectron = 0.5109989461;  // in MeV.  Slightly more precise value.
}

K37FermiFunction::~K37FermiFunction() 
{
	delete CGF;
}

G4double K37FermiFunction::getVFF(G4double kineticEnergyChoice) 
{
	T = kineticEnergyChoice;
	W = (T/massOfElectron) + 1.;
	p = sqrt((W*W) - 1.);
	alphaZWoverP = (alpha*Z*W)/p;
	exponentResult = exp(M_PI*alphaZWoverP);
	twoPR = 2.*p*R;
	powerResult = pow(twoPR, twoGammaMinusTwo);
	ratio = (CGF->absSquaredComplexGamma('g', gamma, alphaZWoverP))/
	  (CGF -> squaredRealGamma((2.*gamma +1.)));
	fermi = twoGammaPlusTwo*ratio*powerResult*exponentResult;
	
	return fermi;
}
