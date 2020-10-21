// Authors: Spencer Behling and Benjamin Fenker 2013

#ifndef K37FermiFunction_h
#define K37FermiFunction_h 1

//#include "K37FermiFunction.hh"
#include "globals.hh"

class K37ComplexGammaFunction;

class K37FermiFunction 
{
public:
	K37FermiFunction(double zChoice = -18., double aChoice = 37.);
	~K37FermiFunction();
	// acronym stands for get (V)alue of (F)ermi (F)unction takes kinetic energy
	// as its argument
	G4double getVFF(G4double kineticEnergyChoice);
	// Kinetic energy should be in units of ???
	// in Ben/Spencer's code it takes units of MeV.  Unless I changed that.

private:
	G4double A;                 // mass number of isotope
	G4double R;                 // Radius of the nucleus
	G4double Z;                 // Charge of isotope negative for beta minus decay
	G4double alpha;
	G4double gamma;
	G4double twoGammaPlusTwo;
	G4double twoGammaMinusTwo;
	G4double massOfElectron;
	G4double W;
	G4double p;
	G4double alphaZWoverP;
	G4double exponentResult;
	G4double twoPR;
	G4double powerResult;
	G4double ratio;
	G4double fermi;
	G4double T;
	
	K37ComplexGammaFunction* CGF;
};


#endif
