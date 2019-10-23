// Authors: Spencer Behling, Benjamin Fenker, Melissa Anholm - 2012

#ifndef K37Cloud_h
#define K37Cloud_h 1

#include "G4ParticleDefinition.hh" // calls:  #include <CLHEP/Units/SystemOfUnits.h>
#include "G4SystemOfUnits.hh" // calls:  #include <CLHEP/Units/SystemOfUnits.h> and then lets us use unit names for free.
#include "G4PhysicalConstants.hh"
#include "G4IonTable.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

//#include "K37AtomicMessenger.hh"
//class K37AtomicMessenger;

class K37Cloud 
{
public:
	K37Cloud();
	~K37Cloud();

	// For the whole cloud:
	G4ThreeVector GetCloudCenter()      { return cloud_center_;               }
	G4ThreeVector GetInitialCloudSize() { return initial_cloud_size_;         }
	G4ThreeVector GetTemperature()      { return temperature_;                }
	G4ThreeVector GetSailVelocity()     { return sail_velocity_;              }
	G4ThreeVector GetVelocitySigma()    { return velocity_sigma_;             }  // must be already set up with temperature?
	G4double GetFreeExpansionTime()     { return expansion_before_polarized_; }
	G4double GetOP_CycleTime()          { return cycleTime;                   }
	
	void SetCloudCenter(G4ThreeVector center);
	void SetTemperature(G4ThreeVector temp);
	void SetTemperature(G4double temp);
	void SetInitialCloudSize(G4ThreeVector size);
	void SetInitialCloudSize(G4double size);
	void SetSailVelocity(G4ThreeVector vel);

	void SetFreeExpansionTime(G4double time);
	void SetOP_CycleTime(G4double time);

private:
// For the whole cloud:
	G4ThreeVector cloud_center_;
	G4ThreeVector initial_cloud_size_;
	G4ThreeVector velocity_sigma_;
	G4ThreeVector temperature_;
	G4ThreeVector sail_velocity_;
	G4double cycleTime;
	G4double expansion_before_polarized_;
	
public:
	void SetupVelocitySigma(G4ThreeVector temperature);  // We might like this to be private, but our instantiation of this class is itself private, so wev.
	G4double CalcSigma(G4double temperature);            // We might like this to be private, but our instantiation of this class is itself private, so wev.
	void Initialize();                                   // We might like this to be private, but our instantiation of this class is itself private, so wev.
private:	
	G4bool initialize_complete_;
	
};

#endif
