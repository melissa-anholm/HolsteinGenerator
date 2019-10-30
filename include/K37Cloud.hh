// Authors: Spencer Behling, Benjamin Fenker, Melissa Anholm - 2012

#ifndef K37Cloud_h
#define K37Cloud_h 1

#include "G4ParticleDefinition.hh" // calls:  #include <CLHEP/Units/SystemOfUnits.h>
#include "G4SystemOfUnits.hh" // calls:  #include <CLHEP/Units/SystemOfUnits.h> and then lets us use unit names for free.
#include "G4ThreeVector.hh"
#include "globals.hh"


class K37Cloud 
{
public:
	K37Cloud();
	~K37Cloud();

public:
	// For the whole cloud:
//	G4ThreeVector GetCloudCenter()        { return initial_position;           }
	G4ThreeVector GetInitialCloudPosition() { return initial_position;           }
	G4ThreeVector GetFinalCloudPosition()   { return final_position;             }
	G4ThreeVector GetInitialCloudSize()     { return initial_cloud_size;         }
	G4ThreeVector GetFinalCloudSize()       { return final_cloud_size;           }
	
	G4double GetFreeExpansionTime()         { return expansion_before_polarized; }
	G4double GetOP_CycleTime()              { return cycleTime;                  }

//	G4ThreeVector GetTemperature()      { return temperature;                }
//	G4ThreeVector GetSailVelocity()     { return sail_velocity;              }
//	G4ThreeVector GetVelocitySigma()    { return velocity_sigma;             }  // must be already set up with temperature?
	
//	void SetCloudCenter(G4ThreeVector center);
	void SetInitialCloudPosition(G4ThreeVector center);
	void SetFinalCloudPosition(G4ThreeVector center);
	
	void SetInitialCloudSize(G4ThreeVector size);
	void SetInitialCloudSize(G4double size);
	
	void SetFinalCloudSize(G4ThreeVector size);
	void SetFinalCloudSize(G4double size);
	
	void SetFreeExpansionTime(G4double time);
	void SetOP_CycleTime(G4double time);

	// Given that I have an initial cloud size and initial cloud position, 
	// I might set up the cloud's evolution by Temperature and Sail Velocity, 
	// or equivalently I might use Final Size and Final Position.  
	// It would be a mistake to try to code both in.
	// Ben used temperature and sail velocity, but I don't really like it.
	// I might change it later.
	/*
	void SetTemperature(G4ThreeVector temp);
	void SetTemperature(G4double temp);
	void SetSailVelocity(G4ThreeVector vel);
	*/


	
private:
	void set_up_sail_velocity();  // initial_position, final_position, and cycleTime must already be set up.
//	void set_up_temperature();    // must have initial_size, final_size, and cycleTime already set up.
	
	
public:  
	void PrintCloud();

private:
// For the whole cloud:
	G4ThreeVector initial_position;
	G4ThreeVector final_position; //
	
	G4ThreeVector initial_cloud_size;
	G4ThreeVector final_cloud_size;
	
//	G4ThreeVector temperature;
//	G4ThreeVector velocity_sigma;
	G4double cycleTime;
	G4double expansion_before_polarized;
	
//	G4ThreeVector sail_velocity;  // legacy.
	
public:
	void SetupVelocitySigma(G4ThreeVector temperature);  // We might like this to be private, but our instantiation of this class is itself private, so wev.
	G4double CalcSigma(G4double temperature);            // We might like this to be private, but our instantiation of this class is itself private, so wev.
//	void Initialize();                                   // We might like this to be private, but our instantiation of this class is itself private, so wev.
//private:	
//	G4bool initialize_complete_;
	
};

#endif
