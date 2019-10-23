// Authors: Spencer Behling, Benjamin Fenker, Melissa Anholm - 2012
//
/// \class K37CloudSetup  K37CloudSetup.hh K37CloudSetup
/// \brief This is a class to change the cloud size based on temperature
/// \author Spencer Behling
/// \date Jan. 5 2012
/// \since version 1.0
///
// 24 September 2019 - Major revision to class structure and functionality by MJA.  
// K37CloudSetup is separated into three classes now -- K37AtomicSetup (this), 
// K37Cloud, and K37SublevelPopulations.
// K37AtomicSetup is just here to consolidate all the atomic stuff into one neat 
// package, for better backward compatibility.  
// K37Cloud will deal with the spatial properties of the cloud itself.  
// eg, size, position, temperature.  
// K37SublevelPopulations will deal with the spin properties of the K37 atoms. 
// eg, polarization, alignment, octopole moment, etc.  
// In particular, K37SublevelPopulations will deal with the messy problems involved
// in attempting to predict higher multipole moments based on our knowledge of the 
// lower moments.

#ifndef K37AtomicSetup_h
#define K37AtomicSetup_h 1

#include "G4ParticleDefinition.hh" // calls:  #include <CLHEP/Units/SystemOfUnits.h>
//#include "G4SystemOfUnits.hh"      // calls:  #include <CLHEP/Units/SystemOfUnits.h> and then lets us use unit names for free.
#include "G4PhysicalConstants.hh"
#include "G4IonTable.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "K37Config.hh"
#include "K37Cloud.hh"
#include "K37SublevelPopulations.hh"

#ifndef SIMPLE_MC
	#include "K37AtomicMessenger.hh"
	class K37AtomicMessenger;
#endif

class K37Cloud;
class K37SublevelPopulations;

class K37AtomicSetup 
{
public:
	K37AtomicSetup();
	~K37AtomicSetup();
	//
	
public:
	K37Cloud * GetCloud()              { return the_cloud; } // Is this ok for K37?
	K37SublevelPopulations * GetPops() { return the_pops;  } // Is this OK for K37?  ... I don't think it's even ok for Holstein.
	
private:
	K37Cloud               * the_cloud;
	K37SublevelPopulations * the_pops;
	#ifndef SIMPLE_MC
		K37AtomicMessenger     * AtomicMessenger;   
	#endif

public:
	// cloud:
	G4ThreeVector GetCloudCenter()      { return the_cloud->GetCloudCenter();       }
	G4ThreeVector GetInitialCloudSize() { return the_cloud->GetInitialCloudSize();  }
	G4ThreeVector GetTemperature()      { return the_cloud->GetTemperature();       }
	G4ThreeVector GetSailVelocity()     { return the_cloud->GetSailVelocity();      }
	G4ThreeVector GetVelocitySigma()    { return the_cloud->GetVelocitySigma();     }  // must be already set up with temperature?
	G4double GetFreeExpansionTime()     { return the_cloud->GetFreeExpansionTime(); }
	G4double GetOP_CycleTime()          { return the_cloud->GetOP_CycleTime();      }
	// pops:
	double GetPolarization()            { return the_pops->GetPolarization(); }
	double GetAlignment()               { return the_pops->GetAlignment();    }
	void SetPolarization(double pol)    { the_pops->AdjustPolarization(pol);  }
	
	/*
	// Setup Multipole:
	void Setup_FromPolarizationOnly(double pol)                                       { the_pops->Setup_FromPolarizationOnly(pol);                        }
	void Setup_FromPolarizationAlignment(double pol, double ali)                      { the_pops->Setup_FromPolarizationAlignment(pol, ali);              }
	void Setup_FromPolarizationAlignmentOctopole(double pol, double ali, double oct)  { the_pops->Setup_FromPolarizationAlignmentOctopole(pol, ali, oct); } // this method is stupid in that it mixes conventions.
	void Setup_FromDipoleQuadrupoleOctopole(double dip, double quad, double oct)      { the_pops->Setup_FromDipoleQuadrupoleOctopole(dip, quad, oct);     }
	*/
	
	void SetCloudCenter(G4ThreeVector center)    { the_cloud->SetCloudCenter(center);     }
	void SetTemperature(G4ThreeVector temp)      { the_cloud->SetTemperature(temp);       }
	void SetTemperature(G4double temp)           { the_cloud->SetTemperature(temp);       }
	void SetInitialCloudSize(G4ThreeVector size) { the_cloud->SetInitialCloudSize(size);  }
	void SetInitialCloudSize(G4double size)      { the_cloud->SetInitialCloudSize(size);  }
	void SetSailVelocity(G4ThreeVector vel)      { the_cloud->SetSailVelocity(vel);       }
	void SetFreeExpansionTime(G4double time)     { the_cloud->SetFreeExpansionTime(time); }
	void SetOP_CycleTime(G4double time)          { the_cloud->SetOP_CycleTime(time);      }

	void SetMatchedRunsetLetter(G4String newRunsetLetter);
	G4String GetMatchedRunsetLetter();

private:
	G4String MatchedRunsetLetter;
		
	void SetupVelocitySigma(G4ThreeVector temperature) { the_cloud->SetupVelocitySigma(temperature); }
	G4double CalcSigma(G4double temperature)           { return the_cloud->CalcSigma(temperature); }
	G4bool initialize_complete_;
	void Initialize();
	
//private:
// these are private to prevent the user from calling them directly.
// Instead, messenger classes will use SetPolarizationAlignment(...)
// and SetPolarizationOnly(...).
//	void SetPolarization(double pol) { the_pops->SetPolarization(pol); }
//	void SetAlignment(double ali)    { the_pops->SetAlignment(ali);    }
//	void SetOctopole(double oct)     { the_pops->SetOctopole(oct);     }
// Utility classes...
//	void MakeAlignmentFromPolarization()         { the_pops->MakeAlignmentFromPolarization();         }
//	void MakeOctopoleFromPolarizationAlignment() { the_pops->MakeOctopoleFromPolarizationAlignment(); }
	
};

#endif
