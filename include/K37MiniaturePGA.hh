// Authors: Spencer Behling, Benjamin Fenker, Melissa Anholm - 2013

// PGA asks the atomic setup where the decay originates,
//     decides if it's a 2% event,
//     determines the recoil charge for this event, 
//     asks the atomic setup for the polarization and alignment,
//     calls one of the event generators to 'create' the event, 
//     sets the beta 'vertex', and if needed the recoil, SOE, and gamma vertices.


#ifndef K37MiniaturePGA_h
#define K37MiniaturePGA_h 1

#include <vector>

#include "globals.hh"

//#include "G4VUserPrimaryGeneratorAction.hh"
#include "K37AtomicSetup.hh"
#include "Holstein52Generator.hh"
//#include "K37EventGenerator.hh"  
//#include "K37EventGeneratorNoRecoilOrderEffects.hh"  // Do I need you?  ... Yes, apparently I do.
//#include "K37PrimaryGeneratorMessenger.hh"


//class K37PrimaryGeneratorMessenger;  // somehow it's not good enough included above.  
//class G4ParticleDefinition;
class G4PrimaryVertex;

using std::vector;

// it's a gun.
//class K37PrimaryGeneratorAction //: public G4VUserPrimaryGeneratorAction 
class K37MiniaturePGA
{
public:
//	K37PrimaryGeneratorAction(K37EventGenerator*, K37EventGeneratorNoRecoilOrderEffects*);
//	K37PrimaryGeneratorAction( Holstein52Generator*, K37EventGenerator*, K37EventGeneratorNoRecoilOrderEffects*);
	K37MiniaturePGA( Holstein52Generator*);
	~K37MiniaturePGA();

public:
//	void GeneratePrimaries(G4Event* anEvent);  // ok, what gets called after this?  what calls this??
	void GeneratePrimaries();  // ok, what gets called after this?  what calls this??
	G4double getVelocity(G4double kineticEnergy, G4double massOfParticle = 0.510998);
	
//	void SetPolarization(G4double pol);
//	void SetAlignment(G4double ali);
	
//	K37EventGenerator* GetEventGenerator()      { return evGenerator; }  // this seems to break things.
	Holstein52Generator* GetHolsteinGenerator() { return holstein52_generator; }  // does this break things too???
	
	void SetUse_Holstein52(bool tf) { use_holstein52 = tf;   }
	bool GetUse_Holstein52()        { return use_holstein52; }
	
	G4ThreeVector GetEventPositionFromCloud(K37AtomicSetup *);
	void GetMultipoleMomentsFromPops(K37AtomicSetup *);

public:
//	G4double ReturnChargeStateThisEvent() { return recoil_charge_this_event; }
	G4double ReturnIs2p()                 { return is_2p;                    }
private:
//	double recoil_charge_this_event;
	double is_2p;
//	G4double GetChargeStateThisEvent(); // sets recoil_charge_this_event from event generation data.  does this still work?
//	G4double GetIs2p();                 // sets is_2p from event generation data.
	
	
private:
	bool use_holstein52;
//	bool use_gps;
	void NormalizeChargeStateRatio();  // should go to cloud?
	void NormalizeBranchingRatio();    // should go to .... idk, isotope??
	void LoadChargeStates();

	G4bool TwoPercentEvent();

//	G4VPrimaryGenerator * particleGun;
//	G4ParticleDefinition* electron;
//	G4ParticleDefinition* positron;
//	G4ParticleDefinition* gamma;
	
	// Ratios of recoil charge distribution...
//	vector<G4double> charge_state_ratio_;      // Ar0 -> Ar+7
//	vector<G4double> branching_ratio; //Just contains mainbrnach and 2% now. 
//	                       //The other branches enter in at the 0.04% level. 

	G4ThreeVector v;
	G4ThreeVector EventVertex;
	
//	K37EventGenerator                     *evGenerator;
	Holstein52Generator                   *holstein52_generator;  // maybe *this* thing is what owns the atomic setup?  ... yes, I think that's true. Even in TheHolstein.
//	K37EventGeneratorNoRecoilOrderEffects *twoPercent;

	G4bool makeTwoPercent;
	G4bool thisEventIsATwoPercent;

private:
	double P_thisevent;
	double T_thisevent;
	double Mz_thisevent;
	double Mz2_thisevent;
	double Mz3_thisevent;
};

#endif


