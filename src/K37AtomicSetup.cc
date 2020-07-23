// Authors: Spencer Behling, Benjamin Fenker, Melissa Anholm 2013

#include <fstream>
#include <sstream>


#include "globals.hh"
#include "K37Config.hh"
//#include "Randomize.hh"
#include "G4SystemOfUnits.hh"  // it's in K37AtomicSetup.hh now.?  or not.
#include "G4UnitsTable.hh"

//#include "K37Cloud.hh"
//#include "K37SublevelPopulations.hh" // already in the header file.
#include "K37AtomicSetup.hh"


class K37Cloud;
class K37SublevelPopulations;

//using std::ifstream;

K37AtomicSetup::K37AtomicSetup() : 
	MatchedRunsetLetter( G4String("0") )
{
	G4cout << "Called K37AtomicSetup::K37AtomicSetup()" << G4endl;
	
	the_cloud = new K37Cloud();
	the_pops  = new K37SublevelPopulations();
	
	#ifndef SIMPLE_MC
	AtomicMessenger = new K37AtomicMessenger(this);  // this might break it?
	#endif
	
//	cout << "K37AtomicSetup has created a new cloud, and a new set of pops.  In fact, we're done creating the K37AtomicSetup." << endl;
}


K37AtomicSetup::~K37AtomicSetup() 
{ 
	G4cout << "Deleting the atomic setup (including the_cloud and the_pops), but maybe not its messenger... " << G4endl;
	delete the_cloud;
	delete the_pops;
}

// void K37AtomicSetup::SetMatchedRunsetLetter(G4String newRunsetLetter) { MatchedRunsetLetter = newRunsetLetter; };
// G4String K37AtomicSetup::GetMatchedRunsetLetter() { return MatchedRunsetLetter; };
