// Authors: Spencer Behling, Benjamin Fenker, Melissa Anholm 2013

#include <fstream>
#include <sstream>


#include "globals.hh"
#include "K37Config.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "K37AtomicSetup.hh"


class K37Cloud;
class K37SublevelPopulations;


K37AtomicSetup::K37AtomicSetup() : 
	MatchedRunsetLetter( G4String("0") )
{
	G4cout << "Called K37AtomicSetup::K37AtomicSetup()" << G4endl;
	
	the_cloud = new K37Cloud();
	the_pops  = new K37SublevelPopulations();
	
	#ifndef SIMPLE_MC
	AtomicMessenger = new K37AtomicMessenger(this);  // this might break it?
	#endif
}


K37AtomicSetup::~K37AtomicSetup() 
{ 
	G4cout << "Deleting the atomic setup (including the_cloud and the_pops), but maybe not its messenger... " << G4endl;
	delete the_cloud;
	delete the_pops;
}