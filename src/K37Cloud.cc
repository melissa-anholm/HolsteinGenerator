// Authors: Spencer Behling, Benjamin Fenker, Melissa Anholm 2013

#include <fstream>
#include <sstream>

#include "K37Cloud.hh"

#include "globals.hh"
//#include "G4SystemOfUnits.hh"  // it's in K37Cloud.hh now.  ...no it's not.  ...yes it is!
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

#include "G4IonTable.hh"  // Don't use it.  These aren't the masses you're looking for.

using std::ifstream;

K37Cloud::K37Cloud() : 
	initial_position(   G4ThreeVector(0., 0., 0.) ),
	final_position(     G4ThreeVector(0., 0., 0.) ),
	initial_cloud_size( G4ThreeVector(0., 0., 0.) ),
	final_cloud_size(   G4ThreeVector(0., 0., 0.) ),
	cycleTime(1.9067*ms),
	expansion_before_polarized(300*microsecond), 
	gen_h1(false),
	gen_h2(false),
	gen_h3(false),
	gen_h4(false), 
	gen_rmcp(false),
	gen_ring(false),
	gen_hoopedge(false),
	adjacent(0.0),
	max_theta(-5.0)
{
	bool verbose=false;
	if(verbose)
	{
		G4cout << "Called K37Cloud::K37Cloud() [without parameters.]" << G4endl;
	}
	
	if(verbose)
	{
		PrintCloud();
	}
}

K37Cloud::~K37Cloud() 
{ 
	G4cout << "Deleting the cloud. " << G4endl;
}

void K37Cloud::SetInitialCloudPosition(G4ThreeVector center) 
{
	bool verbose = false;
	initial_position = center;
	if(verbose)
	{
		G4cout << "Updating (initial) cloud center to (" << G4BestUnit(center.x(), "Length")
		       << ", " << G4BestUnit(center.y(), "Length") << ", "
		       << G4BestUnit(center.z(), "Length") << ")" << G4endl;
	}
}
void K37Cloud::SetFinalCloudPosition(G4ThreeVector center)
{
	bool verbose = false;
	final_position = center;
	if(verbose)
	{
		G4cout << "Updating (initial) cloud center to (" << G4BestUnit(center.x(), "Length")
		       << ", " << G4BestUnit(center.y(), "Length") << ", "
		       << G4BestUnit(center.z(), "Length") << ")" << G4endl;
	}
}
//
void K37Cloud::SetInitialCloudSize(G4ThreeVector size) 
{
	bool verbose = false;
	if (size.x() < 0 || size.y() < 0 || size.z() < 0) 
	{
		G4cout << "ERROR! Cannot set negative cloud size" << G4endl;
	} 
	else 
	{
		if(verbose)
		{
			G4cout << "Updating initial cloud size to (" << G4BestUnit(size.x(), "Length")
			       << ", " << G4BestUnit(size.y(), "Length") << ", "
			       << G4BestUnit(size.z(), "Length") << ")" << G4endl;
		}
		initial_cloud_size = size;
	}
}
void K37Cloud::SetFinalCloudSize(G4ThreeVector size)
{
	bool verbose = false;
	if (size.x() < 0 || size.y() < 0 || size.z() < 0) 
	{
		G4cout << "ERROR! Cannot set negative cloud size" << G4endl;
	} 
	else 
	{
		if(verbose)
		{
			G4cout << "Updating final cloud size to (" << G4BestUnit(size.x(), "Length")
			       << ", " << G4BestUnit(size.y(), "Length") << ", "
			       << G4BestUnit(size.z(), "Length") << ")" << G4endl;
		}
		final_cloud_size = size;
	}
}

void K37Cloud::PrintCloud() 
{
	if(!gen_h3 && !gen_h4 && !gen_h1 && !gen_h2)
	{
		G4cout << "Cloud Parameters--- --- --- --- --- --- --- --- --- --- " << G4endl;
		G4cout << " Initial Position:  " << initial_position << G4endl;                                 // Cloud position at AC-MOT off.
		G4cout << "   Final Position:  " << final_position << G4endl;                                   // Cloud position at AC-MOT off.
		G4cout << "     Initial Size:  " << initial_cloud_size << G4endl;                               // Cloud size at AC-MOT off.
		G4cout << "       Final Size:  " << final_cloud_size << G4endl;                                 // Cloud size at AC-MOT off.
		G4cout << "   Expansion Time:  " << expansion_before_polarized/millisecond << " ms" << G4endl;  // Time from AC-MOT off to OP on.
		G4cout << "       Cycle Time:  " << cycleTime/millisecond << " ms" << G4endl;                   // Time from AC-MOT off to AC-MOT on.
	//	G4cout << "    Sail velocity:  " << sail_velocity << " (calculated) " << G4endl;                 
	//	G4cout << "      Temperature:  " << temperature << " (caculated?) " << G4endl;                    
		G4cout << "--- --- --- --- --- --- --- --- --- --- --- --- --- --- " << G4endl;
	}
	else
	{
		if(gen_h3)
		{
			G4cout << "\"Cloud\" Parameters- --- --- --- --- --- --- --- --- --- " << G4endl;
			G4cout << "Using Hoop 3 !" << G4endl;
			G4cout << "--- --- --- --- --- --- --- --- --- --- --- --- --- --- " << G4endl;
		}
		else if(gen_h4)
		{
			G4cout << "\"Cloud\" Parameters- --- --- --- --- --- --- --- --- --- " << G4endl;
			G4cout << "Using Hoop 4 !" << G4endl;
			G4cout << "--- --- --- --- --- --- --- --- --- --- --- --- --- --- " << G4endl;
		}
		else
		{
			G4cout << "Cloud Parameters:  it's some sort of hoop or something." << G4endl;
		}
	}
	return;
}

void K37Cloud::SetFreeExpansionTime(G4double this_time)
{
	bool verbose = false;
	expansion_before_polarized = this_time;
	if(verbose)
	{
		G4cout << "Updating expansion time before polarized to:  ";
		G4cout << G4BestUnit(expansion_before_polarized, "Time") << G4endl;
	}
}
void K37Cloud::SetOP_CycleTime(G4double this_time)
{
	bool verbose = false;
	cycleTime = this_time;
	if(verbose)
	{
		G4cout << "Updating total OP cycle time to:  " << G4BestUnit(cycleTime, "Time") << G4endl;
	}
	// set up the sail velocity now?
	// set_up_sail_velocity();
	//
	// set up temperature too?
	// set_up_temperature();
}


