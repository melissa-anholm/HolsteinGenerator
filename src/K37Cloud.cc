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
//	temperature(        G4ThreeVector(0., 0., 0.) ),
//	sail_velocity(      G4ThreeVector(0., 0., 0.) ),
	cycleTime(1.9067*ms),
	expansion_before_polarized(300*microsecond)
{
	bool verbose=false;
	if(verbose)
	{
		G4cout << "Called K37Cloud::K37Cloud() [without parameters.]" << G4endl;
	}
	
//	set_up_sail_velocity();
//	set_up_temperature();
//	SetupVelocitySigma(temperature);
	
	if(verbose)
	{
		PrintCloud();
	}
}

K37Cloud::~K37Cloud() 
{ 
	G4cout << "Deleting the cloud. " << G4endl;
}

/*
void K37Cloud::SetupVelocitySigma(G4ThreeVector temperature) 
{
//	G4cout << "Setting velocity sigma." << G4endl;
	
	velocity_sigma.setX(CalcSigma(temperature.x()));
	velocity_sigma.setY(CalcSigma(temperature.y()));
	velocity_sigma.setZ(CalcSigma(temperature.z()));
	
//	G4cout << "Done setting velocity sigma." << G4endl;
}
G4double K37Cloud::CalcSigma(G4double temperature) 
{
	// MB velocity distribution for each dimension is gaussian (normal) centered
	// at zero with standard deviation sigma = sqrt(kT/m)
	// G4 keeps masses in terms of energy so the c*c gives it actual velocity units
	
	// MJA:  Since we're using this to calculate cloud expansion, we want 37K, not 37Ar.
	// G4ParticleTable::GetIon() is obsolete and has been removed.  Use G4IonTable::GetIon() instead.
//	G4IonTable * ionTable = G4IonTable::GetIonTable();
//	G4double mass = ionTable -> GetIon(19, 37, 0) -> GetPDGMass();  // Z, A, J ???	
	
	// Horrible kludge!!
	// Below is the G4-accepted way to get the mass.  (note:  the third argument is spin, not charge!)
	// GetPDGMass() returns something in units of energy -- ie:  MeV, *not* MeV/c^2.
	// I don't know what units GetAtomicMass() comes out in, but it's clearly neither energy nor mass. 
	// In any case, I want the cloud temperature to be calculated initially before calling run/initialize,
	// however the G4IonTable can't be constructed yet at this stage of the program.  We really just need
	// to get the mass for a quick speed calculation, so we'll just kludge in a value instead.  
	//	// G4IonTable * ionTable = G4IonTable::GetIonTable();
	//	// G4double mass = ionTable -> GetIon(19, 37, 0) -> GetPDGMass();     // MeV [not MeV/c^2]
	//	// G4double mass = ionTable -> GetIon(19, 37, 0) -> GetAtomicMass();  // unknown units.
	//	// G4cout << "Mass of 37K is:  " << std::setprecision(32) << mass / (amu) << " amu" << G4endl;  // no it's not.
	G4double mass = 36.97337589 * amu;
	G4double mean_speed = sqrt(temperature * k_Boltzmann / mass); 
	
	bool verbose = false;
	if(verbose)
	{
		G4cout << "Mass of 37K is:  " << std::setprecision(32) << mass / MeV*c_squared << " MeV" << G4endl;
		G4cout << "v_rms = " << mean_speed / (meter/second) << " m/s" << G4endl;
		G4cout << "v_rms = " << G4BestUnit(mean_speed, "Velocity") << G4endl;
	}
	
//	G4IonTable * ionTable = G4IonTable::GetIonTable();  // mass comes out in some kind of energy units, with c=1.  doesn't really match with the other Geant4 units ...
//	mass = ionTable -> GetIon(19, 37, 0) -> GetPDGMass();
//	G4cout << "Mass of 37K is:  " << std::setprecision(32) << mass / MeV << " MeV" << G4endl;
//	mean_speed = sqrt(temperature * k_Boltzmann * c_squared / mass); 
//	G4cout << "v_rms = " << mean_speed / (meter/second) << " m/s" << G4endl;
//	G4cout << "v_rms = " << G4BestUnit(mean_speed, "Velocity") << G4endl;

	return mean_speed;
}
*/

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
void K37Cloud::SetInitialCloudSize(G4double size) 
{
	SetInitialCloudSize(G4ThreeVector(size, size, size));
}
//
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
void K37Cloud::SetFinalCloudSize(G4double size) 
{
	SetFinalCloudSize(G4ThreeVector(size, size, size));
}



void K37Cloud::PrintCloud() 
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
	return;
}

/*
//void K37Cloud::SetSailVelocity(G4ThreeVector vel) 
//{
//	bool verbose = false;
//	if(verbose)
//	{
//		G4cout << "Updating sail velocity to (" << vel.x()/(mm/microsecond) << ", "
//		       << vel.y()/(mm/microsecond) << ", " << vel.z()/(mm/microsecond) << ")"
//		       << " mm / us " << G4endl;
//	}
//	sail_velocity = vel;
//}
*/
void K37Cloud::set_up_sail_velocity()  // initial_position, final_position, and cycleTime must already be set up.
{
	G4ThreeVector sail_velocity = (final_position - initial_position)/cycleTime;
}

/*
void K37Cloud::SetTemperature(G4ThreeVector temp) 
{
	bool verbose=false;
	temperature = temp;
	if (temp.x() < 0 || temp.y() < 0 || temp.z() < 0) 
	{
		G4cout << "ERROR! Cannot set negative temperature." << G4endl;
	} 
	else
	{
		if(verbose)
		{
			G4cout << "Updating cloud temp to (" << G4BestUnit(temp.x(), "Temperature")
			       << ", " << G4BestUnit(temp.y(), "Temperature") << ", "
			       << G4BestUnit(temp.z(), "Temperature") << ")" << G4endl;
		}
		SetupVelocitySigma(temp);
	}
	// go backwards to set up the final size now?
}
void K37Cloud::SetTemperature(G4double temp) 
{
	SetTemperature(G4ThreeVector(temp, temp, temp));
	G4cout << "Set cloud temperature to:  " << temp / kelvin << " Kelvin." << G4endl;
}
*/

/*
void K37Cloud::set_up_temperature()    // must have initial_size, final_size, and cycleTime already set up.
{
	// initial_cloud_size
	// 
}
*/

void K37Cloud::SetFinalCloudPosition(G4ThreeVector center)
{
	bool verbose = false;
	initial_position = center;
	if(verbose)
	{
		G4cout << "Updating (initial) cloud center to (" << G4BestUnit(center.x(), "Length")
		       << ", " << G4BestUnit(center.y(), "Length") << ", "
		       << G4BestUnit(center.z(), "Length") << ")" << G4endl;
	}
	// set up the sail velocity now?
	// set_up_sail_velocity();
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


