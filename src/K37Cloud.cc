// Authors: Spencer Behling, Benjamin Fenker, Melissa Anholm 2013

#include <fstream>
#include <sstream>

#include "K37Cloud.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"  // it's in K37Cloud.hh now.
#include "G4UnitsTable.hh"

using std::ifstream;

K37Cloud::K37Cloud() : 
	cloud_center_(       G4ThreeVector(0., 0., 0.) ),
	initial_cloud_size_( G4ThreeVector(0., 0., 0.) ),
	temperature_(        G4ThreeVector(0., 0., 0.) ),
	sail_velocity_(      G4ThreeVector(0., 0., 0.) ),
	cycleTime(1.9067*ms),
	//  expansion_before_polarized_(332*microsecond), 
	//  expansion_before_polarized_(400*microsecond), 
	expansion_before_polarized_(300*microsecond), 
	initialize_complete_(false) 
{
	G4cout << "Called K37Cloud::K37Cloud() [with parameters!]" << G4endl;
	//
	G4cout << "Initialized cloud with center: " << cloud_center_ << G4endl;
	G4cout << "                  temperature: " << temperature_ << G4endl;
	G4cout << "                     position: " << initial_cloud_size_ << G4endl;
	G4cout << "                sail velocity: " << sail_velocity_ << G4endl;
}


K37Cloud::~K37Cloud() 
{ 
	G4cout << "Deleting the cloud setup. " << G4endl;
}

void K37Cloud::SetupVelocitySigma(G4ThreeVector temperature) 
{
	velocity_sigma_.setX(CalcSigma(temperature.x()));
	velocity_sigma_.setY(CalcSigma(temperature.y()));
	velocity_sigma_.setZ(CalcSigma(temperature.z()));
}

G4double K37Cloud::CalcSigma(G4double temperature) 
{
	// MB velocity distribution for each dimension is gaussian (normal) centered
	// at zero with standard deviation sigma = sqrt(kT/m)
	// G4 keeps masses in terms of energy so the c*c gives it actual velocity units
	
	// MJA:  Since we're using this to calculate cloud expansion, we want 37K, not 37Ar.
	// G4ParticleTable::GetIon() is obsolete and has been removed.  Use G4IonTable::GetIon() instead.
	G4IonTable * ionTable = G4IonTable::GetIonTable();
	G4double mass = ionTable -> GetIon(19, 37, 0) -> GetPDGMass();
	
	G4double mean_speed = sqrt(temperature * k_Boltzmann * c_squared / mass);
	return mean_speed;
}

void K37Cloud::SetCloudCenter(G4ThreeVector center) 
{
	bool verbose = false;
	cloud_center_ = center;
	if(verbose)
	{
		G4cout << "Updating cloud center to (" << G4BestUnit(center.x(), "Length")
		       << ", " << G4BestUnit(center.y(), "Length") << ", "
		       << G4BestUnit(center.z(), "Length") << ")" << G4endl;
	}
}

void K37Cloud::SetTemperature(G4ThreeVector temp) 
{
	bool verbose=false;
	temperature_ = temp;
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
}

void K37Cloud::SetTemperature(G4double temp) 
{
	SetTemperature(G4ThreeVector(temp, temp, temp));
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
			G4cout << "Updating cloud size to (" << G4BestUnit(size.x(), "Length")
			       << ", " << G4BestUnit(size.y(), "Length") << ", "
			       << G4BestUnit(size.z(), "Length") << ")" << G4endl;
		}
		initial_cloud_size_ = size;
	}
}

void K37Cloud::SetInitialCloudSize(G4double size) 
{
	SetInitialCloudSize(G4ThreeVector(size, size, size));
}

void K37Cloud::Initialize() 
{
	SetupVelocitySigma(temperature_);
	initialize_complete_ = true;
}

void K37Cloud::SetSailVelocity(G4ThreeVector vel) 
{
	bool verbose = false;
	if(verbose)
	{
		G4cout << "Updating sail velocity to (" << vel.x()/(mm/microsecond) << ", "
		       << vel.y()/(mm/microsecond) << ", " << vel.z()/(mm/microsecond) << ")"
		       << " mm / us " << G4endl;
	}
	sail_velocity_ = vel;
}

void K37Cloud::SetFreeExpansionTime(G4double this_time)
{
	bool verbose = false;
	expansion_before_polarized_ = this_time;
	if(verbose)
	{
		G4cout << "Updating expansion time before polarized to:  ";
		G4cout << G4BestUnit(expansion_before_polarized_, "Time") << G4endl;
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
}


