// Authors: Spencer Behling, Benjamin Fenker, Melissa Anholm 2013

#include <iomanip>
#include <cmath>
#include <map>
#include <string>


#include "G4IonTable.hh"
#include "G4SingleParticleSource.hh"
//#include "G4GeneralParticleSource.hh"
#include "G4RandomDirection.hh"  //

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

#include "K37MiniaturePGA.hh"
#include "K37Config.hh"
//#include "K37RunManagerMT.hh"
//#include "K37RunAction.hh"
//#include "K37Data.hh"

using std::string;
using std::ifstream;
using std::map;
//using K37_ABC::K37_Data;

//
//K37MiniaturePGA::K37MiniaturePGA(K37EventGenerator* evGen, K37EventGeneratorNoRecoilOrderEffects* twoP) : 
//K37MiniaturePGA::K37MiniaturePGA( Holstein52Generator* holstein52_gen, K37EventGenerator* evGen, K37EventGeneratorNoRecoilOrderEffects* twoP) : 
K37MiniaturePGA::K37MiniaturePGA( Holstein52Generator* holstein52_gen ) : 
	use_holstein52(true),
//	charge_state_ratio_(8, 1.0/8.0),
//	branching_ratio(2, 1.0/2.0),
	v(),
	EventVertex(),
//	evGenerator(evGen),
	holstein52_generator(holstein52_gen),
//	twoPercent(twoP),
	is_2p(false),
	thisEventIsATwoPercent(false),
	makeTwoPercent(false),
	make_monoenergetic(false),
	monoenergetic_energy(2.5*MeV)
	{
	G4cout << "Created a new MiniPGA." << G4endl;
	
//	new  G4UnitDefinition("eV/c",   "eV/c_light",   "Momentum", eV/c_light);
	new  G4UnitDefinition("kg*m/s", "kg*m/s",       "Momentum", kg*m/s);
	//

//	gunMessenger = new K37PrimaryGeneratorMessenger(this);
//	branching_ratio[0] = 0.9789; //Main branch we care about
//	branching_ratio[1] = 0.0207; //Two percent branch.
//	NormalizeBranchingRatio();
	
	
//	if (use_gps) { particleGun = new G4GeneralParticleSource(); } 
//	else         { particleGun = new G4SingleParticleSource();  }
	
//	particleGun = new G4SingleParticleSource();
//	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
//	G4String particleName;
//	// Ar37_ion = particleTable->GetIon(18, 37, 0);
//	positron = particleTable->FindParticle(particleName="e+");
//	electron = particleTable->FindParticle(particleName="e-");
//	gamma = particleTable->FindParticle(particleName="gamma");
	
	/*
//	charge_state_ratio_[0] = 1.0;
//	charge_state_ratio_[1] = 0.69699;
//	charge_state_ratio_[2] = 0.19483;
//	charge_state_ratio_[3] = 0.08173;
//	charge_state_ratio_[4] = 0.01953;
//	charge_state_ratio_[5] = 0.00581;
//	charge_state_ratio_[6] = 0.00110;
//	charge_state_ratio_[7] = 0.0;
	*/
//	recoil_charge_ =  1.0;  // +1 state by default. 
}

K37MiniaturePGA::~K37MiniaturePGA() 
{
//	delete evGenerator;
	delete holstein52_generator;
//	delete twoPercent; 
}


G4ThreeVector K37MiniaturePGA::GetEventPositionFromCloud(K37AtomicSetup * the_atomic_setup) // replaces void K37AtomicSetup::makeEvent() 
{
	G4double decay_time_ = -10;
	while (decay_time_ < the_atomic_setup->GetFreeExpansionTime() ) 
	{ // pick a new thing for "decay_time_" until you get something longer than "expansion_before_polarized_".
		decay_time_ = (the_atomic_setup->GetOP_CycleTime())*G4UniformRand();  // G4UniformRand() ranges from 0 to 1.
	}
	double fraction_expanded = decay_time_ / the_atomic_setup->GetOP_CycleTime();
	
	G4ThreeVector position_initial = the_atomic_setup->GetInitialCloudPosition();
	G4ThreeVector position_final   = the_atomic_setup->GetFinalCloudPosition();
	G4ThreeVector decay_center     = position_initial + fraction_expanded*(position_final - position_initial);
	
	G4ThreeVector sigma_initial  = the_atomic_setup->GetInitialCloudSize();
	G4ThreeVector sigma_final    = the_atomic_setup->GetFinalCloudSize();
	G4ThreeVector decay_trapsize = sigma_initial + fraction_expanded*(sigma_final - sigma_initial);
	
	G4ThreeVector decay_position = 
		G4ThreeVector(G4RandGauss::shoot(decay_center.x(), decay_trapsize.x()),
					  G4RandGauss::shoot(decay_center.y(), decay_trapsize.y()),
					  G4RandGauss::shoot(decay_center.z(), decay_trapsize.z()) );
	
	EventVertex = decay_position;
	return decay_position;
	//	
}

void K37MiniaturePGA::GetMultipoleMomentsFromPops(K37AtomicSetup * the_atomic_setup)
{
	P_thisevent   = the_atomic_setup->GetPops()->get_P();
	T_thisevent   = the_atomic_setup->GetPops()->get_T();
	Mz_thisevent  = the_atomic_setup->GetPops()->get_Mz();
	Mz2_thisevent = the_atomic_setup->GetPops()->get_Mz2();
	Mz3_thisevent = the_atomic_setup->GetPops()->get_Mz3();
}


//void K37MiniaturePGA::GeneratePrimaries(G4Event* anEvent) 
//void K37MiniaturePGA::GeneratePrimaries() 
bool K37MiniaturePGA::GeneratePrimaries() 
{
	bool the_bool = false;
//	cout << "Does K37MiniaturePGA::GeneratePrimaries() ever even get called?!  ... apparently so!" << endl;
	
//	K37AtomicSetup * the_atomic_setup = K37RunManagerMT::GetMasterRunManager()->GetUserActionInitialization()->GetAtomicSetup();
	K37AtomicSetup * the_atomic_setup = GetHolsteinGenerator()->GetAtomicSetup();
	
	GetMultipoleMomentsFromPops(the_atomic_setup);
//	G4ThreeVector the_position = 
	GetEventPositionFromCloud(the_atomic_setup);
	
	
//	cout << "In GeneratePrimaries():  INitial Position is:  " << EventVertex.x()/mm << ", " << EventVertex.y()/mm << ", " << EventVertex.z()/mm << endl;
	
//	double tmp_P = the_atomic_setup->GetPops()->get_P();
//	double tmp_T = the_atomic_setup->GetPops()->get_T();
//	G4cout << "PGA:GeneratePrimaries(...)[98%+2%]:  the_pops says:  P=" << tmp_P << ";\tT=" << tmp_T << G4endl;
	// Sanity check:  (ie, are any of these things NAN?)
	
	
	/*
	// Debug:
	if(the_atomic_setup==NULL) { G4cout << "omg, there's no the_atomic_setup!"    << G4endl; }
	else                       { G4cout << "phew! we found the the_atomic_setup." << G4endl; }
	K37Cloud * the_actual_cloud = the_atomic_setup->GetCloud();
	if(the_actual_cloud==NULL) { G4cout << "we didn't find the actual cloud.  :(" << G4endl; }
	else                       { G4cout << "we found the actual cloud, too!" << G4endl; }
	K37SublevelPopulations * the_pops = the_atomic_setup->GetPops();
	if(!the_pops)              { G4cout << "we didn't find the pops.  :(" << G4endl; }
	else                       { G4cout << "we found the pops..." << G4endl; }
//	G4cout << "PGA:GeneratePrimaries(...)[98%+2%]:  the_atomic_setup says:  P=" << the_atomic_setup->GetPolarization() << ";\tT=" << the_atomic_setup->GetAlignment() << G4endl;
	*/
	
//	G4cout << "K37MiniaturePGA::GeneratePrimaries(...) is NOT inserting data into the two percent branch right now." << G4endl;
//	// Below:  These fuck shit up.  Moved to K37EventAction, which may or may not also break it.  ...Subsequently moved to K37Run, where it seems fine?
//	(*active_channels_)["TWO_PERCENT_BRANCH"] -> InsertData(is_2p);
//	(*active_channels_)["ION_CHARGE"]         -> InsertData(recoil_charge_this_event);
	
	//
//	EventVertex = GetEventPositionFromCloud(the_atomic_setup);
//	GetEventPositionFromCloud(the_atomic_setup);  // function now sets EventVertex itself.
	
	// Set initial position of all particles where the cloud tells you
//	particleGun->SetParticlePosition(EventVertex);
	
	// Setup initial momenta of B-decay particles
	if(thisEventIsATwoPercent)  // it's not.
	{
	//	twoPercent -> MakeEvent(the_atomic_setup -> GetPolarization(), the_atomic_setup -> GetAlignment() );
	}
	else
	{
		if(use_holstein52)
		{
			if(make_monoenergetic)
			{
				the_bool = holstein52_generator -> shoot_decayevent(monoenergetic_energy);
			}
			else
			{
				the_bool = holstein52_generator -> shoot_decayevent();  // this probably shoots from a different atomic position than we just calculated and will now save to the file...  ...wait, seriously?  No, it's fine -- this function literally doesn't deal with the initial position at all.  ...but it *does* re-calculate it in K37Run, before we save it.  .... wait, still?  Or did we fix it?
			}
		}
	//	else
	//	{
	//		evGenerator -> MakeEvent(the_atomic_setup -> GetPolarization(), the_atomic_setup -> GetAlignment() );
	//	}
	}
	
	return the_bool;
}