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
//	vertex(NULL),
	EventVertex(),
//	K37Neutral(NULL),
//	K37Minus(NULL),
//	decayTableAr37Minus(NULL),
//	K37MinusDecayMode(NULL),
//	evGenerator(evGen),
	holstein52_generator(holstein52_gen),
//	twoPercent(twoP),
//	make_beta_(true),
//	make_recoil_(true),
//	make_shakeoff_electrons_(true),
//	make_uniform_E_(false),
//	recoil_charge_file_("charge_state_dist.res")//
	is_2p(false),
	thisEventIsATwoPercent(false),
	makeTwoPercent(false)
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
//	delete particleGun;
//	delete gunMessenger;
//	delete K37Neutral;
//	delete K37Minus;
//	delete decayTableAr37Minus;
//	delete K37MinusDecayMode;
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
void K37MiniaturePGA::GeneratePrimaries() 
{
	if(makeTwoPercent /* && recoil_charge_ != -3 */ ) 
	{ // 2% branch is allowed && not photoions.
//		thisEventIsATwoPercent = TwoPercentEvent();
	}
//	recoil_charge_this_event = GetChargeStateThisEvent();
//	is_2p                    = GetIs2p();
	
	
//	K37AtomicSetup * the_atomic_setup = K37RunManagerMT::GetMasterRunManager()->GetUserActionInitialization()->GetAtomicSetup();
	K37AtomicSetup * the_atomic_setup = GetHolsteinGenerator()->GetAtomicSetup();
	
	GetMultipoleMomentsFromPops(the_atomic_setup);
	GetEventPositionFromCloud(the_atomic_setup);
	
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
			holstein52_generator -> shoot_decayevent();
		}
	//	else
	//	{
	//		evGenerator -> MakeEvent(the_atomic_setup -> GetPolarization(), the_atomic_setup -> GetAlignment() );
	//	}
	}
	/*
	if (make_beta_) 
	{
		this->setBetaVertex();  //
		anEvent->AddPrimaryVertex(vertex);
	}
	if (make_recoil_) 
	{
		this->setDaughterVertex(recoil_charge_this_event);
		anEvent->AddPrimaryVertex(vertex);
	}
	if (make_shakeoff_electrons_) 
	{
		SetSOelectronVertices(anEvent, recoil_charge_this_event + 1);
	}
	if (thisEventIsATwoPercent) 
	{
		SetGammaVertex(anEvent);
	}
	*/
}
/*U
void K37MiniaturePGA::setBetaVertex() 
{
	vertex = new G4PrimaryVertex(EventVertex, 0);
	G4PrimaryParticle* particle = 0;
	if(thisEventIsATwoPercent)
	{
		particle = new G4PrimaryParticle(positron, twoPercent->eMomentumX(),
			twoPercent->eMomentumY(), twoPercent->eMomentumZ());
	}
	else
	{
		particle = new G4PrimaryParticle(positron, evGenerator->eMomentumX(),
			evGenerator->eMomentumY(), evGenerator->eMomentumZ());
	}
	vertex->SetPrimary(particle);
}
void K37MiniaturePGA::SetGammaVertex(G4Event *ev)  // for 2% events.
{
	vertex = new G4PrimaryVertex(EventVertex, 0);
	G4double gamma_energy = 2.796 * MeV;
	G4ThreeVector gamma_direction = G4RandomDirection(); //isotropic
	G4PrimaryParticle *twoPercentGamma = new G4PrimaryParticle(gamma);
	twoPercentGamma -> SetMomentumDirection(gamma_direction);
	twoPercentGamma -> SetKineticEnergy(gamma_energy);

	vertex -> SetPrimary(twoPercentGamma);
	ev -> AddPrimaryVertex(vertex);
}
void K37MiniaturePGA::setDaughterVertex(G4double recoil_charge) 
{
	G4bool debug = false;
	G4ThreeVector momentum(0, 0, 0);
	vertex = new G4PrimaryVertex(EventVertex, 0);
	
	G4IonTable * ionTable = G4IonTable::GetIonTable();
	G4ParticleDefinition *Ar37_ion = ionTable -> GetIon(18, 37, 0);  // 37Ar	
	G4PrimaryParticle* particle = 0;
	
	if(thisEventIsATwoPercent)
	{
		particle = new G4PrimaryParticle(Ar37_ion, twoPercent->dMomentumX(),
			twoPercent->dMomentumY(), twoPercent->dMomentumZ());
	}
	else
	{
		particle = new G4PrimaryParticle(Ar37_ion, evGenerator->dMomentumX(),
			evGenerator->dMomentumY(), evGenerator->dMomentumZ());
	}
	// Simulate any charge state >= +1
	particle -> SetCharge(recoil_charge * eplus);
	if (debug) 
	{
		G4cout << "Ion mass: " << Ar37_ion -> GetPDGMass()/c_squared/kg << " kg"
		       << G4endl;
		G4cout << "Ion velocity: " << G4endl << "v0_x = "
		       << (particle -> GetPx()/(particle -> GetMass())*c_light)/(mm/ns)
		       << " mm/ns" << G4endl << "v0_y = "
		       << (particle -> GetPy()/(particle -> GetMass())*c_light)/(mm/ns)
		       << " mm/ns" << G4endl << "v0_z = "
		       << (particle -> GetPz()/(particle -> GetMass())*c_light)/(mm/ns)
		       << " mm/ns" << G4endl;
		G4cout << "Ion momentum: "
		       << particle -> GetTotalMomentum()/(c_light*kg*m/s) << " kg*m/s"
		       << G4endl;
		G4cout << "Ion kinetic energy: " << particle -> GetKineticEnergy()/joule
		       << " J " << G4endl;
	}
	vertex->SetPrimary(particle);
}
void K37MiniaturePGA::SetSOelectronVertices(G4Event *ev, G4int num_so_electrons) 
{
	for (G4int i = 0; i < num_so_electrons; i++) 
	{
		G4double kinetic_energy = -1.0;
		// Ar binding energy is 15.7 eV and SOE have around this energy
		// The 5.0 eV width is a total guess
		while (kinetic_energy < 0.0) 
		{
			kinetic_energy = CLHEP::RandGauss::shoot(100.0*eV, 5.0*eV);
		}
		G4ThreeVector momentum = 
		GetMomentumIsotropic(kinetic_energy, electron_mass_c2);
		
		G4PrimaryParticle *particle = new G4PrimaryParticle(electron,
		                                               momentum.getX(),
		                                               momentum.getY(),
		                                               momentum.getZ());
		vertex = new G4PrimaryVertex(EventVertex, 0);  // 0 means occurs at t = 0
		vertex -> SetPrimary(particle);
		ev -> AddPrimaryVertex(vertex);
	}
}
*/
/*
G4ThreeVector K37MiniaturePGA::GetMomentumIsotropic( G4double kinetic_energy, G4double mass) // used to make SOE vertex.
{
	G4bool debug = false;
	G4double total_energy = kinetic_energy + mass;
	G4double momentum = sqrt((total_energy*total_energy) - (mass * mass));
	G4double mu = 1.0 - 2.0*G4UniformRand();
	G4double theta = acos(mu);
	G4double phi = 2.0*M_PI*G4UniformRand();
	G4double px, py, pz;
	px = momentum * sin(theta) * cos(phi);
	py = momentum * sin(theta) * sin(phi);
	pz = momentum * cos(theta);
	if (debug) 
	{
		G4cout << "\tT = " << G4BestUnit(kinetic_energy, "Energy") << "\tE = "
		       << G4BestUnit(total_energy, "Energy") << "\tP = "
		       << momentum << G4endl;
		G4cout << "\ttheta = " << theta << "\tphi = " << phi << G4endl;
		G4cout << "\tpx = " << px << "\tpy = " << py << "\tpz = " << pz << G4endl;
	}
	return G4ThreeVector(px, py, pz);
}
*/
/*
void K37MiniaturePGA::NormalizeBranchingRatio()
{
	G4double sum = 0.0;
	for (unsigned int i = 0; i < branching_ratio.size(); i++)
	{
		sum += branching_ratio[i];
	}
	for (unsigned int i = 0; i < branching_ratio.size(); i++)
	{
		branching_ratio[i] = branching_ratio[i] / sum;
	}
}
*/

/*
void K37MiniaturePGA::SetRecoilCharge(G4int charge) 
{
	if (charge < 0 && charge != -2 && charge != -3) 
	{
		G4cout << "Negative ions not supported as primary particles.  ";
		G4cout << "No change.  Recoil charge = " << recoil_charge_ << G4endl;
	} 
	else if (charge == -2) 
	{
		recoil_charge_ = charge;
		G4cout << "Set charge state = -2"  << G4endl;
		G4cout << "Charge state will be randomized based on file:  " << recoil_charge_file_ << G4endl;
		LoadChargeStates();         
		NormalizeChargeStateRatio();
	} 
	else if (charge == -3) 
	{
		G4cout << "Set charge state = -3:  ";
		G4cout << "Photoions."  << G4endl;  // How implemented is this??
		recoil_charge_ = charge;
	} 
	else 
	{ // Positive ions action
		recoil_charge_ = charge;
		G4cout << "Set recoil charge = " << recoil_charge_ << G4endl;
	}
}
void K37MiniaturePGA::NormalizeChargeStateRatio() 
{
	G4double sum = 0.0;
	for (unsigned int i = 0; i < charge_state_ratio_.size(); i++) 
	{
		sum += charge_state_ratio_[i];
	}
	for (unsigned int i = 0; i < charge_state_ratio_.size(); i++) 
	{
		charge_state_ratio_[i] = charge_state_ratio_[i] / sum;
	}
}
G4double K37MiniaturePGA::GetChargeStateThisEvent() // does this still work?
{
	G4double this_charge_state = -5.0;
	if (recoil_charge_ == -2) // randomized charge state based on distribution.
	{
		G4double sum = 0.0;
		unsigned an_index = 0;
		G4double guess = CLHEP::RandFlat::shoot(0.0, 1.0);
		while (this_charge_state < 0 && an_index < charge_state_ratio_.size()) 
		{
			sum += charge_state_ratio_[an_index];
			if (guess < sum) 
			{
				this_charge_state = an_index * eplus;
			}
			an_index++;
		}
		G4cout << "this_charge_state = " << this_charge_state << G4endl;
	}
	else if(recoil_charge_ == -3) // photoions.
	{
		this_charge_state = 1.0 * eplus;
	}
	else // recoil_charge_ != -2, so the charge state is given and not randomized.
	{
		this_charge_state = recoil_charge_ * eplus;
	}
	return this_charge_state;
}
*/

/*
G4double K37MiniaturePGA::GetIs2p()
{
	is_2p = -1;
	if(thisEventIsATwoPercent) { is_2p = 1.0; }
	else                       { is_2p = 0.0; }
	return is_2p;
} 
G4bool K37MiniaturePGA::TwoPercentEvent()  
{ // Decide whether or not this event is a 2% branch event.
	G4double guess = CLHEP::RandFlat::shoot(0.0, 1.0);
	if (guess < branching_ratio[1])
	{
		return true;
	}
	else
	{
		return false;
	}
}
*/
/*
void K37MiniaturePGA::SetChargeStatesFile(G4String filename) 
{
	// Set the parameter...
	recoil_charge_file_ = filename;
	G4cout << "Set the recoil charge states file to " << filename << G4endl;
	// Load up the charge states from the file...
	LoadChargeStates();
	NormalizeChargeStateRatio();
	// Automatically set the recoil charge to -2 so the file distribution is used.
	recoil_charge_ = -2.0; 
	return;  
}
void K37MiniaturePGA::LoadChargeStates()
{
	G4String fname_string = G4String(CONFIGURATION_DIRECTORY)+recoil_charge_file_;
	
	G4cout << "* Loading charge states from file:  " << recoil_charge_file_ << G4endl;
	int nstates = 8;
	// Make sure charge_state_ratio_ is the right size...?

	int dummyline = 1;

	G4int cs;
	G4double value;
	G4String nl;

	ifstream file;
	file.open(fname_string.c_str(), std::ifstream::in);
	if (file.is_open()) 
	{
		for(int i=0; i<dummyline; i++)
		{
			getline(file, nl);
		}
		for(int i=0; i<nstates; i++)
		{
			// check if it's not EOF.
			file >> cs >> value;// >> nl;
			getline(file, nl);
			charge_state_ratio_[i] = value;
		}
	}
	else
	{
		G4cout << "* Could not open file " << fname_string << G4endl;
		G4cout << "* Charge state distribution was NOT loaded. " << G4endl;
		G4cout << "* Setting recoil_charge_ to +1." << G4endl;
		recoil_charge_ = 1.0;
	}
	return;
}
G4String K37MiniaturePGA::GetChargeStatesFile()
{
	return recoil_charge_file_;
} 
void K37MiniaturePGA::GetChargeStates()
{
	int nstates = 8;
	G4cout << "Charge state distribution is:  " << G4endl;
	for(int i=0; i<nstates; i++)
	{
		G4cout << "Charge +" << i << ":   " << charge_state_ratio_[i] << G4endl;
	}
	return;
} 
*/
