// Code by Melissa Anholm
// Feb 2019 - 
// This code is intended to work with the 98% branch of 37K decay.
// this is beta+ decay, and both the parent and daughter have spin 3/2.
// this information is likely baked in somewhere.

// System:
#include <iostream>  // cout, endl
#include <iomanip>   // setw
#include <iterator>

#include <sys/time.h>

#undef NDEBUG
#include<assert.h>


//#include <map>     // for text input function
//#include <cmath>     // pow
//#include "fstream"   // 

// G4:  
#include <Randomize.hh>

//#include<G4ThreeVector.hh> // probably the correct way to include ThreeVector.h.
// #include <CLHEP/Units/SystemOfUnits.h> in HolsteinVars.cpp.  actually G4SystemOfUnits.hh
// #include<G4SystemOfUnits.hh> in HolsteinVars.cpp.  *That* includes <CLHEP/Units/SystemOfUnits.h>.
// #include <G4SystemOfUnits.hh> now in HolsteinDecay.hh

/*
// Root:
#include "TFile.h"
#include <TTree.h>
#include <TBranch.h> // might not need this...
*/

// Project:
#include "Holstein52Isotope.hh"    // formerly HolsteinVars
#include "Holstein52Generator.hh"  // formerly HolsteinDecay
#include "K37SublevelPopulations.hh"
#include "K37FermiFunction.hh"

#include "K37MiniaturePGA.hh"
#include "MiniAggregator.hh"


using std::cout;
using std::endl;

// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // 
// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // 
int main(int argc, char *argv[]) 
{
	cout << "Hello, world!" << endl;
	timeval t1;
	gettimeofday(&t1, NULL);
	G4long randseed = (t1.tv_usec * t1.tv_sec);// + pid;
	CLHEP::HepRandom::setTheSeed(randseed);
	G4cout << "randseed: " << randseed << G4endl;
	
	
	string filename = "output.root";
	
	/*
	K37SublevelPopulations * thepops     = new K37SublevelPopulations(1);
	
//	thepops -> print_pops();
//	thepops -> renormalize();
//	thepops -> print_pops();
//	thepops -> print_moments();
	
//	thepops -> AdjustPolarization(0.9913);  // expect T = -0.9770±0.0017 (+) or T = -0.9761±0.0021 (-).
	
	thepops -> AdjustPolarization(0.9913);
	thepops -> print_pops();
	thepops -> print_moments();
	thepops -> AdjustPolarization(-0.80);
//	thepops -> print_pops();
//	thepops -> print_moments();
	thepops -> AdjustPolarization(0.9913);
	thepops -> print_pops();
	thepops -> print_moments();

//	thepops -> AdjustPolarization(0.95);
//	thepops -> print_pops();
//	thepops -> print_moments();
//	thepops -> AdjustPolarization(0.9);
	
	thepops -> AdjustPolarization(1.0);
	*/
	
	HolsteinVars           * pointervars      = new HolsteinVars();	
	K37AtomicSetup         * the_atomic_setup = new K37AtomicSetup();
//	K37FermiFunction       * the_FF           = new K37FermiFunction();
	
	Holstein52Generator    * the_decay        = new Holstein52Generator(pointervars, the_atomic_setup);
//	Holstein52Generator    * the_decay        = new Holstein52Generator(pointervars, the_atomic_setup, the_FF);
	K37MiniaturePGA        * the_PGA          = new K37MiniaturePGA(the_decay);
	
//	the_PGA->SetMakeMonoenergetic(true);
//	the_PGA->SetMonoenergeticEnergy(2.0*MeV);
	
	the_atomic_setup -> SetPolarization(0.99);
	the_atomic_setup -> print_pops();
	the_atomic_setup -> print_moments();
	
//	the_PGA->GetHolsteinGenerator()->run_fast(true);
	the_PGA->GetHolsteinGenerator()->set_use_cone(true);
	the_PGA->GetHolsteinGenerator()->set_conecostheta( 0.90 );
	
	
	/*
	cout << "I guess we've created the HolsteinDecay." << endl;
	
	K37SublevelPopulations * thepops = the_atomic_setup->GetPops();
	
	cout << "We've got the pops." << endl;
	thepops -> print_pops();
	thepops -> print_moments();
	thepops -> AdjustPolarization(-1.0);
	thepops -> print_pops();
	thepops -> print_moments();
	thepops -> AdjustPolarization(0.99);
	thepops -> print_pops();
	thepops -> print_moments();
	//
	
	the_PGA->GetHolsteinGenerator()->run_fast(true);
	the_PGA->GetHolsteinGenerator()->set_use_cone(true);
//	the_PGA->GetHolsteinGenerator()->set_use_cone(false);
	the_PGA->GetHolsteinGenerator()->set_conecostheta( 0.90 );
	cout << "get_conecostheta() = " << the_PGA->GetHolsteinGenerator()->get_conecostheta() << endl;
	
	
	K37Cloud * the_cloud = the_atomic_setup->GetCloud();
	
	cout << "We've got the cloud." << endl;
	the_cloud->SetInitialCloudPosition( G4ThreeVector(0.1*mm, 1*mm, 1.2*mm) );
	the_cloud->SetInitialCloudSize( G4ThreeVector(2.0*mm, 3.0*mm, 1.0*mm) );
	the_cloud->PrintCloud();
	*/
	
	/*
	cout << "g_A = "  << pointervars->g_A << endl;
	cout << "g_V = "  << pointervars->g_V << endl;
	cout << "M_GT = " << pointervars->M_GT << endl;
	cout << "M_F = "  << pointervars->M_F << endl;

	double rho = pointervars->g_A * pointervars->M_GT / (pointervars->g_V * pointervars->M_F);
	double Abeta = -2.0*rho*( sqrt(3.0/5.0) - rho/5.0 ) / ( 1.0+rho*rho );
	cout << "rho = (g_A*M_GT)/(g_V*M_F) = " << rho << endl;
	cout << "Abeta = -2*rho*(sqrt(3/5) - rho/5) / (1+rho^2) = " << Abeta << endl;
	
	rho *= -1.0;
	cout << "fake rho = " << rho << endl; // 
	Abeta = -2.0*rho*( sqrt(3.0/5.0) - rho/5.0 ) / ( 1.0+rho*rho );
	cout << "Abeta = -2*rho*(sqrt(3/5) - rho/5) / (1+rho^2) = " << Abeta << endl;
	*/
	
	
	cout << "-- -- --" << endl; // << endl;
	/*
	thepops -> Setup_FromPolarizationOnly(0.95);
	thepops -> Setup_FromPolarizationOnly(0.99);
	thepops -> Setup_FromPolarizationOnly(1.0);
//	thepops -> Setup_FromPolarizationOnly(1.5);
	thepops -> Setup_FromPolarizationOnly(-1.0);
	thepops -> Setup_FromPolarizationOnly(-0.95);
	thepops -> Setup_FromPolarizationOnly(-0.9912);
	*/
	/*
	thepops -> Setup_FromPolarizationOnly(0.9912);
	thepops -> print_pops();
	thepops -> print_moments();
	*/
	/*
	thepops -> Setup_FromPolarizationAlignment(1.0, -1.0);
	thepops -> Setup_FromPolarizationAlignment(-0.9912, -0.9736);
	thepops -> Setup_FromPolarizationAlignment(-0.9912, -0.96);
//	thepops -> print_pops();
	*/
	cout << "-- -- --" << endl; // << endl;
	
	
	
	TFile * f = new TFile(filename.c_str(), "RECREATE");
	f -> cd();
	TTree * tree = new TTree("ntuple", "ntuple");
	
	Double_t upper_E;
	TBranch *upper_e_b = tree -> Branch("upper_scint_E", &upper_E);
	Double_t lower_E;
	TBranch *lower_e_b = tree -> Branch("lower_scint_E", &lower_E);
	
	vector<double> * bb1_t_x = 0;
	vector<double> * bb1_t_y = 0;
	TBranch *bb1_t_x_branch = tree -> Branch("bb1_top_x", &bb1_t_x);
	TBranch *bb1_t_y_branch = tree -> Branch("bb1_top_y", &bb1_t_y);
	bb1_t_x -> clear();
	bb1_t_y -> clear();
	
	vector<double> * bb1_b_x = 0;
	vector<double> * bb1_b_y = 0;
	TBranch *bb1_b_x_branch = tree -> Branch("bb1_bottom_x", &bb1_b_x);
	TBranch *bb1_b_y_branch = tree -> Branch("bb1_bottom_y", &bb1_b_y);
	bb1_b_x -> clear();
	bb1_b_y -> clear();
	
	Int_t sigma_plus = -1;
	TBranch *sigma_plus_branch = tree -> Branch("TTLBit_SigmaPlus", &sigma_plus);  
	//
	vector<double> * scint_time_t = 0;  
	vector<double> * scint_time_b = 0;  
	TBranch * tdc_scint_t_b = tree -> Branch("TDC_SCINT_TOP",    &scint_time_t);  
	TBranch * tdc_scint_b_b = tree -> Branch("TDC_SCINT_BOTTOM", &scint_time_b);  
	
//	Bool_t pdf_acceptance = kFALSE;
	Bool_t det_acceptance = kFALSE;
//	TBranch *pdf_acceptance_branch = tree -> Branch("pdf_acceptance", &pdf_acceptance);
	TBranch *det_acceptance_branch = tree -> Branch("det_acceptance", &det_acceptance);
	Bool_t jtw_acceptance = kFALSE;
	Bool_t holstein_acceptance = kFALSE;
	TBranch *jtw_acceptance_branch      = tree -> Branch("jtw_acceptance",      &jtw_acceptance);
	TBranch *holstein_acceptance_branch = tree -> Branch("holstein_acceptance", &holstein_acceptance);
	
//	Double_t the_prob;
//	TBranch *prob_branch           = tree -> Branch("PDF_probability",      &the_prob);
	Double_t jtw_prob;
	Double_t holstein_prob;
	TBranch *jtw_prob_branch       = tree -> Branch("jtw_probability",      &jtw_prob);
	TBranch *holstein_prob_branch  = tree -> Branch("holstein_probability", &holstein_prob);
	
	Double_t the_costheta;
	Double_t the_traveltime;
	TBranch *angle_branch = tree -> Branch("costheta_lab",    &the_costheta);
	TBranch *time_branch  = tree -> Branch("time_to_travel",  &the_traveltime);
	
	Double_t px;
	Double_t py;
	Double_t pz;
	TBranch *px_b = tree -> Branch("initial_px", &px);
	TBranch *py_b = tree -> Branch("initial_py", &py);
	TBranch *pz_b = tree -> Branch("initial_pz", &pz);
	
	Double_t vx;
	Double_t vy;
	Double_t vz;
	TBranch *vx_b = tree -> Branch("initial_vx", &vx);  // what units do these come out in?
	TBranch *vy_b = tree -> Branch("initial_vy", &vy);  
	TBranch *vz_b = tree -> Branch("initial_vz", &vz);  
	
	// ...
	Double_t gen_Ebeta_tot;
	Double_t gen_Ebeta_kin;
	Double_t gen_pbeta;
	Double_t gen_beta_beta;
	TBranch *gen_Ebeta_tot_branch = tree -> Branch("gen_Ebeta_tot", &gen_Ebeta_tot);
	TBranch *gen_Ebeta_kin_branch = tree -> Branch("gen_Ebeta_kin", &gen_Ebeta_kin);
	TBranch *gen_pbeta_branch     = tree -> Branch("gen_pbeta",     &gen_pbeta);
	TBranch *gen_beta_beta_branch = tree -> Branch("gen_betabeta",  &gen_beta_beta);
	
	
	//
	Double_t jtw_Abeta;
	Double_t jtw_xi;
	Double_t jtw_rho;
	Double_t holstein_Abeta;
	TBranch *jtw_Abeta_branch      = tree -> Branch("jtw_Abeta",      &jtw_Abeta);
	TBranch *jtw_xi_branch         = tree -> Branch("jtw_xi",         &jtw_xi);
	TBranch *jtw_rho_branch        = tree -> Branch("jtw_rho",        &jtw_rho);
	TBranch *holstein_Abeta_branch = tree -> Branch("holstein_Abeta", &holstein_Abeta);
	
	
	//
	
	cout << "Printing vars!" << endl;
	the_PGA->GetHolsteinGenerator()->print_vars();
	cout << "Vars are printed." << endl;
	
	bool verbose = false;
	
	the_PGA->GetHolsteinGenerator()->use_roc=true;
	
	bool event_accepted = false;
	int nhalfevents =        10;
//	int nhalfevents =       100;
//	int nhalfevents =      1000;
//	int nhalfevents =   1000000;
//	int nhalfevents = 100000000;
	
//	int mismatch_eventcounter = 0;
//	int jtw_accept_counter = 0;
	int holstein_accept_counter = 0;
	cout << "Let's go!" << endl;
	for(int sigmacount = 0; sigmacount <2; sigmacount++)
	{
		if(sigmacount==0)
		{
			the_PGA->GetHolsteinGenerator()->GetAtomicSetup()->set_sigma_plus();
			cout << "Beginning event generation for sigma+." << endl;
		}
		else if(sigmacount==1)
		{
			the_PGA->GetHolsteinGenerator()->GetAtomicSetup()->set_sigma_minus();
			cout << "Beginning event generation for sigma-." << endl;
		}
		else
		{
			cout << "ono, you broke it  :(" << endl;
			break;
		}
		//
		for(int i=0; i<nhalfevents; i++)
		{
			if( (i % 100000) == 0) 
			{ 
				cout<<"Working on event "<< i << endl; 
			}
			
			event_accepted = false;
		//	int j=0;
			while( !event_accepted ) 
			{ 
			//	cout << "* shot " << j << endl;
			//	event_accepted = the_PGA->GetHolsteinGenerator() -> shoot_decayevent(); 
			//	event_accepted = the_PGA->GetHolsteinGenerator() -> pdf_acceptance;
			//	event_accepted = the_PGA->GetHolsteinGenerator() -> det_acceptance;
				
				if(the_PGA->GetMakeMonoenergetic())
				{
				//	event_accepted = the_PGA->GetHolsteinGenerator() -> shoot_monoenergetic_decayevent( the_PGA->GetMonoenergeticEnergy() );
					event_accepted = the_PGA->GetHolsteinGenerator() -> shoot_decayevent( the_PGA->GetMonoenergeticEnergy() );
				}
				else
				{
					event_accepted = the_PGA->GetHolsteinGenerator() -> shoot_decayevent();
				}
				
				jtw_acceptance = the_PGA->GetHolsteinGenerator()      -> jtw_acceptance;
				holstein_acceptance = the_PGA->GetHolsteinGenerator() -> holstein_acceptance;
				/*
				// below:  this isn't a thing that makes sense anymore.
				if( (jtw_acceptance && !holstein_acceptance) || (!jtw_acceptance && holstein_acceptance) )
				{
					mismatch_eventcounter++;
					if(verbose)
					{
						cout << "** YAY!\ti=" << i << endl;
						cout << "Ebeta = " << the_PGA->GetHolsteinGenerator()->Ebeta_tot_MeV << endl;
						cout << "jtw_acceptance = " << jtw_acceptance << endl;
						cout << "holstein_acceptance = " << holstein_acceptance << endl;
						cout << "pdf_acceptance = " << the_PGA->GetHolsteinGenerator()->pdf_acceptance << endl;
						cout << "event_accepted = " << event_accepted << endl;
					}
				}
				*/
			//	j++;
			}
			if(verbose) 
			{ 
				cout << "* i = " << i << endl;
				the_PGA->GetHolsteinGenerator() -> print_results(); 
			}
			// Fill in the gen data:
			px = the_PGA->GetHolsteinGenerator()->get_beta_Px()/MeV;
			py = the_PGA->GetHolsteinGenerator()->get_beta_Py()/MeV;
			pz = the_PGA->GetHolsteinGenerator()->get_beta_Pz()/MeV;
			
			vx = the_PGA->GetHolsteinGenerator()->initial_velocity.x();
			vy = the_PGA->GetHolsteinGenerator()->initial_velocity.y();
			vz = the_PGA->GetHolsteinGenerator()->initial_velocity.z();
			
		//	pdf_acceptance = the_PGA->GetHolsteinGenerator()->pdf_acceptance;
			det_acceptance = the_PGA->GetHolsteinGenerator()->det_acceptance;  // if runfast is enabled, 
			the_traveltime = the_PGA->GetHolsteinGenerator()->time_to_travel;
			
		//	jtw_acceptance      = the_PGA->GetHolsteinGenerator()->jtw_acceptance;  // not even checked anymore..
			holstein_acceptance = the_PGA->GetHolsteinGenerator()->holstein_acceptance;
			
		//	if(jtw_acceptance)
		//	{
		//		jtw_accept_counter++;
		//	}
			if(holstein_acceptance)
			{
				holstein_accept_counter++;
			}
			
			jtw_rho        = the_PGA->GetHolsteinGenerator()->jtw_rho;
			jtw_xi         = the_PGA->GetHolsteinGenerator()->jtw_xi;
			jtw_Abeta      = the_PGA->GetHolsteinGenerator()->jtw_Abeta;
			holstein_Abeta = the_PGA->GetHolsteinGenerator()->holstein_Abeta;
		//	cout << "JTW:  rho    = " << jtw_rho   << endl;
		//	cout << "      xi     = " << jtw_xi    << endl;
		//	cout << "      A_beta = " << jtw_Abeta << endl;
			
			//
			if( the_PGA->GetHolsteinGenerator()->GetAtomicSetup()->get_sigma() > 0 ) { sigma_plus = 1; }
			else                                                                     { sigma_plus = 0; }
			
		//	the_prob     = the_PGA->GetHolsteinGenerator()->the_probability;
			jtw_prob     = the_PGA->GetHolsteinGenerator()->jtw_probability;
			holstein_prob= the_PGA->GetHolsteinGenerator()->holstein_probability;
			the_costheta = the_PGA->GetHolsteinGenerator()->costheta_lab;
			
			gen_Ebeta_tot = the_PGA->GetHolsteinGenerator()->Ebeta_tot_MeV;
			gen_Ebeta_kin = the_PGA->GetHolsteinGenerator()->Ebeta_kin_MeV;
			gen_pbeta     = the_PGA->GetHolsteinGenerator()->pbeta_MeV;
			gen_beta_beta = the_PGA->GetHolsteinGenerator()->vbeta_over_c;
			
			// Fill in stuff that looks like the data ntuples:
			if(the_PGA->GetHolsteinGenerator()->going_up)
			{
				lower_E = 0.0;
				scint_time_t -> push_back(0.0);
			
				upper_E = the_PGA->GetHolsteinGenerator()->Ebeta/MeV;
				bb1_t_x -> push_back( the_PGA->GetHolsteinGenerator()->hit_position.x()/mm );
				bb1_t_y -> push_back( the_PGA->GetHolsteinGenerator()->hit_position.y()/mm );
			}
			else if( !the_PGA->GetHolsteinGenerator()->going_up )
			{
				upper_E = 0.0;
				scint_time_b -> push_back(0.0);
				lower_E = the_PGA->GetHolsteinGenerator()->Ebeta/MeV;
				
				bb1_b_x -> push_back( the_PGA->GetHolsteinGenerator()->hit_position.x()/mm );
				bb1_b_y -> push_back( the_PGA->GetHolsteinGenerator()->hit_position.y()/mm );
			}
			
			//
			tree -> Fill();
			
			scint_time_t -> clear();
			scint_time_b -> clear();
			bb1_t_x -> clear();
			bb1_t_y -> clear();
			bb1_b_x -> clear();
			bb1_b_y -> clear();
		}
	} // loop over polarizations.
	
	cout << "Done generating events!"  << endl;
	cout << "Total events generated:  " << nhalfevents*2 << endl;
//	cout << "mismatched events:       " << mismatch_eventcounter << endl;
//	cout << "JTW accepted events:     " << jtw_accept_counter << endl;
	cout << "Holstein accepted events:" << holstein_accept_counter << endl;
	
	tree -> GetCurrentFile() -> Write("",TObject::kOverwrite);  
	tree -> GetCurrentFile() -> Close();
	
	cout << "Events are written and file " << filename << " is closed." << endl;
	
	
	return 0;
}










// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // 

