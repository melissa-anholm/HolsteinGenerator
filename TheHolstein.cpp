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

// Root:
#include "TFile.h"
#include <TTree.h>
#include <TBranch.h> // might not need this...

// Project:
#include "HolsteinVars.hh"
#include "HolsteinDecay.hh"
#include "K37SublevelPopulations.hh"


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
	
	K37SublevelPopulations * thepops     = new K37SublevelPopulations(1);
	
	thepops -> print_pops();
	thepops -> renormalize();
	thepops -> print_pops();
	thepops -> print_moments();
	
	
	HolsteinVars           * pointervars = new HolsteinVars();
	HolsteinDecay          * the_decay   = new HolsteinDecay(pointervars, thepops);
	
	the_decay->run_fast(true);
	the_decay->set_use_cone(true);
//	the_decay->set_use_cone(false);
	the_decay->set_conecostheta( 0.90 );
	cout << "get_conecostheta() = " << the_decay->get_conecostheta() << endl;
	
	
	
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
	
	
	/*
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
	
	Double_t the_prob;
	TBranch *prob_branch           = tree -> Branch("PDF_probability",      &the_prob);
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
	the_decay->print_vars();
	cout << "Vars are printed." << endl;
	
	bool verbose = false;
	
	the_decay->use_roc=true;
	
	bool event_accepted = false;
//	int nhalfevents =       100;
//	int nhalfevents =   1000000;
	int nhalfevents = 100000000;
	
	int mismatch_eventcounter = 0;
	cout << "Let's go!" << endl;
	for(int sigmacount = 0; sigmacount <2; sigmacount++)
	{
		if(sigmacount==0)
		{
			the_decay->the_pops->set_sigma_plus();
			cout << "Beginning event generation for sigma+." << endl;
		}
		else if(sigmacount==1)
		{
			the_decay->the_pops->set_sigma_minus();
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
			//	event_accepted = the_decay -> shoot_decayevent(); 
			//	event_accepted = the_decay -> pdf_acceptance;
			//	event_accepted = the_decay -> det_acceptance;
				
				event_accepted = the_decay -> shoot_decayevent();
				
				jtw_acceptance = the_decay->jtw_acceptance;
				holstein_acceptance = the_decay->holstein_acceptance;
				if( (jtw_acceptance && !holstein_acceptance) || (!jtw_acceptance && holstein_acceptance) )
				{
					mismatch_eventcounter++;
					if(verbose)
					{
						cout << "** YAY!\ti=" << i << endl;
						cout << "Ebeta = " << the_decay->Ebeta_tot_MeV << endl;
						cout << "jtw_acceptance = " << jtw_acceptance << endl;
						cout << "holstein_acceptance = " << holstein_acceptance << endl;
						cout << "pdf_acceptance = " << the_decay->pdf_acceptance << endl;
						cout << "event_accepted = " << event_accepted << endl;
					}
				}
			//	j++;
			}
			if(verbose) 
			{ 
				cout << "* i = " << i << endl;
				the_decay -> print_results(); 
			}
			// Fill in the gen data:
			px = the_decay->initial_momentum.x()/MeV;
			py = the_decay->initial_momentum.y()/MeV;
			pz = the_decay->initial_momentum.z()/MeV;
			
			vx = the_decay->initial_velocity.x();
			vy = the_decay->initial_velocity.y();
			vz = the_decay->initial_velocity.z();
			
		//	pdf_acceptance = the_decay->pdf_acceptance;
			det_acceptance = the_decay->det_acceptance;  // if runfast is enabled, 
			the_traveltime = the_decay->time_to_travel;
			
			jtw_acceptance      = the_decay->jtw_acceptance;
			holstein_acceptance = the_decay->holstein_acceptance;
			
			
			jtw_rho        = the_decay->jtw_rho;
			jtw_xi         = the_decay->jtw_xi;
			jtw_Abeta      = the_decay->jtw_Abeta;
			holstein_Abeta = the_decay->holstein_Abeta;
		//	cout << "JTW:  rho    = " << jtw_rho   << endl;
		//	cout << "      xi     = " << jtw_xi    << endl;
		//	cout << "      A_beta = " << jtw_Abeta << endl;
			
			//
			if( the_decay->the_pops->get_sigma() > 0 ) { sigma_plus = 1; }
			else                                       { sigma_plus = 0; }
			
			the_prob     = the_decay->the_probability;
			jtw_prob     = the_decay->jtw_probability;
			holstein_prob= the_decay->holstein_probability;
			the_costheta = the_decay->costheta_lab;
			
			gen_Ebeta_tot = the_decay->Ebeta_tot_MeV;
			gen_Ebeta_kin = the_decay->Ebeta_kin_MeV;
			gen_pbeta     = the_decay->pbeta_MeV;
			gen_beta_beta = the_decay->vbeta_over_c;
			
			// Fill in stuff that looks like the data ntuples:
			if(the_decay->going_up)
			{
				lower_E = 0.0;
				scint_time_t -> push_back(0.0);
			
				upper_E = the_decay->Ebeta/MeV;
				bb1_t_x -> push_back( the_decay->hit_position.x()/mm );
				bb1_t_y -> push_back( the_decay->hit_position.y()/mm );
			}
			else if( !the_decay->going_up )
			{
				upper_E = 0.0;
				scint_time_b -> push_back(0.0);
				lower_E = the_decay->Ebeta/MeV;
				
				bb1_b_x -> push_back( the_decay->hit_position.x()/mm );
				bb1_b_y -> push_back( the_decay->hit_position.y()/mm );
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
	cout << "mismatched events:       " << mismatch_eventcounter << endl;
	
	tree -> GetCurrentFile() -> Write("",TObject::kOverwrite);  
	tree -> GetCurrentFile() -> Close();
	
	cout << "Events are written and file " << filename << " is closed." << endl;
	*/
	
	return 0;
}










// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // 

