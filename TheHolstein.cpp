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

// Project:
#include "Holstein52Isotope.hh"    // formerly HolsteinVars
#include "Holstein52Generator.hh"  // formerly HolsteinDecay
#include "K37SublevelPopulations.hh"
#include "K37FermiFunction.hh"

#include "K37MiniaturePGA.hh"

#include "MiniAggregator.hh"
#include "TFile.h"
#include <TTree.h>
#include <TBranch.h> // might not need this...


using std::cout;
using std::endl;

// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // 
string int_to_string(int thisint)
{
	string outstring;
	std::stringstream ss;
	ss.str( std::string() );
	ss.clear();
	ss << thisint;
	outstring = ss.str();
	return outstring;
}

TTree * makefill_tree(int n_events, HolsteinVars * pointervars, K37AtomicSetup * the_atomic_setup, Holstein52Generator *  the_decay, K37MiniaturePGA * the_PGA, bool verbose=false)
{
	TTree * tree = new TTree("ntuple", "ntuple");
	
	Double_t upper_E;
	TBranch *upper_e_b = tree -> Branch("upper_scint_E", &upper_E);
	Double_t lower_E;
	TBranch *lower_e_b = tree -> Branch("lower_scint_E", &lower_E);

	vector<double> * bb1_t_x = 0;
	vector<double> * bb1_t_y = 0;
	vector<double> * bb1_t_r = 0;
	TBranch *bb1_t_x_branch = tree -> Branch("bb1_top_x", &bb1_t_x);
	TBranch *bb1_t_y_branch = tree -> Branch("bb1_top_y", &bb1_t_y);
	TBranch *bb1_t_r_branch = tree -> Branch("bb1_top_r", &bb1_t_r);
	bb1_t_x -> clear();
	bb1_t_y -> clear();
	bb1_t_r -> clear();

	vector<double> * bb1_b_x = 0;
	vector<double> * bb1_b_y = 0;
	vector<double> * bb1_b_r = 0;
	TBranch *bb1_b_x_branch = tree -> Branch("bb1_bottom_x", &bb1_b_x);
	TBranch *bb1_b_y_branch = tree -> Branch("bb1_bottom_y", &bb1_b_y);
	TBranch *bb1_b_r_branch = tree -> Branch("bb1_bottom_r", &bb1_b_r);
	bb1_b_x -> clear();
	bb1_b_y -> clear();
	bb1_b_r -> clear();
	
	double gen_rhit_t;
	double gen_rhit_b;
	TBranch *gen_rhit_t_branch = tree -> Branch("gen_rhit_t", &gen_rhit_t);
	TBranch *gen_rhit_b_branch = tree -> Branch("gen_rhit_b", &gen_rhit_b);

	Int_t sigma_plus = -1;
	TBranch *sigma_plus_branch = tree -> Branch("TTLBit_SigmaPlus", &sigma_plus);  
	//
	vector<double> * scint_time_t = 0;  
	vector<double> * scint_time_b = 0;  
	TBranch * tdc_scint_t_b = tree -> Branch("TDC_SCINT_TOP",    &scint_time_t);  
	TBranch * tdc_scint_b_b = tree -> Branch("TDC_SCINT_BOTTOM", &scint_time_b);  

//	Bool_t det_acceptance = kFALSE;
//	TBranch *det_acceptance_branch = tree -> Branch("det_acceptance", &det_acceptance);
//	Bool_t jtw_acceptance = kFALSE;
//	Bool_t holstein_acceptance = kFALSE;
//	TBranch *jtw_acceptance_branch      = tree -> Branch("jtw_acceptance",      &jtw_acceptance);
//	TBranch *holstein_acceptance_branch = tree -> Branch("holstein_acceptance", &holstein_acceptance);

	Double_t jtw_prob;
	Double_t holstein_prob;
	TBranch *jtw_prob_branch       = tree -> Branch("jtw_probability",      &jtw_prob);
	TBranch *holstein_prob_branch  = tree -> Branch("holstein_probability", &holstein_prob);

	Double_t the_costheta;
	Double_t the_traveltime;
	TBranch *angle_branch = tree -> Branch("gen_costheta",    &the_costheta);   // formerly named "costheta_lab".
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
	
	Double_t x;
	Double_t y;
	Double_t z;
	TBranch *x_b = tree -> Branch("initial_x", &x);  // gen_X_DECAY
	TBranch *y_b = tree -> Branch("initial_y", &y);
	TBranch *z_b = tree -> Branch("initial_z", &z);

	// ...
	Double_t gen_Ebeta_tot;
	Double_t gen_Ebeta_kin;
	Double_t gen_pbeta;
	Double_t gen_beta_beta;
	TBranch *gen_Ebeta_tot_branch = tree -> Branch("gen_Ebeta",     &gen_Ebeta_tot);
	TBranch *gen_Ebeta_kin_branch = tree -> Branch("gen_Tbeta",     &gen_Ebeta_kin);
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
	double naive_hit_t;
	double naive_hit_b;
	TBranch * naive_hit_t_branch = tree -> Branch("naive_hit_t", &naive_hit_t);
	TBranch * naive_hit_b_branch = tree -> Branch("naive_hit_b", &naive_hit_b);
	
	
//	bool verbose = false;	
	bool event_accepted = false;
	G4ThreeVector the_decaypos;
//	int holstein_accept_counter = 0;
	cout << "Let's go!" << endl;
	for(int sigmacount = 0; sigmacount <2; sigmacount++)
	{
		if(sigmacount==0)  // skip sigma plus.  just leave it alone.
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
		for(int i=0; i<n_events/2; i++)
		{
			if( (i % 100000) == 0) 
			{ 
				cout<<"Working on event "<< i/1000 << "k" << endl; 
			}
			
			event_accepted = false;
			bool holstein_acceptance = false;
			bool detector_acceptance = false;
			naive_hit_t = 0.0;
			naive_hit_b = 0.0;
			
			while( !event_accepted ) 
			{ 
				holstein_acceptance = the_PGA->GeneratePrimaries();  // holstein acceptance only.  haven't checked for detector acceptance yet.
				the_decaypos = the_PGA-> Get_DecayPosition();
				detector_acceptance = the_PGA->GetHolsteinGenerator()->check_detector_acceptance(the_decaypos);
				
				
				if(holstein_acceptance && detector_acceptance)
				{
					event_accepted=true;
				}
			}
			if(verbose) 
			{ 
				cout << "* i = " << i << endl;
				the_PGA->GetHolsteinGenerator() -> print_results(); 
			}
			// Fill in the gen data:
			px = the_PGA->GetHolsteinGenerator()->get_beta_Px()/keV;
			py = the_PGA->GetHolsteinGenerator()->get_beta_Py()/keV;
			pz = the_PGA->GetHolsteinGenerator()->get_beta_Pz()/keV;
			
			vx = the_PGA->GetHolsteinGenerator()->initial_velocity.x();
			vy = the_PGA->GetHolsteinGenerator()->initial_velocity.y();
			vz = the_PGA->GetHolsteinGenerator()->initial_velocity.z();
			
			x = the_decaypos.x()/mm;
			y = the_decaypos.y()/mm;
			z = the_decaypos.z()/mm;
			
		//	cout << "Initial Position:  (x,   y,  z) = (" << x << ", " << y << ", " << z << ")" << endl;
		//	cout << "Initial Velocity:  (vx, vy, vz) = " << vx/(m/s) << ", " << vy/(m/s) << ", " << vz/(m/s) << ")" << endl;
			
			the_traveltime = the_PGA->GetHolsteinGenerator()->time_to_travel;
			
			if( vz > 0.0 )
			{
				naive_hit_t = 1.0;
				naive_hit_b = 0.0;
			}
			else
			{
				naive_hit_t = 0.0;
				naive_hit_b = 1.0;
			}
			
			jtw_rho        = the_PGA->GetHolsteinGenerator()->jtw_rho;
			jtw_xi         = the_PGA->GetHolsteinGenerator()->jtw_xi;
			jtw_Abeta      = the_PGA->GetHolsteinGenerator()->jtw_Abeta;
			
			holstein_Abeta = 0;  // this is broken now.  
			
			//
			if( the_PGA->GetHolsteinGenerator()->GetAtomicSetup()->get_sigma() > 0 ) { sigma_plus = 1; }
			else                                                                     { sigma_plus = 0; }
			
			jtw_prob     = the_PGA->GetHolsteinGenerator()->jtw_probability;
			holstein_prob= the_PGA->GetHolsteinGenerator()->holstein_probability;
			the_costheta = the_PGA->GetHolsteinGenerator()->costheta_lab;
			
			gen_Ebeta_tot = the_PGA->GetHolsteinGenerator()->Ebeta_tot_MeV*1000.0;
			gen_Ebeta_kin = the_PGA->GetHolsteinGenerator()->Ebeta_kin_MeV*1000.0;
			gen_pbeta     = the_PGA->GetHolsteinGenerator()->pbeta_MeV*1000.0;
			gen_beta_beta = the_PGA->GetHolsteinGenerator()->vbeta_over_c;
			
			// Fill in stuff that looks like the data ntuples:
			if(the_PGA->GetHolsteinGenerator()->going_up)
			{
				lower_E = 0.0;
				scint_time_t -> push_back(0.0);
				
			//	upper_E = the_PGA->GetHolsteinGenerator()->Ebeta/MeV;
				upper_E = the_PGA->GetHolsteinGenerator()->Ebeta/keV;
				bb1_t_x -> push_back( the_PGA->GetHolsteinGenerator()->hit_position.x()/mm );
				bb1_t_y -> push_back( the_PGA->GetHolsteinGenerator()->hit_position.y()/mm );
				
				gen_rhit_t = sqrt(bb1_t_x->at(0)*bb1_t_x->at(0) + bb1_t_y->at(0)*bb1_t_y->at(0));
				gen_rhit_b = 80.0;
				bb1_t_r -> push_back( gen_rhit_t );
			}
			else if( !the_PGA->GetHolsteinGenerator()->going_up )
			{
				upper_E = 0.0;
				scint_time_b -> push_back(0.0);
			//	lower_E = the_PGA->GetHolsteinGenerator()->Ebeta/MeV;
				lower_E = the_PGA->GetHolsteinGenerator()->Ebeta/keV;
				
				bb1_b_x -> push_back( the_PGA->GetHolsteinGenerator()->hit_position.x()/mm );
				bb1_b_y -> push_back( the_PGA->GetHolsteinGenerator()->hit_position.y()/mm );
				
				gen_rhit_t = 80.0;
				gen_rhit_b = sqrt(bb1_b_x->at(0)*bb1_b_x->at(0) + bb1_b_y->at(0)*bb1_b_y->at(0));
				bb1_b_r -> push_back( gen_rhit_b );
			}
			
			//
			tree -> Fill();
			
			scint_time_t -> clear();
			scint_time_b -> clear();
			bb1_t_x -> clear();
			bb1_t_y -> clear();
			bb1_t_r -> clear();
			
			bb1_b_x -> clear();
			bb1_b_y -> clear();
			bb1_b_r -> clear();
		}
	} // loop over polarizations.
	
	
	
	
	cout << "Done generating events!"   << endl;
	cout << "Total events generated:  " << n_events/2 * 2 << endl;
//	cout << "Holstein accepted events:" << holstein_accept_counter << endl;
	
	
	MiniAggregator * the_miniagg = new MiniAggregator();
	the_miniagg -> PrintMetaDataHeader();
	cout << "MetaHeader is printed." << endl;
//	the_miniagg -> SaveMetaData();
//	the_miniagg -> SaveMetaData(int n_events, HolsteinVars * pointervars, K37AtomicSetup * the_atomic_setup, Holstein52Generator *  the_decay, K37MiniaturePGA * the_PGA)
//	cout << "MetaData is saved." << endl;
	
	
//	the_miniagg -> SaveMetaData(int runno,  string filename, n_events, int nsaved, double cloud_x, double cloud_y, double cloud_z, double cloud_dx, double cloud_dy, double cloud_dz)
	
	return tree;
}

// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // 
int main(int argc, char *argv[]) 
{
	cout << "Called Holstein." << endl;
	
	int nevents=10;
	G4ThreeVector trap_center_i = G4ThreeVector(0.0*mm, 0.0*mm, 1.0*mm);
	G4ThreeVector trap_center_f = G4ThreeVector(0.0*mm, 0.0*mm, 1.0*mm);
	G4ThreeVector trap_size_i   = G4ThreeVector(1.0*mm, 0.5*mm, 0.0*mm);
	G4ThreeVector trap_size_f   = G4ThreeVector(1.0*mm, 0.5*mm, 0.0*mm);
	
	string fname="Output0.root";
	double x_i, y_i, z_i, x_f, y_f, z_f; 
	double sigma_x_i, sigma_y_i, sigma_z_i, sigma_x_f, sigma_y_f, sigma_z_f;
	double t_tmp;
	G4double t_expand = 300.0*microsecond;
	
	string the_paramfile;
	int the_linenumber;
	
//	cout << "argc = " << argc << endl;
//	cout << "argv[1] = " << argv[1] << endl;

	if(argc==2)
	{
		cout << "This is really just a kludge.  Running from the first line of run_params.txt." << endl;
		//
		the_paramfile = "Output/run_params.txt";
		the_linenumber = atoi(argv[1]);
	}
	if(argc==3)
	{
		the_paramfile = "Output/"+string( argv[1] );
		the_linenumber = atoi(argv[2]);
	}
	if( argc==2 || argc==3 )
	{
		std::ifstream ifs(the_paramfile.c_str(), std::ifstream::in);
		if (!ifs.is_open())  { cout << "Error reading " << the_paramfile << endl; }
		
		for(int i=1; i<=the_linenumber && !ifs.eof() && ifs.good(); i++)
		{
			ifs >> fname >> nevents >> t_tmp >> x_i >> y_i >> z_i >> sigma_x_i >> sigma_y_i >> sigma_z_i >> x_f >> y_f >> z_f >> sigma_x_f >> sigma_y_f >> sigma_z_f;
		}
	}
	
	/*
	// old
	if( argc==2 || argc==3 )
	{
		std::ifstream ifs(the_paramfile.c_str(), std::ifstream::in);
		if (!ifs.is_open()) 
		{
			cout << "Error reading " << the_paramfile << endl;
		}
	
		for(int i=1; i<=the_linenumber && !ifs.eof() && ifs.good(); i++)
		{
			ifs >> fname >> nevents >> x >> y >> z >> dx >> dy >> dz;
		}
	}
	*/
	t_expand = t_tmp*microsecond;

	trap_center_i = G4ThreeVector(x_i*mm, y_i*mm, z_i*mm);
	trap_size_i   = G4ThreeVector(sigma_x_i*mm, sigma_y_i*mm, sigma_z_i*mm);
	
	trap_center_f = G4ThreeVector(x_f*mm, y_f*mm, z_f*mm);
	trap_size_f   = G4ThreeVector(sigma_x_f*mm, sigma_y_f*mm, sigma_z_f*mm);
	//
	
	cout << "fname = " << fname << ";\tnevents = " << nevents << "t_expand = " << t_expand/microsecond << " us" << endl;
	cout << "x_i, y_i, z_i    = " << x_i << ", " << y_i << ", " << z_i << endl;
	cout << "x_f, y_f, z_f    = " << x_f << ", " << y_f << ", " << z_f << endl;
	cout << "sx_i, sy_i, sz_i = " << sigma_x_i << ", " << sigma_y_i << ", " << sigma_z_i << endl;
	cout << "sx_f, sy_f, sz_f = " << sigma_x_f << ", " << sigma_y_f << ", " << sigma_z_f << endl;
	
	/*
	if(argc>=2)  // call with args:  ./holstein runno nevents
	{
		runno = atoi(argv[1]);
		nevents = atoi(argv[2]);
	}
	else if(argc==2)
	{
		runno = 0;
		nevents = atoi(argv[1]);
	}
	else
	{
	//	int nhalfevents =        10;
	//	int nhalfevents =       100;
	//	int nhalfevents =      1000;
	//	int nhalfevents =   1000000;
	//	int nhalfevents = 100000000;
	//	int nhalfevents = 1900000000;
	//	int nhalfevents = 300000000*10;  // no, this overflows.
	//	nevents = nhalfevents*2;
		runno = 0;
		nevents = 20;
	}
	*/
	/*
	if( argc==9 )  //   ./holstein runno nevents x y z dx dy 
	{
		
	}
	else
	{
		cout << "Did not call with trap size/position parameters.  Using defaults." << endl;
	}
	*/
	
	
	timeval t1;
	gettimeofday(&t1, NULL);
	G4long randseed = (t1.tv_usec * t1.tv_sec);// + pid;
	CLHEP::HepRandom::setTheSeed(randseed);
	G4cout << "randseed: " << randseed << G4endl;
	
	string outputdir = "Output/";
	string filename = outputdir+fname; //"output_"+int_to_string(runno)+".root";	
	
	HolsteinVars           * pointervars      = new HolsteinVars();	
	K37AtomicSetup         * the_atomic_setup = new K37AtomicSetup();

	the_atomic_setup->SetFreeExpansionTime( t_expand );
	
	the_atomic_setup->SetInitialCloudPosition( trap_center_i );
	the_atomic_setup->SetFinalCloudPosition( trap_center_f );
	the_atomic_setup->SetInitialCloudSize( trap_size_i );
	the_atomic_setup->SetInitialCloudSize( trap_size_f );


	Holstein52Generator    * the_decay        = new Holstein52Generator(pointervars, the_atomic_setup);
	K37MiniaturePGA        * the_PGA          = new K37MiniaturePGA(the_decay);
	
	
//	void SetInitialCloudPosition(G4ThreeVector center) { the_cloud->SetInitialCloudPosition(center); };
//	void SetFinalCloudPosition(G4ThreeVector center)   { the_cloud->SetFinalCloudPosition(center);   };
	
	the_atomic_setup -> SetPolarization(0.9912);
//	the_atomic_setup -> Setup_FromPolarizationAlignment(-0.9912, -0.9736);
	the_atomic_setup -> print_pops();
	the_atomic_setup -> print_moments();
	the_PGA->GetHolsteinGenerator()->set_use_cone(true);
	the_PGA->GetHolsteinGenerator()->set_conecostheta( 0.90 );
	cout << "-- -- --" << endl; // << endl;
	
	
	
	TFile * f = new TFile(filename.c_str(), "RECREATE");
	f -> cd();
	
	cout << "Printing vars!" << endl;
	the_PGA->GetHolsteinGenerator()->print_vars();
	cout << "Vars are printed." << endl;
	
	the_PGA->GetHolsteinGenerator()->use_roc=true;
	
	
	TTree * tree = makefill_tree(nevents, pointervars, the_atomic_setup, the_decay, the_PGA);
	
	tree -> GetCurrentFile() -> Write("",TObject::kOverwrite);  
	tree -> GetCurrentFile() -> Close();
	
	cout << "Events are written and file " << filename << " is closed." << endl;
	return 0;
}










// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // 

