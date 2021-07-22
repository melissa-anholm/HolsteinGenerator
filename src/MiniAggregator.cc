// Melissa Anholm - 2019


#include "MiniAggregator.hh"

void MiniAggregator::PrintMetaDataHeader()
{
	G4cout << "Called MiniAggregator::PrintMetaDataHeader()." << G4endl;
	
	FILE *io_file;
	io_file = fopen(metadata_filename.c_str(), "a+");
	if(io_file == NULL)
	{
		G4cout << "Couldn't open the MetaData file, for some reason..." << G4endl;
		G4cout << "metadata_filename:  " << metadata_filename << G4endl;
		return;
	}
	
	fprintf(io_file, "%s", "Run/I:");
	fprintf(io_file, "%s", "Filename/C:");
	fprintf(io_file, "%s", "has_been_summed/I:");  // formerly "BadFlag"
	fprintf(io_file, "%s", "is_a_sum/I:");  // formerly "is_summed"
	
	fprintf(io_file, "%s", "SaveEventTypes/C:");
	fprintf(io_file, "%s", "matches_runset/C:");
	
	fprintf(io_file, "%s", "EventsGenerated/I:");
	fprintf(io_file, "%s", "EventsSaved/I:");
	fprintf(io_file, "%s", "MinCosTheta/D:");
	
	// 2% branch?
	
//	fprintf(io_file, "%s", "ChargeState/I:");
	fprintf(io_file, "%s", "Efield/D:");
	
	// Polarization (new notation)
	fprintf(io_file, "%s", "Mz/D:");
	fprintf(io_file, "%s", "Mz2/D:");
	fprintf(io_file, "%s", "Mz3/D:");
	
	fprintf(io_file, "%s", "StepperType/I:");
	fprintf(io_file, "%s", "StepperName/C:");
	fprintf(io_file, "%s", "StepMin_mm/D:"); // is this ... actually the minimum step?
	fprintf(io_file, "%s", "PhysicsListName/C:");
	
	// Cloud parameters:
	fprintf(io_file, "%s", "Trap_x_i_mm/D:");    // initial
	fprintf(io_file, "%s", "Trap_y_i_mm/D:");    // initial
	fprintf(io_file, "%s", "Trap_z_i_mm/D:");    // initial
	fprintf(io_file, "%s", "Trap_sigma_x_i_mm/D:"); // initial
	fprintf(io_file, "%s", "Trap_sigma_y_i_mm/D:"); // initial
	fprintf(io_file, "%s", "Trap_sigma_z_i_mm/D:"); // initial
	
	fprintf(io_file, "%s", "Trap_x_f_mm/D:");    // final
	fprintf(io_file, "%s", "Trap_y_f_mm/D:");    // final
	fprintf(io_file, "%s", "Trap_z_f_mm/D:");    // final
	fprintf(io_file, "%s", "Trap_sigma_x_f_mm/D:"); // final
	fprintf(io_file, "%s", "Trap_sigma_y_f_mm/D:"); // final
	fprintf(io_file, "%s", "Trap_sigma_z_f_mm/D:"); // final
	
	
	fprintf(io_file, "%s", "ExpandBeforePolarized_ms/D:");  // 
	fprintf(io_file, "%s", "OP_CycleTime_ms/D:");  // 
//	fprintf(io_file, "%i\t",   1 );     //  recoil charge goes here.
	
	// Things in old notation that will probably become obsolete soon:
//	fprintf(io_file, "%s", "Rho/D:");
	fprintf(io_file, "%s", "Polarization/D:");
	fprintf(io_file, "%s", "Alignment/D:");
	
	
	fprintf(io_file, "%s", "g_V/D:");
	fprintf(io_file, "%s", "g_A/D:");
	fprintf(io_file, "%s", "g_S/D:");
	fprintf(io_file, "%s", "g_T/D:");
	
	
	fprintf(io_file, "%s", "EventGenerator/C:");
	
	fprintf(io_file, "%s", "MonoEnergy_MeV/D");

	fprintf(io_file, "\n");
	fclose(io_file);
}

/*
void MiniAggregator::SaveMetaData()
{
	G4cout << "Called MiniAggregator::SaveMetaData()." << G4endl;
	
//	if( !setsave_metadata )              { return; } 
	
	if( G4Threading::G4GetThreadId()==0 )  // first worker thread only.  this will be called before it gets called from the master, I think.
	{
		G4RunManager              * runMan           = G4RunManager::GetRunManager();
		K37PrimaryGeneratorAction * kpga             = (K37PrimaryGeneratorAction*)(runMan -> GetUserPrimaryGeneratorAction()); 
		
	//	PGA_eventgenerator_string = kpga->Get_EventGenString();
		PGA_mincostheta = kpga->Get_MinCosTheta();
		
		if( !kpga->GetMakeMonoenergetic() )
		{
			PGA_monoenergy = -10.0*MeV;
		}
		else
		{
			PGA_monoenergy = kpga->GetMonoenergeticEnergy();
		}
		
		G4cout << "PGA_monoenergy = " << PGA_monoenergy << G4endl;
	//	G4cout << "kpga->Get_EventGenString() = " << PGA_eventgenerator_string << G4endl;
		G4cout << "kpga->Get_MinCosTheta() = " << PGA_mincostheta << G4endl;
	}
//	if( G4Threading::IsWorkerThread() )  { return; } // if it's a worker thread, just go away now.
	
	//
	FILE *io_file;
	io_file = fopen(metadata_filename, "a+");  // create the file if it doesn't exist.
	if(io_file==NULL)
	{
		cout << "ERROR:  " << metadata_filename << " could not be opened." << endl;
		return;
	}
	// otherwise, do stuff.
	
	K37RunManagerMT           * MasterRunMan     = K37RunManagerMT::GetMasterRunManager();
//	K37ElectricFieldSetup     * Efield           = MasterRunMan->GetUserDetectorConstruction()->GetField();  // hopefully this doesn't kill it...
//	K37PhysicsList            * physlist         = MasterRunMan->GetPhysicsList();
	
	K37RunAction              * the_runaction    = (K37RunAction*)(G4RunManager::GetRunManager()->GetUserRunAction());
	
	K37AtomicSetup            * the_atomic_setup = (MasterRunMan->GetUserActionInitialization())->GetAtomicSetup();
	K37SublevelPopulations    * the_pops         = the_atomic_setup->GetPops();  // still needed??
	
	HolsteinVars              * the_isotope      = (MasterRunMan->GetUserActionInitialization())->GetHolsteinIsotope();
	
	
	
	if( !the_atomic_setup ) { G4cout << "no the_atomic_setup.  :("  << G4endl; }
	G4cout << "the_atomic_setup->GetPolarization() = " << the_atomic_setup->GetPolarization() << G4endl;
	G4cout << "the_atomic_setup->GetAlignment()    = " << the_atomic_setup->GetAlignment() << G4endl;
	
	fprintf(io_file, "%i\t",   filenumber);
	fprintf(io_file, "%s\t",   GetGlobalMiniName().c_str() );
	fprintf(io_file, "%i\t",   0);  // "has_been_summed/I:"
	fprintf(io_file, "%i\t",   0);  //    "is_a_sum/I:"
	
	fprintf(io_file, "%s\t",   (the_runaction->Get_AcceptanceTypesString()).c_str() ); 
	fprintf(io_file, "%s\t",   (the_atomic_setup->GetMatchedRunsetLetter()).c_str() );     // 
	
	fprintf(io_file, "%i\t",   the_runaction->Get_Nevents_total() );      // number of events generated.
	fprintf(io_file, "%i\t",   the_runaction->Get_Naccepted_total() );    // number of events accepted.
	
	fprintf(io_file, "%f\t",   PGA_mincostheta );  //
	
	// 2% branch?
	fprintf(io_file, "%f\t",   ("EfieldFieldValue" );
	//
	fprintf(io_file, "%f\t",   the_pops->get_Mz() );  // "Mz/D:" 
	fprintf(io_file, "%f\t",   the_pops->get_Mz2() ); // "Mz2/D:" 
	fprintf(io_file, "%f\t",   the_pops->get_Mz3() ); // "Mz3/D:" 
	
	
	fprintf(io_file, "%i\t",   "EfieldStepperTypeIndex" );
	fprintf(io_file, "%s\t",   "EfieldStepperName" );   //
	fprintf(io_file, "%f\t",   "EfieldMinStep" );
	
//	G4cout << "About to look for the EmName from PhysicsList." << G4endl;
	fprintf(io_file, "%s\t",   "physlistEmName");  // this one breaks it.  Or not?

	// true cloud parameters.
	fprintf(io_file, "%f\t",   the_atomic_setup->GetInitialCloudPosition().x()/mm); 
	fprintf(io_file, "%f\t",   the_atomic_setup->GetInitialCloudPosition().y()/mm);  
	fprintf(io_file, "%f\t",   the_atomic_setup->GetInitialCloudPosition().z()/mm);  
	fprintf(io_file, "%f\t",   the_atomic_setup->GetInitialCloudSize().x()/mm);
	fprintf(io_file, "%f\t",   the_atomic_setup->GetInitialCloudSize().y()/mm);
	fprintf(io_file, "%f\t",   the_atomic_setup->GetInitialCloudSize().z()/mm);
	
	fprintf(io_file, "%f\t",   the_atomic_setup->GetFinalCloudPosition().x()/mm); 
	fprintf(io_file, "%f\t",   the_atomic_setup->GetFinalCloudPosition().y()/mm);  
	fprintf(io_file, "%f\t",   the_atomic_setup->GetFinalCloudPosition().z()/mm);  
	fprintf(io_file, "%f\t",   the_atomic_setup->GetFinalCloudSize().x()/mm);
	fprintf(io_file, "%f\t",   the_atomic_setup->GetFinalCloudSize().y()/mm);
	fprintf(io_file, "%f\t",   the_atomic_setup->GetFinalCloudSize().z()/mm);
	
	//
	fprintf(io_file, "%f\t",   the_atomic_setup->GetFreeExpansionTime()/ms);
	fprintf(io_file, "%f\t",   the_atomic_setup->GetOP_CycleTime()/ms);
	
	fprintf(io_file, "%f\t",   the_atomic_setup->GetPolarization() );  // this doesn't seem to work...
	fprintf(io_file, "%f\t",   the_atomic_setup->GetAlignment() );
	
	fprintf(io_file, "%f\t",   the_isotope->get_g_Vector() );
	fprintf(io_file, "%f\t",   the_isotope->get_g_Axial() );
	fprintf(io_file, "%f\t",   the_isotope->get_g_Scalar() );
	fprintf(io_file, "%f\t",   the_isotope->get_g_Tensor() );
	
	fprintf(io_file, "%s\t",   "SimpleHolstein" );  // 
	
	fprintf(io_file, "%f\t",   PGA_monoenergy/MeV );
	
	
	fprintf(io_file, "\n");
	//
	fclose(io_file);
	cout << "Done with SaveMetaData()." << endl;
}
*/

void MiniAggregator::SaveMetaData(int runno,  string filename, int ngen, int nsaved, double cloud_x, double cloud_y, double cloud_z, double cloud_dx, double cloud_dy, double cloud_dz)
{
	string metadata_filename = "MetaInfo.txt";
	FILE *io_file;
	io_file = fopen(metadata_filename.c_str(), "a+");  // create the file if it doesn't exist.
	if(io_file==NULL)
	{
		cout << "ERROR:  " << metadata_filename << " could not be opened." << endl;
		return;
	}
	
	fprintf(io_file, "%s\n",     (filename).c_str() ); 
	fprintf(io_file, "\t%s%i\n", "N_gen=", ngen); 
	
//	fprintf(io_file, "%f\t",   the_isotope->get_g_Tensor() );
//	fprintf(io_file, "%s\t",   (the_runaction->Get_AcceptanceTypesString()).c_str() ); 


	fclose(io_file);
	cout << "Done with GlobalAggregator::SaveMetaData()." << G4endl;
}

/*
MiniAggregator::MiniAggregator(Holstein52Generator * gen)
{
	the_generator = gen;
	//
	full_filename = string("simpleoutput.root");
	f = new TFile(full_filename.c_str(), "RECREATE");
	f -> cd();
	tree = new TTree("ntuple", "ntuple");
}

MiniAggregator::~MiniAggregator()
{
//	tree -> GetCurrentFile() -> Write("",TObject::kOverwrite);  
	tree -> GetCurrentFile() -> Close();
}

void MiniAggregator::InitializeBranches()
{
	TBranch *upper_e_b = tree -> Branch("upper_scint_E", &upper_E);
	TBranch *lower_e_b = tree -> Branch("lower_scint_E", &lower_E);
	
	TBranch *bb1_t_x_branch = tree -> Branch("bb1_top_x", &bb1_t_x);
	TBranch *bb1_t_y_branch = tree -> Branch("bb1_top_y", &bb1_t_y);
	bb1_t_x -> clear();
	bb1_t_y -> clear();
	
	TBranch *bb1_b_x_branch = tree -> Branch("bb1_bottom_x", &bb1_b_x);
	TBranch *bb1_b_y_branch = tree -> Branch("bb1_bottom_y", &bb1_b_y);
	bb1_b_x -> clear();
	bb1_b_y -> clear();
	
	TBranch *sigma_plus_branch = tree -> Branch("TTLBit_SigmaPlus", &sigma_plus);  
	//
	TBranch * tdc_scint_t_b = tree -> Branch("TDC_SCINT_TOP",    &scint_time_t);  
	TBranch * tdc_scint_b_b = tree -> Branch("TDC_SCINT_BOTTOM", &scint_time_b);  
	
	TBranch *det_acceptance_branch = tree -> Branch("det_acceptance", &det_acceptance);
	TBranch *jtw_acceptance_branch      = tree -> Branch("jtw_acceptance",      &jtw_acceptance);
	TBranch *holstein_acceptance_branch = tree -> Branch("holstein_acceptance", &holstein_acceptance);
	
	TBranch *prob_branch           = tree -> Branch("PDF_probability",      &the_prob);
	TBranch *jtw_prob_branch       = tree -> Branch("jtw_probability",      &jtw_prob);
	TBranch *holstein_prob_branch  = tree -> Branch("holstein_probability", &holstein_prob);
	
	TBranch *angle_branch = tree -> Branch("costheta_lab",    &the_costheta);
	TBranch *time_branch  = tree -> Branch("time_to_travel",  &the_traveltime);
	
	TBranch *px_b = tree -> Branch("initial_px", &px);
	TBranch *py_b = tree -> Branch("initial_py", &py);
	TBranch *pz_b = tree -> Branch("initial_pz", &pz);
	
	TBranch *vx_b = tree -> Branch("initial_vx", &vx);  // what units do these come out in?
	TBranch *vy_b = tree -> Branch("initial_vy", &vy);  
	TBranch *vz_b = tree -> Branch("initial_vz", &vz);  
	
	// ...
	TBranch *gen_Ebeta_tot_branch = tree -> Branch("gen_Ebeta_tot", &gen_Ebeta_tot);
	TBranch *gen_Ebeta_kin_branch = tree -> Branch("gen_Ebeta_kin", &gen_Ebeta_kin);
	TBranch *gen_pbeta_branch     = tree -> Branch("gen_pbeta",     &gen_pbeta);
	TBranch *gen_beta_beta_branch = tree -> Branch("gen_betabeta",  &gen_beta_beta);
	//
	TBranch *jtw_Abeta_branch      = tree -> Branch("jtw_Abeta",      &jtw_Abeta);
	TBranch *jtw_xi_branch         = tree -> Branch("jtw_xi",         &jtw_xi);
	TBranch *jtw_rho_branch        = tree -> Branch("jtw_rho",        &jtw_rho);
	TBranch *holstein_Abeta_branch = tree -> Branch("holstein_Abeta", &holstein_Abeta);
}

void MiniAggregator::FillBranches()
{
	
}

void MiniAggregator::WriteToFile()
{
	tree -> GetCurrentFile() -> Write("",TObject::kOverwrite);
}
*/
