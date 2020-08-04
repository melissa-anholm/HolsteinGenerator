// Melissa Anholm - 2019


#include "MiniAggregator.hh"


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
