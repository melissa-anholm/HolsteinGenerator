// Melissa Anholm - 2019

#ifndef HolsteinMiniAggregator_h
#define HolsteinMiniAggregator_h 1

// C/C++:
#include <string>
#include <vector>


// Root:
#include "TFile.h"
#include <TTree.h>
#include <TBranch.h> // might not need this...

#include "Holstein52Generator.hh"

using std::string;
using std::vector;

class MiniAggregator
{
public:
	MiniAggregator():metadata_filename( string("Output/MetaData.txt") ) 
	{
		
	};
	
	void PrintMetaDataHeader();
	void SaveMetaData();
//	void SaveMetaData(int n_events, HolsteinVars * pointervars, K37AtomicSetup * the_atomic_setup, Holstein52Generator *  the_decay, K37MiniaturePGA * the_PGA);
	void SaveMetaData(int runno,  string filename, int ngen, int nsaved, double cloud_x, double cloud_y, double cloud_z, double cloud_dx, double cloud_dy, double cloud_dz);
	~MiniAggregator();
	
private:
	string metadata_filename;
};



/*
class MiniAggregator
{
public:
	MiniAggregator( Holstein52Generator * gen );
	~MiniAggregator();
	
	void InitializeBranches();
	void FillBranches();
	void WriteToFile();
	
	string full_filename;
private:
	TFile * f;
	TTree * tree;
	
private:
	Double_t upper_E;
	Double_t lower_E;
	
	vector<double> * bb1_t_x = 0;
	vector<double> * bb1_t_y = 0;
	vector<double> * bb1_b_x = 0;
	vector<double> * bb1_b_y = 0;
	
	Int_t sigma_plus = -1;
	//
	vector<double> * scint_time_t = 0;
	vector<double> * scint_time_b = 0;
	
//	Bool_t pdf_acceptance = kFALSE;
	Bool_t det_acceptance = kFALSE;
	Bool_t jtw_acceptance = kFALSE;
	Bool_t holstein_acceptance = kFALSE;
	
	Double_t the_prob;
	Double_t jtw_prob;
	Double_t holstein_prob;
	
	Double_t the_costheta;
	Double_t the_traveltime;
	
	Double_t px;
	Double_t py;
	Double_t pz;
	
	Double_t vx;
	Double_t vy;
	Double_t vz;
	// ...
	Double_t gen_Ebeta_tot;
	Double_t gen_Ebeta_kin;
	Double_t gen_pbeta;
	Double_t gen_beta_beta;
	//
	Double_t jtw_Abeta;
	Double_t jtw_xi;
	Double_t jtw_rho;
	Double_t holstein_Abeta;
	
	Holstein52Generator * the_generator;
};
*/



#endif

