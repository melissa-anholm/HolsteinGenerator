#include <iostream> // cout, endl
#include <vector>   // 
using std::cout;
using std::endl;
using std::string;
using std::vector;

#include <TFile.h>
#include <TTree.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "TPaveText.h"


#include "HistExtras.cpp"
#include "ColorExtras.cpp"
#include "AsymmetryCanvasLibs.cpp"
#include "MakeSomeTF1s.cpp"

#undef NDEBUG
#include<assert.h>

#define __SHORT_FORM_OF_FILE__ \
(strrchr(__FILE__,'/') ? strrchr(__FILE__,'/')+1 : __FILE__ )


// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
void printresults_A(TFitResultPtr thisfit)
{
	cout << thisfit -> GetName() << ":   " << endl;
	cout << "A = " << thisfit->Parameter( thisfit->Index("Abeta") );
	cout << " +/- " << thisfit->ParError( thisfit->Index("Abeta") );
	cout << "\tchi^2/NDF = " << thisfit->Chi2() / double(thisfit->Ndf() - thisfit->NFreeParameters());
	cout << "\t[chi^2 = " << thisfit->Chi2() << ", NDF = " << thisfit->Ndf() - thisfit->NFreeParameters() << "]" << endl;
}
void printresults_A_b(TFitResultPtr thisfit)
{
	cout << thisfit -> GetName() << ":   " << endl;
	cout << "A = " << thisfit->Parameter( thisfit->Index("Abeta") );
	cout << " +/- " << thisfit->ParError( thisfit->Index("Abeta") );
	cout << "\tb = " << thisfit->Parameter( thisfit->Index("bFierz") );
	cout << " +/- " << thisfit->ParError( thisfit->Index("bFierz") );
	cout << "\tchi^2/NDF = " << thisfit->Chi2() / double(thisfit->Ndf() - thisfit->NFreeParameters());
	cout << "\t[chi^2 = " << thisfit->Chi2() << ", NDF = " << thisfit->Ndf() - thisfit->NFreeParameters() << "]" << endl;
}
void printresults_betaphasespace(TFitResultPtr thisfit)
{
	cout << thisfit -> GetName() << ":   " << endl;
	cout << "scale = " << thisfit->Parameter( thisfit->Index("scale") );
	cout << " +/- " << thisfit->ParError( thisfit->Index("scale") );
	cout << "\tE0 = " << thisfit->Parameter( thisfit->Index("E0") );
	cout << " +/- " << thisfit->ParError( thisfit->Index("E0") );
	cout << "\tchi^2/NDF = " << thisfit->Chi2() / double(thisfit->Ndf() - thisfit->NFreeParameters());
	cout << "\t[chi^2 = " << thisfit->Chi2() << ", NDF = " << thisfit->Ndf() - thisfit->NFreeParameters() << "]" << endl;
}
void printresults_linefunc(TFitResultPtr thisfit)
{
	cout << thisfit -> GetName() << ":   " << endl;
	cout << "slope = " << thisfit->Parameter( thisfit->Index("slope") );
	cout << " +/- " << thisfit->ParError( thisfit->Index("slope") );
	cout << "\tx_intercept = " << thisfit->Parameter( thisfit->Index("x_intercept") );
	cout << " +/- " << thisfit->ParError( thisfit->Index("x_intercept") );
	cout << "\tchi^2/NDF = " << thisfit->Chi2() / double(thisfit->Ndf() - thisfit->NFreeParameters());
	cout << "\t[chi^2 = " << thisfit->Chi2() << ", NDF = " << thisfit->Ndf() - thisfit->NFreeParameters() << "]" << endl;
}
double pbeta(double E)  // E in keV, p in keV/c
{
	double m_e = 0.5109989461*1000.0;
	double pbetac = sqrt( pow(E,2) - pow(m_e,2) );  // MeV for everything...
	return pbetac;
}
bool HistsHaveSameBinning2(TH1D *a, TH1D *b, bool verbose=false) 
{
	bool same = true;
	if (!a || !b) 
	{
		cout << "ERROR:  Histogram doesn't exist" << endl;
		cout << "a=" << a << ", b=" << b << endl;
		same = false;
	//	return same;
	}
	else if ( a -> GetNbinsX() != b -> GetNbinsX() ) 
	{
		cout << "ERROR:  Histograms have different numbers of bins." << endl;
		same = false;
	//	return same;
	}
	double eps = 1.E-3;
	if (same) 
	{
		for (int i = 1; i <= a -> GetNbinsX(); i++) 
		{
			if (fabs(a->GetBinCenter(i) - b->GetBinCenter(i)) > eps)
			{
				same = false;
			}
		}
	}
	//
	if(same && verbose)
	{
		cout << "Histograms " << a->GetName() << " and ";
		cout << b->GetName() << " have the same binning." << endl;
	}
	else if(!same)
	{
		cout << "ERROR:  bin centres are different." << endl;
	}
	return same;
}
TH1D * tf1_to_hist_like(TF1 * this_tf1, TH1D * this_hist, int this_color=int(kBlack) )
{
	TH1D * new_hist = (TH1D*)this_hist -> Clone( this_tf1->GetName() );
	new_hist -> GetListOfFunctions() -> Clear();
	
	new_hist -> Sumw2(kFALSE);
	new_hist -> SetMarkerColor(this_color);
	new_hist -> SetLineColor(this_color);
	int n_bins = new_hist->GetNbinsX();
	// well shit, it probably needs some parameters or something.
	double xmax, xmin;
	xmax = this_tf1->GetXmax();
	xmin = this_tf1->GetXmin();
	
	for (int i=1; i<n_bins; i++)  // Bins i=0, i=n_bins are the underflow and overflow?
	{
	//	this_tf1->Eval(this_hist->GetBinCenter(i));
	//	cout << "i = " << i << endl;
		if( new_hist->GetBinCenter(i)>=xmin &&  new_hist->GetBinCenter(i)<=xmax )
		{
			new_hist -> SetBinContent(i, this_tf1->Eval(this_hist->GetBinCenter(i)) );
		}
		else
		{
			new_hist -> SetBinContent(i, 0);
		}
	}
//	cout << "Returning the hist." << endl;
	return new_hist;
}

// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
int main(int argc, char *argv[]) 
{
	cout << "Hello, world!" << endl;
	TApplication* rootapp = 0;
	char ** mychar = NULL;
	rootapp = new TApplication("blarg",0, mychar);
	
	string filename = "output.root";
//	string filename = "output_partial.root";
//	string filename = "output_10M.root";
//	string filename = "output.root";
	TFile * f = new TFile(filename.c_str(), "READ");
	cout << "Reading from file:  " << filename << endl;
	f -> cd();
	TTree * tree = (TTree*)f->Get("ntuple");
	
	double m_e = 0.5109989461;
	double E0 = 5.636004;  // calculated E0? 
	double speed_of_light = 299792458.0;
	
	
	vector<double> * bb1_t_x = 0;
	vector<double> * bb1_t_y = 0;
	vector<double> * bb1_b_x = 0;
	vector<double> * bb1_b_y = 0;
	tree -> SetBranchAddress("bb1_top_x",   &bb1_t_x);
	tree -> SetBranchAddress("bb1_top_y",   &bb1_t_y);
	tree -> SetBranchAddress("bb1_bottom_x",&bb1_b_x);
	tree -> SetBranchAddress("bb1_bottom_y",&bb1_b_y);
	Double_t upper_E;
	Double_t lower_E;
	tree -> SetBranchAddress("upper_scint_E", &upper_E);
	tree -> SetBranchAddress("lower_scint_E", &lower_E);
	Int_t sigma_plus = -1;
	tree -> SetBranchAddress("TTLBit_SigmaPlus", &sigma_plus);  
	vector<double> * scint_time_t = 0;  
	vector<double> * scint_time_b = 0;  
	tree -> SetBranchAddress("TDC_SCINT_TOP",    &scint_time_t);  
	tree -> SetBranchAddress("TDC_SCINT_BOTTOM", &scint_time_b);  

//	Double_t probability;
//	tree -> SetBranchAddress("PDF_probability", &probability);
	Double_t holstein_probability;
	Double_t jtw_probability;
	tree -> SetBranchAddress("holstein_probability", &holstein_probability);
	tree -> SetBranchAddress("jtw_probability", &jtw_probability);
	Double_t the_costheta;
	tree -> SetBranchAddress("costheta_lab", &the_costheta);
	Double_t time_to_travel;
	tree -> SetBranchAddress("time_to_travel", &time_to_travel);
	
	Bool_t holstein_acceptance = kFALSE;
	Bool_t jtw_acceptance = kFALSE;
	Bool_t det_acceptance = kFALSE;
	tree -> SetBranchAddress("holstein_acceptance", &holstein_acceptance);
	tree -> SetBranchAddress("jtw_acceptance", &jtw_acceptance);
	tree -> SetBranchAddress("det_acceptance", &det_acceptance);
	
	Double_t px;
	Double_t py;
	Double_t pz;
	tree -> SetBranchAddress("initial_px", &px);
	tree -> SetBranchAddress("initial_py", &py);
	tree -> SetBranchAddress("initial_pz", &pz);
	Double_t vx;
	Double_t vy;
	Double_t vz;
	tree -> SetBranchAddress("initial_vx", &vx);
	tree -> SetBranchAddress("initial_vy", &vy);
	tree -> SetBranchAddress("initial_vz", &vz);
	
	// 

//	TH1D* hist_E_t = CreateHist(string("Top Beta Kinetic Energy"),    string("Ben_Ebeta"), int(mDarkRed), 1);
//	TH1D* hist_E_b = CreateHist(string("Bottom Beta Kinetic Energy"), string("Ben_Ebeta"), int(kGreen), 1);
	/*
	TH1D* hist_E_tp = CreateHist(string("Beta Kinetic Energy (Top+)"),    string("Ben_Ebeta"), int(kGreen),   1);
	TH1D* hist_E_bp = CreateHist(string("Beta Kinetic Energy (Bottom+)"), string("Ben_Ebeta"), int(kBlue), 1);
	
	TH1D* hist_E_tm = CreateHist(string("Beta Kinetic Energy (Top-)"),    string("Ben_Ebeta"), int(kOrange),  1);
	TH1D* hist_E_bm = CreateHist(string("Beta Kinetic Energy (Bottom-)"), string("Ben_Ebeta"), int(kMagenta), 1);
	*/
	
	TH2D * bb1_top_position    = CreateHist2d(string("Hit Position Top"),    string("bb1_x"), string("bb1_y"));
	TH2D * bb1_bottom_position = CreateHist2d(string("Hit Position Bottom"), string("bb1_x"), string("bb1_y"));
	
//	TH2D * prob_v_E_more = CreateHist2d( string("Scaled Probability? (more)"), string("Ben_Ebeta_4096"), string("scaled_probability"), 4, 2);
//	TH2D * prob_v_E_less = CreateHist2d( string("Scaled Probability? (less)"), string("Ben_Ebeta_4096"), string("scaled_probability"), 4, 2);

	/*
	// hist_costheta is still kinda cool.  maybe useful later.
	TH1D* hist_costheta = CreateHist(string("Cos(theta_lab)"),    string("costheta"), int(kBlack),   1);
	hist_costheta -> GetXaxis() -> SetRangeUser(-1.0, -0.96);
	*/
	
//	TH2D * costheta_v_E = CreateHist2d( string("cos theta v. E, and probability."), string("Ben_Ebeta_4096"), string("costheta"), 8, 2);
//	TH2D * prob_v_costheta = CreateHist2d( string("Probability(?) of Costheta"), string("costheta"), string("scaled_probability"), 2, 2);
	
	// travel_time_v_costheta is stupid, but I need to leave it because I still have to figure out what units I'm even using.
	TH2D * travel_time_v_costheta = CreateHist2d( string("traveltime vs Costheta"), string("costheta"), string("holstein_betatof"), 1, 1);
	
	int color_holstein = int(kRed);
	int color_jtw      = int(kBlue);
	
	int kurie_rebin = 16*2;
	TH1D * hist_E_both_holstein = CreateHist(string("Beta Total Energy (Holstein)"), string("Ben_Ebeta_4096"), color_holstein, kurie_rebin);
	hist_E_both_holstein -> GetXaxis() -> SetTitle("Beta Energy [keV]");
	hist_E_both_holstein -> GetXaxis() -> SetRangeUser(0.0, 6000.0);
	hist_E_both_holstein -> Sumw2();
	
	TH1D * hist_E_both_jtw = CreateHist(string("Beta Total Energy (JTW)"), string("Ben_Ebeta_4096"), color_jtw, kurie_rebin);
	hist_E_both_jtw -> GetXaxis() -> SetTitle("Beta Energy [keV]");
	hist_E_both_jtw -> GetXaxis() -> SetRangeUser(0.0, 6000.0);
	hist_E_both_jtw -> Sumw2();
	
	//
	TH1D * kurie_plot_holstein = CreateHist(string("Kurie Plot (Holstein)"), string("Ben_Ebeta_4096"), color_holstein, kurie_rebin);
	kurie_plot_holstein -> GetXaxis() -> SetTitle("Beta Energy [keV]");
	kurie_plot_holstein -> GetYaxis() -> SetTitle("sqrt( P(E) / (p*E) ) ");
	kurie_plot_holstein -> GetXaxis() -> SetRangeUser(0.0, 6000.0);
	kurie_plot_holstein -> Sumw2();
	
	TH1D * kurie_plot_jtw = CreateHist(string("Kurie Plot (JTW)"), string("Ben_Ebeta_4096"), color_jtw, kurie_rebin);
	kurie_plot_jtw -> GetXaxis() -> SetTitle("Beta Energy [keV]");
	kurie_plot_jtw -> GetYaxis() -> SetTitle("sqrt( P(E) / (p*E) ) ");
	kurie_plot_jtw -> GetXaxis() -> SetRangeUser(0.0, 6000.0);
	kurie_plot_jtw -> Sumw2();
	
	// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
	// Look through events.
	Long64_t nentries = tree->GetEntries();
	double rhit = -1.0;
	bool accept_holstein=false;
	bool accept_jtw=false;
	
	int nmaxevents = 200000;
	cout << "nentries = " << nentries << ";\tnmaxevents = " << nmaxevents << endl;
//	for(int i=0; (i<nentries) && (i<nmaxevents); i++)
	for(int i=0; (i<nentries); i++)
	{
		if( (i % 100000) == 0) { cout<<"Reached entry "<< i << endl; }
		tree -> GetEntry(i);
		
		if(scint_time_t->size()>0 && scint_time_b->size() > 0 )
			{ cout << "That's bad.  i = "     << i << endl;  assert(0);  return 0; }
		else if(scint_time_t->size()==0 && scint_time_b->size()==0 )
			{ cout << "That's bad too.  i = " << i << endl;  assert(0);  return 0; }
		
		//
		if(scint_time_t->size()>0) { rhit = sqrt( pow((*bb1_t_x)[0], 2) + pow((*bb1_t_y)[0], 2)); }
		if(scint_time_b->size()>0) { rhit = sqrt( pow((*bb1_b_x)[0], 2) + pow((*bb1_b_y)[0], 2)); }
		
	//	accept=true;  // all events.
	//	accept = pdf_acceptance;
	//	accept = det_acceptance;
		accept_holstein = ( holstein_acceptance&&det_acceptance &&(i%2==0) );
		accept_jtw      = ( jtw_acceptance&&det_acceptance &&(i%2==1) );
	//	accept_holstein = ( holstein_acceptance&&det_acceptance);
	//	accept_jtw      = ( jtw_acceptance&&det_acceptance);
		if(rhit>15.5) { accept_holstein=false;  accept_jtw=false; }
		if(accept_holstein)
		{
			hist_E_both_holstein -> Fill( (upper_E+lower_E)*1000.0 );
		}
		if(accept_jtw)
		{
			hist_E_both_jtw -> Fill( (upper_E+lower_E)*1000.0 );
		}
			/*
		//	hist_E_t -> Fill( (upper_E - m_e)*1000.0 );
		//	hist_E_b -> Fill( (lower_E - m_e)*1000.0);
			hist_costheta -> Fill(the_costheta);
			
	//		prob_v_costheta -> Fill(the_costheta, probability);
			travel_time_v_costheta -> Fill(the_costheta, time_to_travel);
			
			// top
			if(scint_time_t->size()>0)
			{
				bb1_top_position -> Fill( (*bb1_t_x)[0], (*bb1_t_y)[0] );
				if(sigma_plus)
				{
			//		hist_E_tp -> Fill((upper_E)*1000.0);  // for beta spectrum fits.
			//		prob_v_E_less -> Fill( (upper_E - m_e)*1000.0, probability);
				}
				else
				{
			//		hist_E_tm -> Fill((upper_E)*1000.0);
			//		prob_v_E_more -> Fill( (upper_E - m_e)*1000.0, probability);
				}
			}
			// bottom
			if(scint_time_b->size()>0)
			{
				bb1_bottom_position -> Fill( (*bb1_b_x)[0], (*bb1_b_y)[0] );
				if(sigma_plus)
				{
			//		hist_E_bp -> Fill((lower_E)*1000.0);
			//		prob_v_E_more -> Fill( (lower_E - m_e)*1000.0, probability);
				}
				else
				{
			//		hist_E_bm -> Fill((lower_E)*1000.0);
			//		prob_v_E_less -> Fill( (lower_E - m_e)*1000.0, probability);
				}
			}
			*/
	}
	cout << "Done looking through events." << endl;
	
	// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
	// do the kurie plot:
	// Holstein:
	double the_last_energy=-1;
	if( !HistsHaveSameBinning2(hist_E_both_holstein, kurie_plot_holstein) )
	{
		cout << "Bad!  You're done." << endl;
		assert(0);
		return 0;
	}
	double the_p, the_P, the_E, the_dE;
	double D_y, the_N, the_n;
	//
	the_N = hist_E_both_holstein->Integral();
	int nbins = hist_E_both_holstein->GetNbinsX();
	TH1D* normalized_probability = (TH1D*)hist_E_both_holstein->DrawNormalized();
	for (int i=1; i<nbins; i++)  // Bins i=0, i=n_bins are the underflow and overflow?
	{
		the_E = hist_E_both_holstein->GetBinCenter(i);
		the_P = normalized_probability->GetBinContent(i);
		the_n = hist_E_both_holstein->GetBinContent(i);
		the_p = pbeta(the_E);
		
		the_dE= 0.5*hist_E_both_holstein->GetBinWidth(i);
		//
		if(the_n == 0) 
		{ 
			if( the_E>500 && the_last_energy==-1)
			{
			//	cout << "** last energy is now!  i=" << i << ";\tthe_P=" << the_P << ";\tthe_n=" << the_n << endl;
				the_last_energy = hist_E_both_holstein->GetBinCenter(i-1);
			}
			the_n = 1.0; // *only* for the purpose of finding an error.
		} 
		D_y = sqrt( 1.0/ (4.0*the_E*the_p*the_N) ) \
			* sqrt( 1.0 + 0 + ( the_n*pow(the_E, -2)/the_p \
				+ 4.0*the_n*pow(the_E, 2)*pow(the_p, -5) + 4.0*the_n*pow(the_p, -3))*the_dE*the_dE );
		
		if(the_E >= 511)
		{
			kurie_plot_holstein -> SetBinContent(i, sqrt(the_P / (the_p*the_E)) );
			kurie_plot_holstein -> SetBinError(i, D_y);
		}
	}
	cout << "Holstein:  the_last_energy = " << the_last_energy << endl;
	cout << "Holstein:  bin width = " << the_dE << endl;
	
	
	
	// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
	// Fit some fits...
	/*
	TF1 * Abeta_func_keV = make_Abeta_func();
	TFitResultPtr thisfit;
	TH1D* hist_asym = make_asymmetry_histogram(hist_E_tp, hist_E_tm, hist_E_bp, hist_E_bm, string("Asymmetry"), int(kBlack) );
	thisfit = hist_asym -> Fit(Abeta_func_keV, "+S", "same",  500, 4500);
	printresults_A(thisfit);
	*/
	//
	
	TF1 * betaphasespace_func_h = make_betaEphasespace_func();
	betaphasespace_func_h -> SetParameter(0, 5.48989e-09);
	betaphasespace_func_h -> SetParameter(1, 5636.46);
	betaphasespace_func_h -> SetLineColor(kMagenta);
	TFitResultPtr thisfit_holstein_energy;
	thisfit_holstein_energy = hist_E_both_holstein -> Fit(betaphasespace_func_h, "+S", "same", m_e*1000.0, E0*1000.0);
	printresults_betaphasespace(thisfit_holstein_energy);

	TF1 * betaphasespace_func_j = make_betaEphasespace_func();
	betaphasespace_func_j -> SetParameter(0, 5.48989e-09);
	betaphasespace_func_j -> SetParameter(1, 5636.46);
	betaphasespace_func_j -> SetLineColor(kBlue);
	TFitResultPtr thisfit_jtw_energy;
	thisfit_jtw_energy = hist_E_both_jtw -> Fit(betaphasespace_func_j, "+S", "same", m_e*1000.0, E0*1000.0);
	printresults_betaphasespace(thisfit_jtw_energy);
	//
	
	/*
	TF1* linefunc = make_linefunc();
	TFitResultPtr thisfit3;
//	thisfit3 = kurie_plot_holstein -> Fit(linefunc, "+SWL", "same", 520.0, the_last_energy);
	thisfit3 = kurie_plot_holstein -> Fit(linefunc, "+S", "same", 520.0, the_last_energy);  // first bin:  515.45 - 535.52
	printresults_linefunc(thisfit3);
	*/
	//
	
	// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
	TH1D * the_tf1_hist_h = tf1_to_hist_like(betaphasespace_func_h, hist_E_both_holstein);
	TH1D * residuhist_holstein  = get_residuals(hist_E_both_holstein, the_tf1_hist_h);
	residuhist_holstein -> GetListOfFunctions() -> Clear();
	
	TH1D * the_tf1_hist_j = tf1_to_hist_like(betaphasespace_func_j, hist_E_both_jtw);
	TH1D * residuhist_jtw  = get_residuals(hist_E_both_jtw, the_tf1_hist_j);
	residuhist_jtw -> GetListOfFunctions() -> Clear();

	
	TH1D * residuhist2h = (TH1D*)hist_E_both_holstein->Clone("residuhist2h");
	residuhist2h -> GetListOfFunctions() -> Clear();
	residuhist2h -> Add(betaphasespace_func_h, -1.0);
	residuhist2h -> SetTitle("(Counts) - (Normalized Phasespace Fit)");
	residuhist2h -> GetYaxis() -> SetTitle("Counts - Normalized Phasespace Fit");
	residuhist2h -> SetStats(0);
	TH1D* acloneh = (TH1D*)residuhist2h->Clone("acloneh");
	acloneh -> Sumw2(false);
	acloneh -> SetMarkerColor(kRed);
	acloneh -> SetMarkerStyle(20);
	acloneh -> SetMarkerSize(0.7);

	TH1D * residuhist2j = (TH1D*)hist_E_both_jtw->Clone("residuhist2j");
	residuhist2j -> GetListOfFunctions() -> Clear();
	residuhist2j -> Add(betaphasespace_func_j, -1.0);
	residuhist2j -> SetTitle("(Counts) - (Normalized Phasespace Fit)");
	residuhist2j -> GetYaxis() -> SetTitle("Counts - Normalized Phasespace Fit");
	residuhist2j -> SetStats(0);
	TH1D* aclonej = (TH1D*)residuhist2j->Clone("aclonej");
	aclonej -> Sumw2(false);
	aclonej -> SetMarkerColor(kBlue);
	aclonej -> SetMarkerStyle(20);
	aclonej -> SetMarkerSize(0.7);


	TH1D * residuhist3h = (TH1D*)hist_E_both_holstein->Clone("residuhist3h");
	residuhist3h -> GetListOfFunctions() -> Clear();
	residuhist3h -> Divide(betaphasespace_func_h);
	residuhist3h -> SetTitle("(Counts) / (Normalized Phasespace Fit)");
	residuhist3h -> GetYaxis() -> SetTitle("Counts / Normalized Phasespace Fit");
	residuhist3h -> GetYaxis() -> SetRangeUser(0.98, 1.02);
	residuhist3h -> SetStats(0);
	TH1D* a_nothercloneh = (TH1D*)residuhist3h->Clone("a_nothercloneh");
	a_nothercloneh -> Sumw2(false);
	a_nothercloneh -> SetMarkerColor(kRed);
	a_nothercloneh -> SetMarkerStyle(20);
	a_nothercloneh -> SetMarkerSize(0.7);

	TH1D * residuhist3j = (TH1D*)hist_E_both_jtw->Clone("residuhist3j");
	residuhist3j -> GetListOfFunctions() -> Clear();
	residuhist3j -> Divide(betaphasespace_func_j);
	residuhist3j -> SetTitle("(Counts) / (Normalized Phasespace Fit)");
	residuhist3j -> GetYaxis() -> SetTitle("Counts / Normalized Phasespace Fit");
	residuhist3j -> GetYaxis() -> SetRangeUser(0.98, 1.02);
	residuhist3j -> SetStats(0);
	TH1D* a_notherclonej = (TH1D*)residuhist3j->Clone("a_notherclonej");
	a_notherclonej -> Sumw2(false);
	a_notherclonej -> SetMarkerColor(kBlue);
	a_notherclonej -> SetMarkerStyle(20);
	a_notherclonej -> SetMarkerSize(0.7);

	// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
	TText *datalabel = new TText();
	datalabel -> SetNDC();
	datalabel -> SetTextColor(1);
	datalabel -> SetTextSize(0.018);
	
	TText *datalabel2 = new TText();
	datalabel2 -> SetNDC();
	datalabel2 -> SetTextColor(1);
	datalabel2 -> SetTextSize(0.018*2);
	
	gStyle->SetPalette(kDeepSea);
	gStyle->SetOptStat("eiuo");
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerSize(0.7);
	
	// Setup global stat print options?
	
	TLegend * myLegend;
	
	// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
	// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
	
	TCanvas * c12 = new TCanvas("Another Ebeta Phase Space Canvas", "Another Ebeta Phase Space Canvas", 100, 0, 900, 700);
//	hist_E_both_holstein -> SetOptStat("e");
	hist_E_both_holstein -> SetTitle("Beta Total Energy");
	hist_E_both_holstein -> SetStats(0);
	hist_E_both_holstein -> Draw("");
	hist_E_both_jtw -> Draw("same");
	
	myLegend = new TLegend(.80,.82,.98,.98,"");
	myLegend -> SetFillColor(0);
	myLegend -> AddEntry (hist_E_both_holstein, "Holstein","l");  
	myLegend -> AddEntry (hist_E_both_jtw, "JTW","l");  
	myLegend -> Draw();
	gPad->Update();

	TCanvas * c13 = new TCanvas("A ResiduCanvas", "A ResiduCanvas", 100, 0, 900, 700);
	c13 -> Divide(1,3);
	c13 -> cd(1);
//	hist_E_both_holstein -> Draw();
	residuhist_holstein -> GetYaxis() -> SetTitle(residuhist_holstein ->GetTitle() );
	residuhist_holstein -> SetTitle("Holstein");
	residuhist_holstein -> SetStats(0);
	residuhist_holstein -> Draw();
	residuhist_holstein -> DrawCopy("P same");
	datalabel2 -> DrawText(0.14, 0.908, __SHORT_FORM_OF_FILE__);
	c13 -> cd(2);
	residuhist2h -> SetStats(0);
	residuhist2h -> Draw();
	acloneh -> Draw("P same");
	c13 -> cd(3);
	residuhist3h -> Draw();
	a_nothercloneh -> Draw("P same");
	gPad->Update();
	
	///////
//	TH1D* normalized_holstein = (TH1D*)hist_E_both_holstein->DrawNormalized();
//	TH1D* normalized_jtw      = (TH1D*)hist_E_both_jtw->DrawNormalized();
	TH1D* normalized_holstein = manual_normalize(hist_E_both_holstein);
	normalized_holstein->SetLineColor(kBlack);
	normalized_holstein->SetMarkerColor(kBlack);
	normalized_holstein->SetMarkerStyle(20);
	normalized_holstein->SetMarkerSize(0.7);
	TH1D* normalized_jtw      = manual_normalize(hist_E_both_jtw);
	
	new TCanvas();
	normalized_holstein -> Draw();
	new TCanvas();
	normalized_jtw -> Draw();
	
	TH1D * residuhist2_subtract = (TH1D*)normalized_holstein->Clone("residuhist2_subtract");
//	residuhist2_subtract -> SetLineColor(kBlack);
	residuhist2_subtract -> SetMarkerColor(kGreen);
//	residuhist2_subtract -> GetListOfFunctions() -> Clear();
	residuhist2_subtract -> Add(normalized_jtw, -1.0);
	residuhist2_subtract -> SetTitle("(Holstein) - (JTW)");
	residuhist2_subtract -> GetYaxis() -> SetTitle("(Holstein) - (JTW)");
	residuhist2_subtract -> SetStats(0);


	TCanvas * c14 = new TCanvas("Holstein-JTW Compare Canvas", "Holstein-JTW Compare Canvas", 100, 0, 900, 700);
	residuhist2_subtract -> Draw("P");
	gPad->Update();
	
	TH1D * residuhist3_divide = (TH1D*)normalized_holstein->Clone("residuhist3_divide");
//	residuhist3_divide -> SetLineColor(kBlack);
	residuhist3_divide -> SetMarkerColor(kMagenta);
//	residuhist3_divide -> GetListOfFunctions() -> Clear();
	residuhist3_divide -> Divide(normalized_jtw);
	residuhist3_divide -> SetTitle("(Holstein) / (JTW)");
	residuhist3_divide -> GetYaxis() -> SetTitle("(Holstein) / (JTW)");
	residuhist3_divide -> GetYaxis() -> SetRangeUser(0.98, 1.02);
	residuhist3_divide -> SetStats(0);

	TCanvas * c15 = new TCanvas("Holstein-JTW Compare Canvas2", "Holstein-JTW Compare Canvas2", 100, 0, 900, 700);
	residuhist3_divide -> Draw("P");
	gPad->Update();

	// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
	/*
	// leave the canvas structure and hists alone, but don't pop it up anymore.
	TCanvas * c = new TCanvas("E_beta Canvas", "E_beta Canvas", 100.0, 0.0, 900, 700);
	c -> Divide(1,2);
	c -> cd(1);
	hist_E_t -> Draw();
	hist_E_t -> SetStats(0);
	gPad->BuildLegend(.76,.90,.98,.98,"");
	datalabel2 -> DrawText(0.10, 0.908, __SHORT_FORM_OF_FILE__);
	datalabel2 -> DrawText(0.10, 0.938, "Simple Simulation E_beta" );
	c -> cd(2);
	hist_E_b -> Draw();
	hist_E_b -> SetStats(0);
	gPad->BuildLegend(.76,.90,.98,.98,"");
	gPad->Update();
	*/
	/*
	// leave the bb1 hit position histograms, but don't pop up the canvas.
	TCanvas * c3 = new TCanvas("BB1 Canvas", "BB1 Canvas", 100.0, 0.0, 900, 450);
	c3 -> Divide(2,1);
	c3 -> cd(1);
	bb1_top_position -> Draw("colz");
	datalabel -> DrawText(0.10, 0.908, __SHORT_FORM_OF_FILE__);
	datalabel -> DrawText(0.10, 0.928, "Simple Simulation A_beta" );
	c3 -> cd(2);
	bb1_bottom_position -> Draw("colz");
	gPad->Update();
	*/
//	TCanvas * c7 = new TCanvas("costheta v E Canvas", "costheta v E Canvas", 100.0, 0.0, 900, 700);
//	costheta_v_E -> Draw("colz");
//	gPad->Update();
	
	// leave travel_time_v_costheta, but don't pop up.
//	TCanvas * c8 = new TCanvas("time costheta canvas", "time costheta canvas", 100, 0, 900, 700);
//	travel_time_v_costheta -> Draw("colz");
//	gPad -> Update();
	//
		
	/*
	TCanvas * c1 = new TCanvas("E_beta Canvas2", "E_beta Canvas2", 100.0, 0.0, 900, 700);
	c1 -> Divide(2,2);
	c1 -> cd(1);
	gPad -> SetLogy();
	hist_E_tp -> Draw();
	datalabel2 -> DrawText(0.10, 0.908, __SHORT_FORM_OF_FILE__);
	datalabel2 -> DrawText(0.10, 0.938, "Simple Simulation E_beta" );
	c1 -> cd(2);
	gPad -> SetLogy();
	hist_E_tm -> Draw();
	c1 -> cd(3);
	gPad -> SetLogy();
	hist_E_bp -> Draw();
	c1 -> cd(4);
	gPad -> SetLogy();
	hist_E_bm -> Draw();
//	hist_E_bm -> SetStats(0);
	gPad->Update();
	*/
	
	//
	/*
	TCanvas * c2 = new TCanvas("Asymmetry Canvas", "Asymmetry Canvas", 100.0, 0.0, 900, 700);
	hist_asym -> Draw();
	hist_asym -> SetStats(0);
	datalabel -> DrawText(0.10, 0.908, __SHORT_FORM_OF_FILE__);
	datalabel -> DrawText(0.10, 0.928, "Simple Simulation A_beta" );
	gPad->Update();
	*/
	
	/*
	TCanvas * c4 = new TCanvas("pdf Canvas", "pdf Canvas", 100.0, 0.0, 900, 700);
	c4->Divide(1,2);
	c4->cd(1);
	prob_v_E_more -> Draw("colz");
	datalabel2 -> DrawText(0.10, 0.908, __SHORT_FORM_OF_FILE__);
	datalabel2 -> DrawText(0.10, 0.938, "Simple Simulation E_beta" );
	c4->cd(2);
	prob_v_E_less -> Draw("colz");
	gPad->Update();
	*/
	
	/*
	TCanvas * c6 = new TCanvas("costheta Canvas", "costheta Canvas", 100.0, 0.0, 900, 700);
	hist_costheta -> SetStats("e");
	hist_costheta -> Draw();
	datalabel -> DrawText(0.14, 0.908, __SHORT_FORM_OF_FILE__);
	datalabel -> DrawText(0.14, 0.928, "Simple Simulation A_beta" );
	gPad->Update();
	*/
	
	/*
	TCanvas * c16 = new TCanvas("Other costheta Canvas", "Other costheta Canvas", 100.0, 0.0, 900, 700);
	prob_v_costheta -> Draw();
	datalabel -> DrawText(0.10, 0.908, __SHORT_FORM_OF_FILE__);
	gPad->Update();
	*/
	
	/*
	TCanvas * c11 = new TCanvas("Kurie Plot Canvas", "Kurie Plot Canvas", 100, 0, 900, 700);
	kurie_plot_holstein -> SetStats(0);
	kurie_plot_holstein -> Draw();
//	datalabel -> DrawText(0.10, 0.908, __SHORT_FORM_OF_FILE__);
//	datalabel -> DrawText(0.10, 0.928, "Simple Simulation A_beta" );
	gPad->Update();
	*/
	/*
	TCanvas * rescan = new TCanvas("Residuals Canvas", "Residuals Canvas", 100.0, 0.0,900,700);
	vector<TPad *> rescan_residupad = make_residupad(hist_E_both_holstein, residuhist, string("E1") );
	rescan_residupad.at(0) -> cd();
	rescan_residupad.at(0) -> SetGridx();
	residuhist -> Draw("hist same");
	hist_E_both_holstein -> GetYaxis()->SetTitle("Counts");
	hist_E_both_holstein -> GetXaxis()->SetNdivisions(10);  // 
	gPad->Update();
//	datalabelr -> DrawText(0.10, 0.908, __SHORT_FORM_OF_FILE__);
//	gPad -> BuildLegend(.76,.88,.98,.98,"");
	gPad -> Modified();
	hist_E_both_holstein -> Draw("same");
	gPad -> Update();
	*/
	
	/*
	TCanvas * rescan2 = new TCanvas("Residuals Canvas2", "Residuals Canvas2", 100.0, 0.0,900,700);
	vector<TPad *> rescan_residupad2 = make_residupad(hist_E_both_holstein, residuhist, string("E1") );
	rescan_residupad2.at(0) -> cd();
	rescan_residupad2.at(0) -> SetGridx();
	residuhist2 -> Draw("hist same");
	hist_E_both_holstein -> GetXaxis()->SetNdivisions(10);  // 
	gPad->Update();
//	datalabelr -> DrawText(0.10, 0.908, __SHORT_FORM_OF_FILE__);
	gPad -> BuildLegend(.76,.88,.98,.98,"");
	gPad -> Modified();
	hist_E_both_holstein -> Draw("same");
	gPad -> Update();
	*/
	
	
	rootapp->Run();
	
	return 0;
}











