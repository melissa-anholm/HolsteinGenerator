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
	
	string filename1 = "output_withroc.root";
	string filename2 = "output_noroc.root";
	TFile * f = new TFile(filename1.c_str(), "READ");
	cout << "Reading from file:  " << filename1 << endl;
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
	
	Bool_t pdf_acceptance = kFALSE;
	Bool_t det_acceptance = kFALSE;
	tree -> SetBranchAddress("pdf_acceptance", &pdf_acceptance);
	tree -> SetBranchAddress("det_acceptance", &det_acceptance);
	
	int kurie_rebin = 16*2;
	TH1D * hist_E_both_withroc = CreateHist(string("Beta Total Energy"), string("Ben_Ebeta_4096"), int(kBlack), kurie_rebin);
	hist_E_both_withroc -> GetXaxis() -> SetTitle("Beta Energy [keV]");
	hist_E_both_withroc -> GetXaxis() -> SetRangeUser(0.0, 6000.0);
	hist_E_both_withroc -> Sumw2();
	
	TH1D * kurie_plot_withroc = CreateHist(string("Kurie Plot"), string("Ben_Ebeta_4096"), int(kBlack), kurie_rebin);
	kurie_plot_withroc -> GetXaxis() -> SetTitle("Beta Energy [keV]");
	kurie_plot_withroc -> GetYaxis() -> SetTitle("sqrt( P(E) / (p*E) ) ");
	kurie_plot_withroc -> GetXaxis() -> SetRangeUser(0.0, 6000.0);
	kurie_plot_withroc -> Sumw2();
	
	// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
	// Look through events.
	Long64_t nentries = tree->GetEntries();
	double rhit = -1.0;
	bool accept=false;
	cout << "nentries = " << nentries << endl;
	for(int i=0; i<nentries; i++)
	{
		if( (i % 100000) == 0) { cout<<"Reached entry "<< i << endl; }
		tree -> GetEntry(i);
		
		if(scint_time_t->size()>0 && scint_time_b->size() > 0 )
			{ cout << "That's bad.  i = "     << i << endl;  assert(0);  return 0; }
		else if(scint_time_t->size()==0 && scint_time_b->size()==0 )
			{ cout << "That's bad too.  i = " << i << endl;  assert(0);  return 0; }
		
		//
		if(scint_time_t->size()>0) 
			{ rhit = sqrt( pow((*bb1_t_x)[0], 2) + pow((*bb1_t_y)[0], 2)); }
		if(scint_time_b->size()>0)
			{ rhit = sqrt( pow((*bb1_b_x)[0], 2) + pow((*bb1_b_y)[0], 2)); }
		
	//	accept=true;  // all events.
	//	accept = pdf_acceptance;
	//	accept = det_acceptance;
		accept = (pdf_acceptance&&det_acceptance);
		if(rhit>15.5) { accept=false; }
		if(accept)
		{
			hist_E_both_withroc -> Fill( (upper_E+lower_E)*1000.0 );
			// top
			if(scint_time_t->size()>0)
			{
			}
			// bottom
			if(scint_time_b->size()>0)
			{
			}
		}
	}
	cout << "Done looking through events." << endl;
	
	// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
	// do the kurie plot:
	double the_last_energy=-1;
	if( !HistsHaveSameBinning2(hist_E_both_withroc, kurie_plot_withroc) )
	{
		cout << "Bad!  You're done." << endl;
		assert(0);
		return 0;
	}
	double the_p, the_P, the_E, the_dE;
	double D_y, the_N, the_n;
	int nbins;
	//
	the_N = hist_E_both_withroc->Integral();
	nbins = hist_E_both_withroc->GetNbinsX();
	TH1D* normalized_probability_withroc = (TH1D*)hist_E_both_withroc->DrawNormalized();
	for (int i=1; i<nbins; i++)  // Bins i=0, i=n_bins are the underflow and overflow?
	{
		the_E = hist_E_both_withroc->GetBinCenter(i);
		the_P = normalized_probability_withroc->GetBinContent(i);
		the_n = hist_E_both_withroc->GetBinContent(i);
		the_p = pbeta(the_E);
		
		the_dE= 0.5*hist_E_both_withroc->GetBinWidth(i);
		//
		if(the_n == 0) 
		{ 
			if( the_E>500 && the_last_energy==-1)
			{
			//	cout << "** last energy is now!  i=" << i << ";\tthe_P=" << the_P << ";\tthe_n=" << the_n << endl;
				the_last_energy = hist_E_both_withroc->GetBinCenter(i-1);
			}
			the_n = 1.0; // *only* for the purpose of finding an error.
		} 
		D_y = sqrt( 1.0/ (4.0*the_E*the_p*the_N) ) \
			* sqrt( 1.0 + 0 + ( the_n*pow(the_E, -2)/the_p \
				+ 4.0*the_n*pow(the_E, 2)*pow(the_p, -5) + 4.0*the_n*pow(the_p, -3))*the_dE*the_dE );
		
		if(the_E >= 511)
		{
			kurie_plot_withroc -> SetBinContent(i, sqrt(the_P / (the_p*the_E)) );
			kurie_plot_withroc -> SetBinError(i, D_y);
		}
	}
	cout << "the_last_energy = " << the_last_energy << endl;
	cout << "bin width = " << the_dE << endl;
	// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
	// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
	// ok, next file...
	
	
	f = new TFile(filename2.c_str(), "READ");
	cout << "Reading from file:  " << filename2 << endl;
	f -> cd();
	tree = (TTree*)f->Get("ntuple");	
	
//	vector<double> * bb1_t_x = 0;
//	vector<double> * bb1_t_y = 0;
//	vector<double> * bb1_b_x = 0;
//	vector<double> * bb1_b_y = 0;
	tree -> SetBranchAddress("bb1_top_x",   &bb1_t_x);
	tree -> SetBranchAddress("bb1_top_y",   &bb1_t_y);
	tree -> SetBranchAddress("bb1_bottom_x",&bb1_b_x);
	tree -> SetBranchAddress("bb1_bottom_y",&bb1_b_y);
//	Double_t upper_E;
//	Double_t lower_E;
	tree -> SetBranchAddress("upper_scint_E", &upper_E);
	tree -> SetBranchAddress("lower_scint_E", &lower_E);
//	Int_t sigma_plus = -1;
	tree -> SetBranchAddress("TTLBit_SigmaPlus", &sigma_plus);  
//	vector<double> * scint_time_t = 0;  
//	vector<double> * scint_time_b = 0;  
	tree -> SetBranchAddress("TDC_SCINT_TOP",    &scint_time_t);  
	tree -> SetBranchAddress("TDC_SCINT_BOTTOM", &scint_time_b);  
	
//	Bool_t pdf_acceptance = kFALSE;
//	Bool_t det_acceptance = kFALSE;
	tree -> SetBranchAddress("pdf_acceptance", &pdf_acceptance);
	tree -> SetBranchAddress("det_acceptance", &det_acceptance);
	
//	int kurie_rebin = 16*2;
	TH1D * hist_E_both_noroc = CreateHist(string("Beta Total Energy (no roc)"), string("Ben_Ebeta_4096"), int(kBlack), kurie_rebin);
	hist_E_both_noroc -> GetXaxis() -> SetTitle("Beta Energy [keV]");
	hist_E_both_noroc -> GetXaxis() -> SetRangeUser(0.0, 6000.0);
	hist_E_both_noroc -> Sumw2();
	
	/*
	TH1D * kurie_plot_noroc = CreateHist(string("Kurie Plot"), string("Ben_Ebeta_4096"), int(kBlack), kurie_rebin);
	kurie_plot_noroc -> GetXaxis() -> SetTitle("Beta Energy [keV]");
	kurie_plot_noroc -> GetYaxis() -> SetTitle("sqrt( P(E) / (p*E) ) ");
	kurie_plot_noroc -> GetXaxis() -> SetRangeUser(0.0, 6000.0);
	kurie_plot_noroc -> Sumw2();
	*/
	
	// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
	// Look through events.
	nentries = tree->GetEntries();
	rhit = -1.0;
	accept=false;
	cout << "nentries = " << nentries << endl;
	for(int i=0; i<nentries; i++)
	{
		if( (i % 100000) == 0) { cout<<"Reached entry "<< i << endl; }
		tree -> GetEntry(i);
		
		if(scint_time_t->size()>0 && scint_time_b->size() > 0 )
			{ cout << "That's bad.  i = "     << i << endl;  assert(0);  return 0; }
		else if(scint_time_t->size()==0 && scint_time_b->size()==0 )
			{ cout << "That's bad too.  i = " << i << endl;  assert(0);  return 0; }
		
		//
		if(scint_time_t->size()>0) 
			{ rhit = sqrt( pow((*bb1_t_x)[0], 2) + pow((*bb1_t_y)[0], 2)); }
		if(scint_time_b->size()>0)
			{ rhit = sqrt( pow((*bb1_b_x)[0], 2) + pow((*bb1_b_y)[0], 2)); }
		
	//	accept=true;  // all events.
	//	accept = pdf_acceptance;
	//	accept = det_acceptance;
		accept = (pdf_acceptance&&det_acceptance);
		if(rhit>15.5) { accept=false; }
		if(accept)
		{
			hist_E_both_noroc -> Fill( (upper_E+lower_E)*1000.0 );
		}
	}
	cout << "Done looking through events." << endl;

	
	// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
	/*
	// Fit some fits...
	TF1 * Abeta_func_keV = make_Abeta_func();
	TFitResultPtr thisfit;
	TH1D* hist_asym = make_asymmetry_histogram(hist_E_tp, hist_E_tm, hist_E_bp, hist_E_bm, string("Asymmetry"), int(kBlack) );
	thisfit = hist_asym -> Fit(Abeta_func_keV, "+S", "same",  500, 4500);
	printresults_A(thisfit);
	//
	*/
	
	TF1 * betaphasespace_func = make_betaEphasespace_func();
	betaphasespace_func -> SetParameter(0, 5.48989e-09);
//	betaphasespace_func -> FixParameter(1, 5636.46);
	betaphasespace_func -> SetParameter(1, 5636.46);
	TFitResultPtr thisfit2;
	thisfit2 = hist_E_both_withroc -> Fit(betaphasespace_func, "+S", "same", m_e*1000.0, E0*1000.0);
	printresults_betaphasespace(thisfit2);
	//
	
	TF1* linefunc = make_linefunc();
	TFitResultPtr thisfit3;
//	thisfit3 = kurie_plot_withroc -> Fit(linefunc, "+SWL", "same", 520.0, the_last_energy);
	thisfit3 = kurie_plot_withroc -> Fit(linefunc, "+S", "same", 520.0, the_last_energy);  // first bin:  515.45 - 535.52
	printresults_linefunc(thisfit3);
	//
	
	// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
	TH1D * the_tf1_hist = tf1_to_hist_like(betaphasespace_func, hist_E_both_withroc);
	TH1D * residuhist_withroc  = get_residuals(hist_E_both_withroc, the_tf1_hist);
	residuhist_withroc -> GetListOfFunctions() -> Clear();
	
	TH1D * residuhist2_withroc = (TH1D*)hist_E_both_withroc->Clone("residuhist2");
	residuhist2_withroc -> GetListOfFunctions() -> Clear();
	residuhist2_withroc -> Add(betaphasespace_func, -1.0);
	residuhist2_withroc -> SetTitle("(Counts) - (Normalized Phasespace Fit)");
	residuhist2_withroc -> GetYaxis() -> SetTitle("Counts - Normalized Phasespace Fit");
	residuhist2_withroc -> SetStats(0);
	TH1D* aclone_withroc = (TH1D*)residuhist2_withroc->Clone("aclone");
	aclone_withroc -> Sumw2(false);
	aclone_withroc -> SetMarkerColor(kRed);
	aclone_withroc -> SetMarkerStyle(20);
	aclone_withroc -> SetMarkerSize(0.7);

	TH1D * residuhist3_withroc = (TH1D*)hist_E_both_withroc->Clone("residuhist3");
	residuhist3_withroc -> GetListOfFunctions() -> Clear();
	residuhist3_withroc -> Divide(betaphasespace_func);
	residuhist3_withroc -> SetTitle("(Counts) / (Normalized Phasespace Fit)");
	residuhist3_withroc -> GetYaxis() -> SetTitle("Counts / Normalized Phasespace Fit");
	residuhist3_withroc -> GetYaxis() -> SetRangeUser(0.98, 1.02);
	residuhist3_withroc -> SetStats(0);
	TH1D* a_notherclone_withroc = (TH1D*)residuhist3_withroc->Clone("a_notherclone");
	a_notherclone_withroc -> Sumw2(false);
	a_notherclone_withroc -> SetMarkerColor(kRed);
	a_notherclone_withroc -> SetMarkerStyle(20);
	a_notherclone_withroc -> SetMarkerSize(0.7);
	
	
	TH1D * diffhist = (TH1D*)hist_E_both_withroc->Clone("diffhist");
	diffhist->SetTitle("Beta Energy:  (ROC) - (No ROC)");
	diffhist->SetName("Beta Energy:  (ROC) - (No ROC)");
	diffhist->Add(hist_E_both_noroc, -1.0);
	diffhist -> GetListOfFunctions() -> Clear();
	
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
	// Setup global stat print options?
	
	
	// -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // -- * -- // 
	
	new TCanvas("name", "name", 100, 0, 900, 700);
	diffhist -> Draw();
	gPad->Update();
	
	
	TCanvas * c12 = new TCanvas("Another Ebeta Phase Space Canvas", "Another Ebeta Phase Space Canvas", 100, 0, 900, 700);
	hist_E_both_withroc -> Draw();
	gPad->Update();
	
	TCanvas * c11 = new TCanvas("Kurie Plot Canvas", "Kurie Plot Canvas", 100, 0, 900, 700);
	kurie_plot_withroc -> SetStats(0);
	kurie_plot_withroc -> Draw();
//	datalabel -> DrawText(0.10, 0.908, __SHORT_FORM_OF_FILE__);
//	datalabel -> DrawText(0.10, 0.928, "Simple Simulation A_beta" );
	gPad->Update();
	
	TCanvas * c13 = new TCanvas("A ResiduCanvas", "A ResiduCanvas", 100, 0, 900, 700);
	c13 -> Divide(1,3);
	c13 -> cd(1);
	hist_E_both_withroc -> Draw();
	hist_E_both_withroc -> SetStats(0);
	datalabel2 -> DrawText(0.14, 0.908, __SHORT_FORM_OF_FILE__);
	c13 -> cd(2);
	residuhist2_withroc -> SetStats(0);
	residuhist2_withroc -> Draw();
	aclone_withroc -> Draw("P same");
	c13 -> cd(3);
	residuhist3_withroc -> Draw();
	a_notherclone_withroc -> Draw("P same");
	gPad->Update();
	
	
	TCanvas * rescan = new TCanvas("Residuals Canvas", "Residuals Canvas", 100.0, 0.0,900,700);
	vector<TPad *> rescan_residupad = make_residupad(hist_E_both_withroc, residuhist_withroc, string("E1") );
	rescan_residupad.at(0) -> cd();
	rescan_residupad.at(0) -> SetGridx();
	residuhist_withroc -> Draw("hist same");
	hist_E_both_withroc -> GetXaxis()->SetNdivisions(10);  // 
	gPad->Update();
//	datalabelr -> DrawText(0.10, 0.908, __SHORT_FORM_OF_FILE__);
	gPad -> BuildLegend(.76,.88,.98,.98,"");
	gPad -> Modified();
	hist_E_both_withroc -> Draw("same");
	gPad -> Update();
	
	
	rootapp->Run();
	
	return 0;
}











