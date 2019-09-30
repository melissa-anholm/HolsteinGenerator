// Code by Melissa Anholm
// Feb 2019 - 
// This code is intended to work with the 98% branch of 37K decay.
// this is beta+ decay, and both the parent and daughter have spin 3/2.
// this information is likely baked in somewhere.


#include <iostream>  // cout, endl
#include <cmath>     // pow
#include "fstream"   // 
#include <iomanip>   // setw


#include <sys/time.h>

//#include <map>     // for text input function

#undef NDEBUG
#include<assert.h>
#include<iterator>

#include<Randomize.hh>
#include<G4ThreeVector.hh> // probably the correct way to include ThreeVector.h.
// #include <CLHEP/Units/SystemOfUnits.h> in HolsteinVars.cpp.  actually G4SystemOfUnits.hh
// #include<G4SystemOfUnits.hh> in HolsteinVars.cpp.  *That* includes <CLHEP/Units/SystemOfUnits.h>.

#include "TFile.h"
#include <TTree.h>
#include <TBranch.h> // might not need this...

//
#include "HolsteinVars.cpp"

// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- //
	
//  HolsteinDecay class relies on other classes: 
//  	HolsteinVars
//  		isotope_values
//  		SS (splitstring)
//  	sublevel_populations
class HolsteinDecay
{
public:
	HolsteinDecay();    /// this one seems broken?!?!?
	HolsteinDecay(HolsteinVars* HV, sublevel_populations * pops);
//	HolsteinDecay(HolsteinVars  HV) { HolsteinDecay( (HolsteinVars*)(&HV) ); };
	
	// The cone!
	void set_use_cone(bool useit) { use_cone = useit; };
	bool get_use_cone() { return use_cone; };
	bool use_roc;  // don't even bother protecting it yet...
	
	// probably should put some sanity checks into functions that set cone_costheta...
	void set_conecostheta(double costheta_max_)  { cone_costheta = costheta_max_; }; 
//	void set_conetheta(double theta_max_) { set_conecostheta( cos(theta_max_) ); }; 
	double get_conecostheta() { return cone_costheta; };
	
	void randomize_direction();
	
	bool check_PDF_acceptance();         // uses initial_momentum
	bool check_detector_acceptance();    // creates hit_position from initial_momentum and initial_position.  then checks.
	void randomize_nuclear(bool doit=true);
	void randomize_atomic(bool doit=true);
	void randomize_start(bool doit=true);   // initial_position, in G4 mm.
	bool shoot_decayevent();
	void print_results();
	void print_vars() { this->Params->print_vars(); };
	
	// specific to the individual decay:
	G4double Ebeta;
	// Variables for export only...
	double Ebeta_tot_MeV;
	double pbeta_MeV;
	double Ebeta_kin_MeV;
	double vbeta_over_c;
	
	double costheta_lab;
	double time_to_travel;
	double the_probability; // obsolete.
	double jtw_probability, holstein_probability; // specific to this value of Ebeta, costheta.
	
	G4ThreeVector initial_momentum;  // needs units of energy.
	G4ThreeVector initial_velocity;  // needs units.  not unitless.
	G4ThreeVector initial_position;  // needs units of position.  ...done in G4units.
	G4ThreeVector hit_position;      // needs units of position.
	bool pdf_acceptance; // obsolete.
	bool jtw_acceptance, holstein_acceptance;
	bool det_acceptance;
	bool going_up;
	void SetAcceptanceMode(string);
	
	
	G4double pbeta(G4double E);  
	G4double get_Ebeta(G4double pbeta);
	double get_v_from_p(G4double pbeta);  
	G4double get_p_from_v(double vbeta);
	
	double FermiFunction(double Z_, G4double E) { return 1.0; };
	
	sublevel_populations * the_pops;
	chamber_geometry the_geometry;

	double Mz, Mz2, Mz3;
	void run_fast(bool doit) { runfast=doit; };
	
	// For export only:
	double jtw_xi, jtw_Abeta, jtw_rho; // not a function of E.  const. for a given set of parameters.
	double holstein_Abeta; // event-specific.  function of E.

private:
	// run parameters..
	bool use_cone;
	double cone_costheta;
	bool runfast; // don't look at detector acceptance if it's not accepted to the PDF.  Should check detector acceptance first, maybe...
//	bool 
	double prob_max;
	
	//, event_P, event_T, event_Mzcubed;
	double M_F_isospin, I_spin, A_nucleons, Z_daughter;  // these are exact, and we're just copying the values over because poor coding practice but better readability.
	double hbarc_eV_nm;  // this value of hbarc has its units baked right in.  Because I want matrix elements to come out unitless.  
	
	G4double E0, M, Delta, m_e;  // G4units propagated for M, Delta, m_e.  E0(??)
	double rc; // unitless.
	double speed_of_light;  // no G4units.  converts to m/s.
	
	double deltaC;
	double quad_parent, quad_daughter;
	double mu_parent, mu_daughter;
	double g_V, g_A, g_II, g_S, g_P, g_M;
	
	double M_F, M_GT;
	double M_r2, M_sr2, M_sL, M_srp, M_rdotp, M_1y, M_2y, M_3y; // "exact".
	double M_Q, M_rp;
	
	G4double delta_M1, delta_M2;  // comes out in MeV from HolsteinVars' sigma_M1, sigma_M2.
	double delta_deltaC;
	double delta_mu_parent, delta_mu_daughter;
	double delta_quad_parent, delta_quad_daughter;
	
	double a1, a2, c1, c2;
	double b, d, e, g, f, h;
	double j2, j3;
	
	double jtw_a1, jtw_c1;
	
	// Need to evaluate the prob. dist. from within the daughter class..
	double F_0(G4double);  // monopole
	double F_1(G4double);  // dipole
	double F_2(G4double);  // quadrupole
	double F_3(G4double);  // octopole
	
	double Lambda0, Lambda1, Lambda2, Lambda3, Lambda4;
	void initialize_lambdafuncs();
	double get_lambda0();
	double get_lambda1();
	double get_lambda2();
	double get_lambda3();
	double get_lambda4();
	
	double get_probability(G4double, double);  // E, costheta
	double get_jtw_probability(G4double, double);
	
	HolsteinVars * Params;
//	TTree * the_tree;
};

void HolsteinDecay::SetAcceptanceMode(string the_mode)
{
	cout << "SetAcceptanceMode(...) doesn't work yet." << endl;
	return;
}


HolsteinDecay::HolsteinDecay() // broken?
{
	cout << "Initializing HolsteinDecay without input parameters.  " << endl;
	cout << "Creating a new set of HolsteinVars, and a new set of sublevel_populations.  " << endl;
	cout << "I think this is broken and will cause segfaults when we try to look for the HolsteinVars or the sublevel_populations." << endl;
	cout << "Hard kill it now." << endl;
	assert(0);
	
	//
	/*
	HolsteinVars * HV           = new HolsteinVars();
	sublevel_populations * pops = new sublevel_populations(1);
	cout << "things are created." << endl;
	HolsteinDecay( (HolsteinVars*)HV, (sublevel_populations*)pops );
	*/
	Params   = new HolsteinVars();
	the_pops = new sublevel_populations(1);
	HolsteinDecay( (HolsteinVars*)Params, (sublevel_populations*)the_pops );
}

HolsteinDecay::HolsteinDecay(HolsteinVars * HV, sublevel_populations * pops):
	// nuclear parameters:
	delta_M1(0), delta_M2(0), delta_deltaC(0), 
	delta_mu_parent(0), delta_mu_daughter(0), delta_quad_parent(0), delta_quad_daughter(0), 
	// run parameters:
	use_roc(true),
	use_cone(false), cone_costheta(0),
	runfast(true), prob_max(0.2)
{
	// Nuclear:
	Params = HV;
	
	M_F_isospin    = Params->M_F_isospin;    // exact
	I_spin         = Params->I_spin;         // exact
	A_nucleons     = Params->A_nucleons;     // exact
	hbarc_eV_nm          = Params->hbarc_eV_nm;          // exact.  unitless. ... it is *not* fucking unitless.  ... it needs to be unitless for now though, because who fucking knows what units the matrix elements are supposed to have.
	Z_daughter     = Params->Z_daughter;     // exact
	speed_of_light = Params->speed_of_light; // exact.  double.  no G4units.
	
	E0    = Params->E0;    // uncertain (propagation) MeV.  comes out in MeV from HolsteinVars.
	M     = Params->M;     // uncertain (propagation) MeV.  comes out in MeV from HolsteinVars.
	Delta = Params->Delta; // uncertain (propagation) MeV.  comes out in MeV from HolsteinVars.  units propagated.
	m_e   = Params->m_e;   // uncertain (direct).     MeV.  comes out in MeV from HolsteinVars.
	
	deltaC = Params->deltaC; // uncertain (direct)
	
	quad_parent   = Params->quad_parent;   // uncertain (direct) 
	quad_daughter = Params->quad_daughter; // uncertain (direct) 
	mu_parent     = Params->mu_parent;     // uncertain (direct) 
	mu_daughter   = Params->mu_daughter;   // uncertain (direct) 
	
	g_V  = Params->g_V;  // exact
	g_A  = Params->g_A;  // exact-ish??
	g_II = Params->g_II; // exact
	g_S  = Params->g_S;  // exact
	g_P  = Params->g_P;  // exact
	g_M  = Params->g_M;  // exact
	
	M_F     = Params->M_F;     // uncertain (propagation, deltaC) 
	M_GT    = Params->M_GT;    // uncertain (propagation, deltaC) 
	M_r2    = Params->M_r2;    // exact
	M_sr2   = Params->M_sr2;   // exact
	M_sL    = Params->M_sL;    // exact
	M_srp   = Params->M_srp;   // exact
	M_rdotp = Params->M_rdotp; // exact
	M_1y    = Params->M_1y;    // exact
	M_2y    = Params->M_2y;    // exact
	M_3y    = Params->M_3y;    // exact
	M_Q     = Params->M_Q;     // uncertain (propagation) -- I think this doesn't get used directly...
	M_rp    = Params->M_rp;    // uncertain (propagation)
	// M_L = Params->M_L;  // exact
	
	a1 = Params->a1; // uncertain (propagation, M1 M2, M_F)
	a2 = Params->a2; // uncertain (propagation, M1 M2)
	c1 = Params->c1; // uncertain (propagation, M1 M2, M_GT, (g_A) )
	c2 = Params->c2; // uncertain (propagation, M1 M2, (g_A) )
	
	b  = Params->b;  // uncertain (propagation, mu_parent mu_daughter)
	d  = Params->d;  // uncertain (propagation, M1 M2, M_GT, (g_A) )
	e  = Params->e;  // exact by isospin(?)+CVC.
	g  = Params->g;  // uncertain (propagation, M1 M2, quad_parent, quad_daughter)
	f  = Params->f;  // uncertain (propagation, M1 M2, quad_parent, quad_daughter)
	h  = Params->h;  // uncertain (propagation, M1 M2, M_GT, (g_A) )
	
	j2 = Params->j2; // uncertain (propagation, M1 M2, (g_A) )
	j3 = Params->j3; // uncertain (propagation, M1 M2, (g_A) )
	
	jtw_a1 = g_V*Params->FindValue("M_F");  // hasn't had isospin correction turned on yet.
	jtw_c1 = g_A*Params->FindValue("M_GT"); // hasn't had isospin correction turned on yet.
	jtw_xi    = 2.0*(jtw_a1*jtw_a1 + jtw_c1*jtw_c1);
	jtw_Abeta = 2.0*(jtw_c1*jtw_c1*(1.0/(3.0/2.0 + 1.0)) + 2.0*jtw_a1*jtw_c1*sqrt((3.0/2.0)/(3.0/2.0 + 1.0)) ) / (2.0*(jtw_a1*jtw_a1 + jtw_c1*jtw_c1));
	jtw_rho   = jtw_c1/jtw_a1;
	
	// Atomic:
	// by default, the_pops is set to be fully polarized sigma+.
	the_pops = pops;
	Mz  = the_pops->get_Mz();   // 1.5
	Mz2 = the_pops->get_Mz2();  // 2.25
	Mz3 = the_pops->get_Mz3();  // 3.375
	
	initialize_lambdafuncs();
	
	// Decay-specific Parameters:
	initial_momentum = G4ThreeVector(0.0*MeV,0.0*MeV,0.0*MeV);  // 
	
}


bool HolsteinDecay::shoot_decayevent() // uses param. runfast to decide whether to check *both* acceptances if one has already failed.
{
	bool verbose = false;
	pdf_acceptance = false;  // obsolete?
	jtw_acceptance = false;
	holstein_acceptance = false;
	
	det_acceptance = false;
	
	this -> randomize_nuclear(false);
	this -> randomize_atomic(false);
	this -> randomize_direction();           // pick Ebeta, costheta (possibly from the cone.)
	
	if(runfast)  // true by default, but I'm pretty sure I disable it in main(...).  this is useful to turn off if I want to record *all* the events generated and check acceptance later.
	{
		this -> check_detector_acceptance(); // set:  det_acceptance.
		if(det_acceptance)
		{
			this -> check_PDF_acceptance();  // set:  the_probability, pdf_acceptance.  jtw_probability, jtw_acceptance.  holstein_probability, holstein_acceptance.
		}
	}
	else
	{
		this -> check_PDF_acceptance();      // set:  the_probability, pdf_acceptance.
		this -> check_detector_acceptance(); // set:  det_acceptance.
	}
	
	if(verbose)
	{
		cout << "pdf_acceptance  = " << int(pdf_acceptance) << endl;
		cout << "holstein_acceptance = " << int(holstein_acceptance) << endl;
		cout << "jtw_acceptance = " << int(jtw_acceptance) << endl;
		cout << "det_acceptance  = " << int(det_acceptance) << endl;
		cout << "Event accepted? = " << int(det_acceptance && (holstein_acceptance || jtw_acceptance)) << endl;
	}
	// permissive:
	return ( det_acceptance && (holstein_acceptance || jtw_acceptance) );
}


void HolsteinDecay::randomize_nuclear(bool doit)
{
	bool verbose = false;
	//use_roc ?
	
	if(verbose)
	{
		if(doit) { cout << "Randomizing nuclear parameters!" << endl; }
		else     { cout << "Not randomizing nuclear parameters." << endl; }
	}
	
	// shoot() vs. fire():
	//   "shoot" *uses* 'the static generator', while "fire" bypasses 
	//   the static generator and instead uses the 'localEngine'.
	delta_M1 = 0.0*MeV;
	delta_M2 = 0.0*MeV;
	if(doit)
	{
		delta_M1 = ( G4RandGauss::shoot(0, Params->sigma_M1/MeV) )*MeV;  // MeV
		delta_M2 = ( G4RandGauss::shoot(0, Params->sigma_M2/MeV) )*MeV;  // MeV
	}
	
	// masses and energies.  These are all in MeV.
	G4double delta_M = (delta_M1 + delta_M2)/2.0;  // both things are already in G4units of MeV.
	G4double delta_Delta = (delta_M1 - delta_M2);  // both things are already in G4units of MeV.
	
	M     = Params->M     + delta_M;      // both things in G4units of MeV already.
	Delta = Params->Delta + delta_Delta;  // both things in G4units of MeV already.
	m_e   = Params->m_e;     // don't even bother with this uncertainty right now.  comes out in G4units from HolsteinVars.
	E0    = (Delta/MeV)*(1.0+(m_e/MeV)*(m_e/MeV)/(2.0*(M/MeV)*(Delta/MeV)) )/(1.0+(Delta/MeV)/(2.0*(M/MeV))) * MeV;  // error is propagated already;  M1, M2 --> M, Delta --> E0.
	rc    = 1.0 + (Delta/MeV)/(2.0*(M/MeV));  // unitless.
	
//	// delta_E0 isn't really needed.  It's more vestigial.
//	double delta_E0 = (0.5 + ( pow(M-0.5*Delta, 2)-m_e*m_e )/( 2.0*pow(M+0.5*Delta,2) ))*delta_M1 + (-M+0.5*Delta)/(M+0.5*Delta)*delta_M2;
	
	// isospin correction (goes into M_F and M_GT):
	double delta_deltaC = 0;
	if(doit)
	{
		delta_deltaC = G4RandGauss::shoot(0, Params->sigma_deltaC);
	}
	deltaC = Params->deltaC + delta_deltaC;
	M_F  = Params->M_F  * sqrt(1.0-Params->deltaC) / sqrt(1.0-deltaC);  // uncorrect the prev. correction, they re-correct with error.
	M_GT = Params->M_GT * sqrt(1.0-Params->deltaC) / sqrt(1.0-deltaC);  // this probably doesn't capture all the physics of M_GT and deltaC ...
	/// !!!  MELISSA LOOK HERE!  
	/// !!!  I THINK THESE VALUES FOR M_F AND M_GT HAVE THE WRONG ISOSPIN SYMMETRY ERROR PROPAGATION, 
	/// !!!  AND NOT JUST BECAUSE OF SUBTLETIES!
	
	
	// magnetic moments:
	if(doit)
	{
		delta_mu_parent   = G4RandGauss::shoot(0, Params->sigma_mu_parent);
		delta_mu_daughter = G4RandGauss::shoot(0, Params->sigma_mu_daughter);
	}
	mu_parent         = Params->mu_parent   + delta_mu_parent;   
	mu_daughter       = Params->mu_daughter + delta_mu_daughter;
	
	// quadrupole moments:
	if(doit)
	{
		delta_quad_parent   = G4RandGauss::shoot(0, Params->sigma_quad_parent);    // MeV
		delta_quad_daughter = G4RandGauss::shoot(0, Params->sigma_quad_daughter);  // MeV
	}
	quad_parent   = Params->quad_parent   + delta_quad_parent;
	quad_daughter = Params->quad_daughter + delta_quad_daughter;
	
	// Use experimental quadrupole moments for 'g'-->M_Q, rather than calculated M_Q.
	double jterm = sqrt( (I_spin + 1.0)*(2.0*I_spin + 3.0)/( I_spin*(2.0*I_spin - 1.0) ));
	M_Q = (2.0/3.0) * M_F_isospin * jterm* (quad_daughter-quad_parent) / (hbarc_eV_nm*hbarc_eV_nm) / ( -4.0/3.0 * g_V );
	// Set M_rp by M_Q:
	M_rp = (E0/MeV)*M_Q/sqrt(6);
	
	// Calculate everything else:
	a1  = g_V*(M_F - (Delta/MeV)*(Delta/MeV)/6.0*M_r2 + (Delta/MeV)/3.0*M_rdotp) / rc;
	a2  = g_V*((M/MeV)*(M/MeV))*(M_r2/6.0) / rc;  
	c1  = g_A*(M_GT - (Delta/MeV)*(Delta/MeV)/6.0*M_sr2 + M_1y*2.0*(Delta/MeV)*(Delta/MeV)/(6.0*sqrt(10.)) + A_nucleons*(Delta/MeV)/(2.0*(M/MeV))*M_sL + (Delta/MeV)/2.0*M_srp ) / rc;
	c2  = g_A*((M/MeV)*(M/MeV))*( M_sr2/6.0 + M_1y/(6.0*sqrt(10.)) ) / rc;
	
	b   = A_nucleons*M_F_isospin * sqrt((I_spin+1.0)/I_spin) * (mu_parent-mu_daughter);
	
	double d_I  = g_A*(-1.0*M_GT + (Delta/MeV)*(Delta/MeV)*M_sr2/6.0 + M_1y*(Delta/MeV)*((M/MeV)+(Delta/MeV)/6.0)/sqrt(10.) + (A_nucleons*M_sL) + ((M/MeV)*M_srp) ) / rc;
	double d_II = -1.0*g_II*A_nucleons*M_GT;
	d   = d_I + d_II;
	
	e   = 0; // Eq. (22).  assumes CVC.
	g   = (2.0/3.0) * (M/MeV)*(M/MeV) * M_F_isospin * jterm*(quad_daughter-quad_parent) / (hbarc_eV_nm*hbarc_eV_nm);
	f   = 2.0*g_V*(M/MeV)*M_rp; // for f, use new value of M_rp
	h   = -1.0*(g_A*2.0*((M/MeV)*(M/MeV))*M_1y/sqrt(10.) + g_P*(A_nucleons*A_nucleons)*M_GT) / rc;
	
	j2  = -g_A*2.0*((M/MeV)*(M/MeV))*M_2y/3.0;
	j3  = -g_A*2.0*((M/MeV)*(M/MeV))*M_3y/3.0;
	
	// must still update the jtw parameters, because I might have changed g_A or M_GT.
	jtw_a1 = g_V*M_F *sqrt(1.0-deltaC);  // turn isospin correction back off for this.  it's already included in M_F.
	jtw_c1 = g_A*M_GT*sqrt(1.0-deltaC);  // turn isospin correction back off for this.  it's already included in M_GT.
	
	jtw_xi    = 2.0*(jtw_a1*jtw_a1 + jtw_c1*jtw_c1);
	jtw_Abeta = 2.0*(jtw_c1*jtw_c1*(1.0/(3.0/2.0 + 1.0)) + 2.0*jtw_a1*jtw_c1*sqrt((3.0/2.0)/(3.0/2.0 + 1.0)) ) / (2.0*(jtw_a1*jtw_a1 + jtw_c1*jtw_c1));
	jtw_rho   = jtw_c1/jtw_a1;
	
	// ok, but actually turn a bunch of this shit off...
	if(!use_roc)
	{
	//	Delta = 0.0;  // it *shouldn't* be zero !  ... I think.
	//	// blah.  don't make any of them zero.  it's not needed.
	//	m_e   = 0.0;
	///	rc    = 1.0;
	//	M     = M*1.0e30;
		
		// could just directly set terms to zero?  shouldn't be needed... I hope.  
		// affects on a2, c2, g, d_I, h, j2, j3 
		// Yeah, let's do it.
		
		// turn isospin mixing off too?  It shouldn't have an effect 
		// on the final result, because all the terms that aren't 
		// corrected go away, and then shit just scales like rho.
		//
		// M_F         = FindValue("M_F")  / sqrt(1.0-deltaC); 
		// M_GT        = FindValue("M_GT") / sqrt(1.0-deltaC);
		// Load them right from the uncorrected text file...
		
	//	M_F  = Params->FindValue("M_F");
	//	M_GT = Params->FindValue("M_GT");
		
	//	jtw_a1 = g_V*M_F *sqrt(1.0-deltaC);  // turn isospin correction back off for this.  it's already included in M_F.
	//	jtw_c1 = g_A*M_GT*sqrt(1.0-deltaC);  // turn isospin correction back off for this.  it's already included in M_GT.
		
		a1 = jtw_a1;
		a2 = 0.0;
		c1 = jtw_c1;
		c2 = 0.0;
		
		b = 0.0;
		d = 0.0;
		e = 0.0;
		g = 0.0;
		
		f = 0.0;
		h = 0.0;
		j2= 0.0;
		j3= 0.0;
	}
	
	if(verbose)
	{
		cout << std::setprecision(12);
		cout << "use_roc = " << use_roc << endl;
		cout << endl;
		cout << "M1     : " << endl;
		cout << "  sigma/MeV: " << Params->sigma_M1/MeV << endl;
		cout << "  delta/MeV: " << delta_M1/MeV << endl;
		cout << "M2     : " << endl;
		cout << "  sigma/MeV: " << Params->sigma_M2/MeV << endl;
		cout << "  delta/MeV: " << delta_M2/MeV << endl;
		cout << endl;
		cout << "-> M/MeV      : " << M/MeV << endl;
		cout << "->    best: " << Params->M/MeV << endl;
		cout << "->   delta: " << (M - Params->M)/MeV << endl;
		cout << "-> Delta/MeV  : " << Delta/MeV << endl;
		cout << "->    best: " << Params->Delta/MeV << endl;
		cout << "->   delta: " << (Delta - Params->Delta)/MeV << endl;
		cout << "-> E0/MeV     : " << E0/MeV << endl;
		cout << "->    best: " << Params->E0/MeV << endl;
		cout << "->   delta: " << (E0 - Params->E0)/MeV << endl;
		cout << "-> rc     : " << rc << endl;
		cout << "->    best: " << 1.0 + Params->Delta/(2.0*Params->M/MeV) << endl;
		cout << "->   delta: " << rc - (1.0 + Params->Delta/(2.0*Params->M/MeV) ) << endl;
		cout << endl;
		cout << "deltaC : " << deltaC << endl;
		cout << "   best: " << Params->deltaC << endl;
		cout << "  sigma: " << Params->sigma_deltaC << endl;
		cout << "  delta: " << delta_deltaC << endl;
		cout << endl;
		cout << "mu_parent  : " << mu_parent << endl;
		cout << "       best: " << Params->mu_parent << endl;
		cout << "      sigma: " << Params->sigma_mu_parent << endl;
		cout << "      delta: " << delta_mu_parent << endl;
		cout << "mu_daughter: " << mu_daughter << endl;
		cout << "       best: " << Params->mu_daughter << endl;
		cout << "      sigma: " << Params->sigma_mu_daughter << endl;
		cout << "      delta: " << delta_mu_daughter << endl;
		cout << endl;
		cout << "quad_parent  : " << quad_parent << endl;
		cout << "         best: " << Params->quad_parent << endl;
		cout << "        sigma: " << Params->sigma_quad_parent << endl;
		cout << "        delta: " << delta_quad_parent << endl;
		cout << "quad_daughter: " << quad_daughter << endl;
		cout << "         best: " << Params->quad_daughter << endl;
		cout << "        sigma: " << Params->sigma_quad_daughter << endl;
		cout << "        delta: " << delta_quad_daughter << endl;

		cout << "- - -" << endl;
		cout << "a1  = " << a1 << endl;
		cout << "best: " << Params->a1 << endl;
		cout << "a2 = " << a2 << endl;
		cout << "best:" << Params->a2 << endl;
		cout << "c1 = " << c1 << endl;
		cout << "best:" << Params->c1 << endl;
		cout << "c2 = " << c2 << endl;
		cout << "best:" << Params->c2 << endl;
		
		cout << "- - -" << endl;
		cout << "b  = " << b << endl;
		cout << "best:" << Params->b << endl;
		cout << "d  = " << d << endl;
		cout << "best:" << Params->d << endl;
	//	cout << "e  = " << e << endl;
	//	cout << "best:" << Params->e << endl;
		cout << "g  = " << g << endl;
		cout << "best:" << Params->g << endl;
		cout << "f  = " << f << endl;
		cout << "best:" << Params->f << endl;
		cout << "h  = " << h << endl;
		cout << "best:" << Params->h << endl;
		cout << endl;
		cout << "j2 = " << j2 << endl;
		cout << "best:" << Params->j2 << endl;
		cout << "j3 = " << j3 << endl;
		cout << "best:" << Params->j3 << endl;
		
		cout << endl;
		cout << "Nuclear parameters have been randomized." << endl;
	}
	
	return;
}
void HolsteinDecay::randomize_atomic(bool doit)
{
	bool verbose=false;
	if(verbose)
	{
		if(doit) { cout << "Randomizing atomic parameters! (or I would be if this function did anything yet.)" << endl; }
		else     { cout << "Not randomizing atomic parameters.  It's not implemented yet anyway." << endl; }
	}
	
	// Do this step whether I randomize atomic params or not, because I might have adjusted the polarization before this event.
	Mz  = this->the_pops->get_Mz();
	Mz2 = this->the_pops->get_Mz2();
	Mz3 = this->the_pops->get_Mz3();
	
	initialize_lambdafuncs();  // this depends on the sublevel populations, which should be set up above.
	return;
}

void HolsteinDecay::initialize_lambdafuncs()
{
	get_lambda0();  // 1.0
	get_lambda1();
	get_lambda2();
	get_lambda3();
	get_lambda4();  // 0.0
}
double HolsteinDecay::get_lambda0()
{
	Lambda0=1; 
	return Lambda0;
}
double HolsteinDecay::get_lambda1()
{
	Lambda1 = Mz/I_spin; // Polarization.
	return Lambda1;
}
double HolsteinDecay::get_lambda2()
{
	double T = 5.0/4.0 - Mz2;
	Lambda2 = T * (2.0*I_spin - 1.0)/(I_spin+1.0);
	return Lambda2;
}
double HolsteinDecay::get_lambda3()
{
	double numerator =   5.0*Mz3;
	double denominator = I_spin*(3.0*I_spin*I_spin + 3.0*I_spin - 1.0);
	Lambda3 = Mz/I_spin - numerator/denominator;

	return Lambda3;
}
double HolsteinDecay::get_lambda4()
{
	Lambda4=0; 
	return Lambda4; 
}

double HolsteinDecay::F_0(double E) // a1 a2 c1 c2 (M/MeV) (E0/MeV) (m_e/MeV) b d h
{
	bool verbose=false;
	double t1, t2, t3, t4, t5;
	t1 = a1*a1 + 2.0*a1*a2/(3.0*(M/MeV)*(M/MeV))*(     (m_e/MeV)*(m_e/MeV) +  4.0*E*(E0/MeV) + 2.0*(m_e/MeV)*(m_e/MeV)*(E0/MeV)/E -  4.0*E*E);
	t2 = c1*c1 + 2.0*c1*c2/(9.0*(M/MeV)*(M/MeV))*(11.0*(m_e/MeV)*(m_e/MeV) + 20.0*E*(E0/MeV) - 2.0*(m_e/MeV)*(m_e/MeV)*(E0/MeV)/E - 20.0*E*E);
	t3 = -2.0*(E0/MeV)/(3.0*(M/MeV))*(c1*c1 + c1*d - c1*b);
	t4 =  2.0*E /(3.0*(M/MeV))*(3.0*a1*a1 + 5.0*c1*c1 -2.0*c1*b);
	t5 = -1.0*(m_e/MeV)*(m_e/MeV)/(3.0*(M/MeV)*E)*(-3.0*a1*e + 2.0*c1*c1 + c1*d - 2.0*c1*b - c1*h*((E0/MeV)-E)/(2.0*(M/MeV)));
	
	double F0 = t1 + t2 + t3 + t4 + t5;
	if(verbose) 
	{ 
		cout << "F0 = " << F0 << endl;
		cout << "a1*a1 + c1*c1 = " << a1*a1 + c1*c1 << endl;   // JTW xi*factor.
	}
	return F0;
}
double HolsteinDecay::F_1(double E) // a1 c1 c2 a2 (E0/MeV) (m_e/MeV) (M/MeV) b d j2 g f 
{
	bool verbose=false;
	
	double f1, f2, f3;
	double t1a, t1b, t2a, t2b, t3a, t3b;
	f1 = Params->delta_uv*sqrt(Params->u/(Params->u+1)) * 2.0;
	t1a = a1*c1 - ((E0/MeV)/(3.0*(M/MeV)))*(a1*c1 + a1*d - a1*b) + (E/(3.0*(M/MeV)))*(7.0*a1*c1 - a1*b + a1*d);  
	t1b = (a1*c2 + c1*a2)*(4.0*E*((E0/MeV)-E)+3.0*(m_e/MeV)*(m_e/MeV))/(3.0*(M/MeV)*(M/MeV));
	
	f2 = Params->gamma_uv/(Params->u+1.0);
	t2a = c1*c1 + 2.0*c1*c2*(8.0*E*((E0/MeV)-E) + 3.0*(m_e/MeV)*(m_e/MeV))/(3.0*(M/MeV)*(M/MeV));
	t2b = -2.0*(E0/MeV)/(3.0*(M/MeV))*(c1*c1 + c1*d - c1*b) + E/(3.0*(M/MeV))*(11.0*c1*c1 - c1*d - 5.0*c1*b);
	
	f3 = Params->lambda_uv/(Params->u+1.0);
	t3a = -5.0*E/(M/MeV)*c1*f + sqrt(3.0/2.0)*c1*g*((E0/MeV)*(E0/MeV) - 11.0*(E0/MeV)*E + 6.0*(m_e/MeV)*(m_e/MeV) + 4.0*E*E)/(6.0*(M/MeV)*(M/MeV));
	t3b = -3.0*c1*j2*(8.0*E*E - 5.0*(E0/MeV)*E - 3.0*(m_e/MeV)*(m_e/MeV))/(6.0*(M/MeV)*(M/MeV));
	
	double F1 = f1*(t1a + t1b) + f2*(t2a + t2b) + f3*(t3a + t3b);
	if(verbose) 
	{ 
		cout << "F1 = " << F1 << endl;
		cout << "c1*c1*(1.0/(3.0/2.0 + 1.0)) + 2*a1*c1*sqrt( (3.0/2.0) /(3.0/2.0 + 1.0)) = " << c1*c1*(1.0/(3.0/2.0 + 1.0)) + 2.0*a1*c1*sqrt( (3.0/2.0) /(3.0/2.0 + 1.0)) << endl; // JTW xi*factor*A_beta.
	} 
	return F1;
}
double HolsteinDecay::F_2(double E) // c1 (E0/MeV) (M/MeV) a1 d b g j2 f j3 c2 (no dep. on a2 ??)
{
	bool verbose=false;
	
	// did I double-check this?
	double f2;
	double t1, t2, t3, t4;
	t1 = Params->theta_uv*E/(2.0*(M/MeV))*(c1*c1 + 8.0*c1*c2*((E0/MeV)-E)/(3.0*(M/MeV)) - c1*d - c1*b);
	f2 = -1.0*Params->delta_uv*E/(M/MeV)*sqrt( (Params->u*(Params->u+1.0))/((2.0*Params->u-1.0)*(2.0*Params->u+3.0)) );
	t2 = sqrt(3.0/2.0)*a1*f + a1*g*(E+2.0*(E0/MeV))/(4.0*(M/MeV)) - sqrt(3.0/2.0)*a1*j2*((E0/MeV)-E)/(2.0*(M/MeV));
	t3 = Params->kappa_uv*E/(2.0*(M/MeV))*(-3.0*c1*f -sqrt(3.0/2.0)*c1*g*((E0/MeV)-E)/(M/MeV) + 3.0*c1*j2*((E0/MeV)-2.0*E)/(2.0*(M/MeV)));
	t4 = Params->epsilon_uv*c1*j3*(21.0*E*E)/(8.0*(M/MeV)*(M/MeV));
	
	double F2 = t1 + f2*t2 + t3 + t4;
	if(verbose) { cout << "F2 = " << F2 << endl; }
	return F2;
}
double HolsteinDecay::F_3(double E) // a1 j3 (M/MeV) c1 g j2 j3
{
	bool verbose=false;
	
	// did I double-check this?
	double f1;
	double t1, t2, t3;
	
	f1 = -1.0*Params->delta_uv*(3.0*Params->u*Params->u + 3.0*Params->u-1.0)*sqrt(Params->u)/sqrt( (Params->u-1.0)*(Params->u+1.0)*(Params->u+2.0)*(2.0*Params->u-1.0)*(2.0*Params->u+3.0) );
	t1 = a1*j3*sqrt(15.0)*E*E/(4.0*(M/MeV)*(M/MeV));
	t2 = Params->rho_uv/(Params->u+1.0)*(5.0*E*E)/(4.0*(M/MeV)*(M/MeV))*(sqrt(3)*c1*g + sqrt(2)*c1*j2);
	t3 = -1.0*Params->sigma_uv/(Params->u+1.0)*(5.0*E*E)/(2.0*(M/MeV)*(M/MeV))*c1*j3;
	
	double F3 = f1*t1 + t2 + t3;
	if(verbose) { cout << "F3 = " << F3 << endl; }
	return F3;
}

G4double HolsteinDecay::pbeta(G4double E)  // *only* in get_probability so far...
{
	G4double pbetac = sqrt( pow(E/MeV,2) - pow((m_e/MeV),2) );  // MeV for everything...
	return pbetac;  // comes out in G4units of MeV too.
}
G4double HolsteinDecay::get_Ebeta(G4double pbeta)  // not used???
{
	G4double E = sqrt( pow(pbeta/MeV,2) + pow((m_e/MeV), 2) );
	return E;
}
double HolsteinDecay::get_v_from_p(G4double pbeta) // [v] = unitless m/s, pbeta in G4units of MeV.
{
	double beta = (pbeta/MeV)* sqrt( (1.0)/((m_e/MeV)*(m_e/MeV) + (pbeta/MeV)*(pbeta/MeV)) );
	double v = beta*speed_of_light;  // meters/second 
	return v;
}
G4double HolsteinDecay::get_p_from_v(double vbeta)  // ...
{
	double gamma = 1.0 / sqrt(1.0 - vbeta*vbeta/(speed_of_light*speed_of_light)); // unitless....
	G4double pbeta = gamma*(m_e/MeV)*vbeta/speed_of_light * MeV;
	return pbeta;
}

double HolsteinDecay::get_probability(G4double E, double costheta) // want to use costheta_lab
{
	int verbose = 0;
	
	double Gv2       = 1.0; // kludge
	double cos2theta = 1.0; // kludge
	double prefactor = 2.0*Gv2*cos2theta/( pow(2.0*pi, 4) );
	
	double the_scaling = prefactor*FermiFunction(Params->Z_daughter, (E/MeV))*pow((E0/MeV)-(E/MeV), 2)*(E/MeV)*pbeta(E);
	
	if(verbose>0)
	{
		double the_F0, the_F1;  // just to check on ROC.
		the_F0 = F_0((E/MeV));
		the_F1 = F_1((E/MeV));
		cout << "F_0     = " << the_F0        << ";\tfake_xi/2  = " << (a1*a1 + c1*c1) << endl;
		cout << "F_1/F_0 = " << the_F1/the_F0 << ";\tfake_Abeta = " << 2.0*(c1*c1*(1.0/(3.0/2.0 + 1.0)) + 2.0*a1*c1*sqrt((3.0/2.0)/(3.0/2.0 + 1.0)) ) / (2.0*(a1*a1 + c1*c1)) << endl; // sign weirdness.
		cout << "fake-er xi/2  = " << (jtw_a1*jtw_a1 + jtw_c1*jtw_c1) << ";\t";
		cout << "fake-er Abeta = " << 2.0*(jtw_c1*jtw_c1*(1.0/(3.0/2.0 + 1.0)) + 2.0*jtw_a1*jtw_c1*sqrt((3.0/2.0)/(3.0/2.0 + 1.0)) ) / (2.0*(jtw_a1*jtw_a1 + jtw_c1*jtw_c1)) << endl;
		cout << "a1 = " << a1 << ";\tc1 = " << c1 << endl;
	}
	
	double t0 = F_0((E/MeV));
	double t1 = F_1((E/MeV)) * Lambda1 * costheta*(pbeta(E)/MeV)/(E/MeV);
	double t2 = F_2((E/MeV)) * Lambda2 * (pbeta(E)/MeV)*(pbeta(E)/MeV) / ((E/MeV)*(E/MeV)) * (costheta*costheta - 1.0/3.0);
	double t3 = F_3((E/MeV)) * Lambda3 * ( pow(costheta*(pbeta(E)/MeV)/(E/MeV), 3) - 3.0/5.0*costheta*pow((pbeta(E)/MeV)/(E/MeV), 3) );
	
	double the_prob = the_scaling*(t0+t1+t2+t3);
	holstein_probability = the_prob;
	holstein_Abeta       = F_1((E/MeV))/t0;
	
	if(verbose>1)
	{
		cout << "HolsteinDecay::get_probability:  the_prob = " << holstein_probability << endl;
		cout << "\tthe_scaling = " << the_scaling << endl;
		cout << "\tholstein_Abeta = " << holstein_Abeta << endl;
		cout << "\tt0 = " << t0 << endl;
		cout << "\tt1 = " << t1 << endl;
		cout << "\tt2 = " << t2 << endl;
		cout << "\tt3 = " << t3 << endl;
	}
	return holstein_probability;
}


double HolsteinDecay::get_jtw_probability(G4double E, double costheta)
{
	double Gv2       = 1.0; // kludge
	double cos2theta = 1.0; // kludge
	double prefactor = 2.0*Gv2*cos2theta/( pow(2.0*pi, 4) );
	double the_scaling = prefactor*FermiFunction(Params->Z_daughter, (E/MeV))*pow((E0/MeV)-(E/MeV), 2)*(E/MeV)*pbeta(E);
	
	double fake_F0 = (jtw_a1*jtw_a1 + jtw_c1*jtw_c1);
	double fake_F1 = jtw_c1*jtw_c1*(1.0/(3.0/2.0 + 1.0)) + 2.0*jtw_a1*jtw_c1*sqrt((3.0/2.0)/(3.0/2.0 + 1.0));
	
	double the_term = fake_F0 + fake_F1 * Lambda1 * costheta*(pbeta(E)/MeV)/(E/MeV);
	double the_prob = the_scaling*the_term;
	
	jtw_probability = the_prob;
	return jtw_probability;
}


void HolsteinDecay::randomize_direction() // saves initial_momentum and initial_velocity in class.  Also Ebeta.  
{
	bool verbose=false;
	
	// pick E
//	double E = G4RandFlat::shoot(m_e, (E0/MeV));
//	Ebeta = E;    // set Ebeta *here*.
	Ebeta = G4RandFlat::shoot(m_e/MeV, E0/MeV)*MeV;  // m_e and E0 have to be in the same units.  also, in units of MeV.
	
	// pick costheta
	costheta_lab = G4RandFlat::shoot(1.0, -1.0);
	while( use_cone && (abs(costheta_lab) <= cone_costheta) ) // the max. is on theta.  cos(theta_max) is the min.
	{
		costheta_lab = G4RandFlat::shoot(1.0, -1.0);
	}
	
	// pick phi
	double phi = G4RandFlat::shoot(0.0, 2.0*pi);
	
	double theta_lab = acos(costheta_lab);
	double the_velocity = get_v_from_p( pbeta(Ebeta) );  // comes out as just a regular double.
	
	initial_velocity.setRThetaPhi(the_velocity, theta_lab, phi);
	
	// Don't set momentum in cartesian coordinates.  spherical coordinates are like a rotation to the best axis.
	// initial_momentum.setRThetaPhi( get_p_from_v(the_velocity), theta_lab, phi); 
	initial_momentum.setRThetaPhi( pbeta(Ebeta), theta_lab, phi); // initial_momentum comes out in G4units of MeV, because those are the units that pbeta has.
	
	// for export...
	Ebeta_tot_MeV = Ebeta/MeV;
	pbeta_MeV     = initial_momentum.mag()/MeV;
	Ebeta_kin_MeV = (Ebeta - m_e)/MeV;
	vbeta_over_c  = initial_velocity.mag();

	if(verbose)
	{
		// test:
		cout << "*" << endl;
		cout << "(Ebeta/MeV)                                                     = " << (Ebeta/MeV) << endl;
		cout << "pbeta(Ebeta)                                              = " << pbeta(Ebeta) << endl;
		cout << "get_v_from_p( pbeta(Ebeta) )                              = " << get_v_from_p( pbeta(Ebeta) ) << endl;
		cout << "    initial_velocity.mag()                                = " << initial_velocity.mag() << endl;
		cout << "    vx/c=" << initial_velocity.x()/speed_of_light << ";\tvy/c=" << initial_velocity.y()/speed_of_light << ";\tvz/c=" << initial_velocity.z()/speed_of_light << endl;
		cout << "    sqrt(vx^2 + vy^2 + vz^2)                              = ";
		cout << sqrt( initial_velocity.x()*initial_velocity.x() 
					+ initial_velocity.y()*initial_velocity.y() 
					+ initial_velocity.z()*initial_velocity.z() ) << endl;

		cout << "get_p_from_v( get_v_from_p( pbeta(Ebeta) ) )/MeV              = " << get_p_from_v( get_v_from_p( pbeta(Ebeta)) )/MeV << endl;
		cout << "    initial_momentum.mag()/MeV                                = " << initial_momentum.mag()/MeV << endl;
		cout << "    px/MeV=" << initial_momentum.x()/MeV << ";\tpy/MeV=" << initial_momentum.y()/MeV << ";\tpz/MeV=" << initial_momentum.z()/MeV  << endl;
		cout << "    sqrt(px^2 + py^2 + pz^2)                              = ";
		cout << sqrt( initial_momentum.x()/MeV*initial_momentum.x()/MeV 
					+ initial_momentum.y()/MeV*initial_momentum.y()/MeV 
					+ initial_momentum.z()/MeV*initial_momentum.z()/MeV ) << endl;
		cout << "get_Ebeta( get_p_from_v( get_v_from_p( pbeta(Ebeta) ) ) )/MeV = " << get_Ebeta( get_p_from_v( get_v_from_p( pbeta(Ebeta)) ) )/MeV << endl;
		
	}
	
	return;
}


bool HolsteinDecay::check_PDF_acceptance()  // uses initial_momentum .
{
	bool verbose = false;
	
	pdf_acceptance      = false;
	jtw_acceptance      = false;
	holstein_acceptance = false;
	
	if(initial_momentum.mag() == 0)
	{
		cout << "No, this probably hasn't been initialized.  Rejected!"  << endl;
		assert(0);
	//	pdf_acceptance = false;
	//	return pdf_acceptance;
	}
//	double prob_max = 0.2;  
	double test_prob = G4RandFlat::shoot(0.0, prob_max);
	
	holstein_probability = get_probability(Ebeta, costheta_lab ); // this is redundant, but whatevs.
	jtw_probability = get_jtw_probability(Ebeta, costheta_lab);
	the_probability = holstein_probability;
	
	if(verbose)
	{
		cout << "the_prob = " << the_probability << ";  holstein_prob = " << holstein_probability << ";  jtw_prob = " << jtw_probability << ";  max = " << prob_max << ";  test_prob = " << test_prob << endl;
	}
	//
	if(holstein_probability > prob_max)
	{
		cout << "Bad!  We've broken our acceptance/rejection method!" << endl;
		cout << "You should increase prob_max within the code.  " << endl;
		cout << "prob_max = " << prob_max << endl;
		cout << "prob     = " << the_probability << endl;
		assert(0); // hard kill.
		return true;
	}
	if(jtw_probability > prob_max)
	{
		cout << "You broke JTW." << endl;
		// don't kill.  
	}
	
	
	// Do we accept this event?  Check if the random number falls *under* the probability distribution.
	// return the Holstein acceptance, even though we calculated the jtw acceptance above.
	
	// Check JTW Acceptance:
	if( test_prob <= jtw_probability)
	{
		jtw_acceptance = true;
		if(verbose) { cout << "\tEvent accepted by JTW PDF!" << endl; }
	}
	else
	{
		jtw_acceptance = false;
		if(verbose) { cout << "\tEvent rejected by JTW PDF." << endl; }
	}
	// Check Holstein Acceptance:
	if( test_prob <= holstein_probability)
	{
		holstein_acceptance = true;
		if(verbose) { cout << "\tEvent accepted by Holstein PDF!" << endl; }
	}
	else
	{
		holstein_acceptance = false;
		if(verbose) { cout << "\tEvent rejected by Holstein PDF." << endl; }
	}
	
	// permissive.
	pdf_acceptance = (holstein_acceptance || jtw_acceptance);
//	cout << "holstein_acceptance=" << holstein_acceptance << endl;
//	cout << "pdf_acceptance=" << pdf_acceptance << endl;
	return pdf_acceptance;
}

void HolsteinDecay::randomize_start(bool doit)
{
	bool verbose=false;
	// pick start position.
	if(verbose)
	{
		if(doit) { cout << "Trying to randomize start position, but this feature isn't actually implemented." << endl; }
		else     { cout << "Selecting non-randomized start position.  It doesn't randomize yet anyway." << endl; }
	}
	initial_position.set(0.0*mm,0.0*mm,0.0*mm);
	if(doit)
	{
		// pick some other position.
	}
	return;
}

bool HolsteinDecay::check_detector_acceptance()
{
	bool verbose=false;
	// Should have already randomized start, randomized direction, and checked PDF acceptance.
	det_acceptance = false;
	// propagate through the vacuum...
	
	// going up, or down?
	if(initial_velocity.z() > 0)
		{ going_up = true; }
	else if( initial_velocity.z() < 0)
		{ going_up = false; }
	else if( initial_velocity.z() == 0 )
	{ 
		cout << endl;
		cout << "* No vertical velocity!  It's a mathematical miracle!" << endl;
		cout << endl;
		det_acceptance = false;
		return det_acceptance;  // I mean, it's not going to hit...
	}
	
	// How far does it have to go before it hits?
	double vertical_distance_to_travel;
	if(going_up)
	{
		vertical_distance_to_travel = the_geometry.vdistance_center_to_dssd/mm - initial_position.z()/mm;
	}
	if(!going_up)
	{
		vertical_distance_to_travel = initial_position.z()/mm - the_geometry.vdistance_center_to_dssd/mm;  // negative.
	}
	
//	double time_to_travel = vertical_distance_to_travel / initial_velocity.z(); 
	time_to_travel = vertical_distance_to_travel / initial_velocity.z(); 
	// no negative times.
	if(time_to_travel<=0) { cout << "Bad.  Negative times aren't allowed." << endl; assert(0); }
	
	hit_position.set(
		(initial_position.x()/mm + time_to_travel*initial_velocity.x() )*mm,
		(initial_position.y()/mm + time_to_travel*initial_velocity.y() )*mm,
		(initial_position.z()/mm + time_to_travel*initial_velocity.z() )*mm );  // these are all the wrong units!!!
	
	if(verbose)
	{
		cout << "\'Hit\' position:  " << endl;
		cout << "\tx/mm = " << hit_position.x()/mm << endl;
		cout << "\ty/mm = " << hit_position.y()/mm << endl;
		cout << "\tz/mm = " << hit_position.z()/mm << endl;
	}
	
	if( abs(hit_position.x()/mm) <= (the_geometry.dssd_width_x/mm)/2.0 && abs(hit_position.y()/mm) <= (the_geometry.dssd_width_y/mm)/2.0 )
	{
		det_acceptance = true;
	}
	return det_acceptance;
}

void HolsteinDecay::print_results()
{
	cout << "holstein_acceptance: " << holstein_acceptance << endl; 
	cout << "pdf_acceptance: " << pdf_acceptance << endl; 
	cout << "jtw_acceptance: " << jtw_acceptance << endl; 
	cout << "det_acceptance: " << det_acceptance << endl; 
	
	//
	cout << "Detector:  ";
	if(going_up) { cout << "Top."    << endl; }
	else         { cout << "Bottom." << endl; }
	
	cout << "(Ebeta/MeV) = " << (Ebeta/MeV) << endl;
	cout << "pbeta/MeV = " << initial_momentum.mag()/MeV << endl;
	cout << "costheta = " << costheta_lab << endl;
	cout << "time_to_travel = " << time_to_travel << " (in some units)" << endl;
	
	cout << "Initial momentum:  " << endl;
	cout << "    px/MeV = " << initial_momentum.x()/MeV << ";\tvx/c = " << initial_velocity.x()/speed_of_light << endl;
	cout << "    py/MeV = " << initial_momentum.y()/MeV << ";\tvy/c = " << initial_velocity.y()/speed_of_light << endl;
	cout << "    pz/MeV = " << initial_momentum.z()/MeV << ";\tvz/c = " << initial_velocity.z()/speed_of_light << endl;

	cout << "'Hit' position:  " << endl;
	cout << "    x/mm  = " << hit_position.x()/mm << endl;
	cout << "    y/mm  = " << hit_position.y()/mm << endl;
	cout << "    z/mm  = " << hit_position.z()/mm << endl;
	cout << "--" << endl;
}

// // // // 
int main(int argc, char *argv[]) 
{
	cout << "Hello, world!" << endl;
	timeval t1;
	gettimeofday(&t1, NULL);
	G4long randseed = (t1.tv_usec * t1.tv_sec);// + pid;
	CLHEP::HepRandom::setTheSeed(randseed);
	G4cout << "randseed: " << randseed << G4endl;
	
	
	string filename = "output.root";
	
	sublevel_populations * thepops     = new sublevel_populations(1);
	HolsteinVars         * pointervars = new HolsteinVars();
	HolsteinDecay        * the_decay   = new HolsteinDecay(pointervars, thepops);
	
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
	
	
	return 0;
}











