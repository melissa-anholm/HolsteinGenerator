// Author: Melissa Anholm - 2019

#ifndef Holstein52Generator_h
#define Holstein52Generator_h 1

#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh> // probably the correct way to include ThreeVector.h.


#include "Holstein52Isotope.hh"  // formerly HolsteinVars
// #include "K37SublevelPopulations.hh"
#include "K37AtomicSetup.hh"


// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- //
//  Holstein52Generator class relies on other classes: 
//  	HolsteinVars
//  		isotope_values
//  		SS (splitstring)
//  	K37SublevelPopulations

//class HolsteinDecay
class Holstein52Generator
{
public:
//	Holstein52Generator();    /// this one seems broken?!?!?
	Holstein52Generator(HolsteinVars * HV, K37AtomicSetup * atomic_setup);
//	Holstein52Generator(HolsteinVars  HV) { Holstein52Generator( (HolsteinVars*)(&HV) ); };
	
	// The cone!
	void set_use_cone(bool useit) { use_cone = useit; };
	bool get_use_cone() { return use_cone; };
	bool use_roc;  // don't even bother protecting it yet...
	
	// probably should put some sanity checks into functions that set cone_costheta...
	void set_conecostheta(double costheta_max_)  { cone_costheta = costheta_max_; }; 
	double get_conecostheta() { return cone_costheta; };
	
	void randomize_direction();
	void randomize_direction_monoenergetic(G4double the_monoenergy);

	bool check_PDF_acceptance();         // uses initial_momentum
	bool check_detector_acceptance();    // creates hit_position from initial_momentum and initial_position.  then checks.
	void randomize_nuclear(bool doit=true);
//	void randomize_atomic(bool doit=true);
//	void randomize_start(bool doit=true);   // initial_position, in G4 mm.

	bool shoot_decayevent();
	bool shoot_monoenergetic_decayevent(G4double the_monoenergy);
	
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
	
private:
	G4ThreeVector initial_momentum;  // needs units of energy.
public:
	G4double get_beta_Px() { return initial_momentum.x(); };  // must return something in units of energy.  does it actually?  who the fuck knows.  ... I think it's ok.
	G4double get_beta_Py() { return initial_momentum.y(); };  // must return something in units of energy.
	G4double get_beta_Pz() { return initial_momentum.z(); };  // must return something in units of energy.
	
	// functions like from K37EventGenerator:  
	G4double get_electron_T_MeV()      { return Ebeta_kin_MeV; }
	G4double get_electron_E_MeV()      { return Ebeta_tot_MeV; }  // ok, but what does the rest of the G4 code think this thing is in units of?
	double get_electron_costheta()     { return costheta_lab; }
	
public:
	G4ThreeVector initial_velocity;  // needs units.  not unitless.
	G4ThreeVector initial_position;  // needs units of position.  ...done in G4units.
	G4ThreeVector hit_position;      // needs units of position.
	
	bool pdf_acceptance; // obsolete.
	bool jtw_acceptance, holstein_acceptance;
	bool det_acceptance;
	bool going_up;
	void SetAcceptanceMode(string);
	
	K37AtomicSetup * GetAtomicSetup() { return the_atomic_setup; };
	
	G4double pbeta(G4double E);  
	G4double get_Ebeta(G4double pbeta);
	double get_v_from_p(G4double pbeta);  
	G4double get_p_from_v(double vbeta);
	
	double FermiFunction(double Z_, G4double E) { return 1.0; };
	
	
	chamber_geometry the_geometry;

	double Mz, Mz2, Mz3;
//	void run_fast(bool doit) { runfast=doit; };
	
	// For export only:
	double jtw_xi, jtw_Abeta, jtw_rho; // not a function of E.  const. for a given set of parameters.
	double holstein_Abeta; // event-specific.  function of E.

private:
	// run parameters..
	bool use_cone;
	double cone_costheta;
//	bool runfast; // don't look at detector acceptance if it's not accepted to the PDF.  Should check detector acceptance first, maybe...
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
	//
	K37AtomicSetup * the_atomic_setup;
//	K37Cloud * the_cloud;
//	K37SublevelPopulations * the_pops;  // Do I *really* want this to be public??
};

#endif
