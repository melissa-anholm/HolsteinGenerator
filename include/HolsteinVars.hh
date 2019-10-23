#ifndef HolsteinVars_h
#define HolsteinVars_h 1

#include <iostream>  // cout, endl
#include <string> 
#include <cmath>     // pow, atan
#include <map>
//#include "fstream"   // 

#include <G4Types.hh>
#include <G4SystemOfUnits.hh>
#include "G4PhysicalConstants.hh"


#include "IsotopeValues.hh"
//#include "SplitString.hh"

using std::map;

//const double pi = std::atan(1.0)*4.0;

// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- //
struct chamber_geometry // numbers indirectly imported from G4.
{
	// distance to the *front* of the dssd is 102.127357.
	G4double vdistance_center_to_dssd = 103.627357*mm; // this is probably a bit too precise.  wev.  it's in mm.
	G4double dssd_width_x = 40.0*mm; // mm.
	G4double dssd_width_y = 40.0*mm; // mm.
	G4double dssd_cut_radius = 15.5*mm; // mm.  probably it's only good for post-processing...  unused.
};


// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- //


//  class HolsteinVars relies on:
//  	class isotope_values
//  	class SplitString
class HolsteinVars
{
public:
	HolsteinVars();
		
	bool initialize_spinfuncs(double u_, double v_);
	double delta_uv, gamma_uv, lambda_uv, theta_uv, kappa_uv, epsilon_uv, rho_uv, sigma_uv, phi_uv;
	
	// Loading the text file.  Implementing code I inherited from Spencer.
	// string paramfilename;
	string nuclear_filename;
//	string atomic_filename;
	
	string isotopeName;      // gets a value loaded into it, but probably doesn't ever get used.
	map<string, isotope_values * > theInputs;
//	map<string, isotope_values * > loadup_textfile(string paramfilename);
//	bool loadup_textfile(string paramfilename);
	
	void print_isotope_values();  // for debugging.  prints from 'theInputs'.
	double FindValue(const std::string &key_) const;
	double FindUncertainty(const std::string &key_) const;
	
	void initialize_physics_parameters();
	bool is_twopercent;
	
	void print_vars();
	void print_matrixelements();
	void print_couplingconstants();
	void print_calculatedJTW();
	void print_holsteinalphabet();
	
	//
	double I_spin, u, v;
	double Z_parent, N_parent, A_nucleons;
	double Z_daughter, N_daughter;
	double T_isospin, T3_parent, T3_daughter, M_F_isospin;

	double hbarc_eV_nm;  // don't use baked in G4 "hbarc" value, because we want things multiplied by this to be unitless.  At least some of the time...
	
	double amu_to_mev, sigma_amu_to_mev;  // this is literally to convert units already.  it'll probably become obsolete.
	double speed_of_light;
	
	// mass, energy, all in MeV.
	G4double E0, M, Delta;  // HolsteinVars: G4units propagated for Delta, E0, M.
	G4double sigma_M1, sigma_M2;  // HolsteinVars:  G4 units propagated for both.  in MeV.
	G4double m_e;//, sigma_m;  // HolsteinVars:  G4units propagated for m_e.  and sigma_m, even though it's unused.
	
	double mu_parent, sigma_mu_parent, mu_daughter, sigma_mu_daughter;
	double quad_parent, sigma_quad_parent, quad_daughter, sigma_quad_daughter;
	double deltaC, sigma_deltaC;
	
	// coupling constants (Ian Towner).
	double g_V, g_A;
	double g_II, g_S, g_P, g_M;
	// lower-case decay parameters:
	double a1, a2, c1, c2;
	double b, d, e, f, g, h;
	double j2, j3;
	
	// matrix elements (Ian Towner).
	double M_F, M_GT;
	double M_r2, M_sr2, M_1y, M_sL, M_srp, M_rdotp, M_rp, M_2y, M_3y;
	double M_L, M_Q; 
	
//	// JTW parameters, with no isospin corrections.
//	double jtw_a1, jtw_c1;
};


#endif
