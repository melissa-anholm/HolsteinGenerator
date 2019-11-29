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

using std::map;


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
	
private:
	// to make HolsteinVars non-copyable
	HolsteinVars(const HolsteinVars &);
	// to make HolsteinVars non-copyable
	const HolsteinVars & operator=(const HolsteinVars &);

public:
	bool initialize_spinfuncs(double u_, double v_);
	
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
	bool is_twopercent;  // it's never two percent, within this class.  this set of stuff only gets called in the 98% case.
	
	void print_vars();
	void print_matrixelements();
	void print_couplingconstants();
	void print_calculatedJTW();
	void print_holsteinalphabet();
	
	//
private:
	// spin funcs:
	double delta_uv, gamma_uv, lambda_uv, theta_uv, kappa_uv, epsilon_uv, rho_uv, sigma_uv, phi_uv; 

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
	
public:  // getter methods, so that we can copy the variables over to the generator without breaking things.
	// spin funcs:
	double get_delta_uv()   { return delta_uv; };
	double get_gamma_uv()   { return gamma_uv; };
	double get_lambda_uv()  { return lambda_uv; };
	double get_theta_uv()   { return theta_uv; };
	double get_kappa_uv()   { return kappa_uv; };
	double get_epsilon_uv() { return epsilon_uv; };
	double get_rho_uv()     { return rho_uv; };
	double get_sigma_uv()   { return sigma_uv; };
	double get_phi_uv()     { return phi_uv; };
	
	//
	double get_I_spin() { return I_spin; };
	double get_u()      { return u; };
	double get_v()      { return v; };

	double get_Z_parent()    { return Z_parent; };
	double get_N_parent()    { return N_parent; };
	double get_A_nucleons()  { return A_nucleons; };
	
	double get_Z_daughter()  { return Z_daughter; };
	double get_N_daughter()  { return N_daughter; };
	
	double get_T_isospin()   { return T_isospin; };
	double get_T3_parent()   { return T3_parent; }; 
	double get_T3_daughter() { return T3_daughter; }; 
	double get_M_F_isospin() { return M_F_isospin; };

	double get_hbarc_eV_nm()      { return hbarc_eV_nm; };  // don't use baked in G4 "hbarc" value, because we want things multiplied by this to be unitless.  At least some of the time...
	double get_amu_to_mev()       { return amu_to_mev; };
	double get_sigma_amu_to_mev() { return sigma_amu_to_mev; };  // this is literally to convert units already.  it'll probably become obsolete.
	double get_speed_of_light()   { return speed_of_light; };

	// mass, energy, all in MeV.
	G4double get_E0()                { return E0; };
	G4double get_M()                 { return M; };
	G4double get_Delta()             { return Delta; };  // HolsteinVars: G4units propagated for Delta, E0, M.
	
	G4double get_sigma_M1()          { return sigma_M1; };
	G4double get_sigma_M2()          { return sigma_M2; };  // HolsteinVars:  G4 units propagated for both.  in MeV.
	G4double get_m_e()               { return m_e; }; // HolsteinVars:  G4units propagated for m_e.
	
	double get_mu_parent()           { return mu_parent; };
	double get_sigma_mu_parent()     { return sigma_mu_parent; };
	double get_mu_daughter()         { return mu_daughter; };
	double get_sigma_mu_daughter()   { return sigma_mu_daughter; };
	
	double get_quad_parent()         { return quad_parent; };
	double get_sigma_quad_parent()   { return sigma_quad_parent; }; 
	double get_quad_daughter()       { return quad_daughter; };
	double get_sigma_quad_daughter() { return sigma_quad_daughter; };
	
	double get_deltaC()              { return deltaC; };
	double get_sigma_deltaC()        { return sigma_deltaC; };
	
	// coupling constants (Ian Towner).
	double get_g_V()  { return g_V; };
	double get_g_A()  { return g_A; };
	double get_g_II() { return g_II; };
	double get_g_S()  { return g_S; };
	double get_g_P()  { return g_P; };
	double get_g_M()  { return g_M; };
	// lower-case decay parameters:
	double get_a1() { return a1; }; 
	double get_a2() { return a2; };
	double get_c1() { return c1; };
	double get_c2() { return c2; };
	
	double get_b()  { return b; };
	double get_d()  { return d; };
	double get_e()  { return e; };
	double get_f()  { return f; };
	double get_g()  { return g; };
	double get_h()  { return h; };
	
	double get_j2() { return j2; };
	double get_j3() { return j3; };
	
	// matrix elements (Ian Towner).
	double get_M_F()     { return M_F; };
	double get_M_GT()    { return M_GT; };
	
	double get_M_r2()    { return M_r2; };
	double get_M_sr2()   { return M_sr2; };
	double get_M_1y()    { return M_1y; };
	double get_M_sL()    { return M_sL; };
	double get_M_srp()   { return M_srp; };
	double get_M_rdotp() { return M_rdotp; }; 
	double get_M_rp()    { return M_rp; };
	double get_M_2y()    { return M_2y; };
	double get_M_3y()    { return M_3y; };
	double get_M_L()     { return M_L; };
	double get_M_Q()     { return M_Q; };


//	// JTW parameters, with no isospin corrections.
//	double jtw_a1, jtw_c1;
	
};


#endif
