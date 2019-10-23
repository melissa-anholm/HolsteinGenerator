// Author: Melissa Anholm - 2019

#ifndef K37SublevelPops_h
#define K37SublevelPops_h 1


#include <string>  
#include <sstream> 
#include <vector>
#include <map>
#include "fstream"   // 

#include "IsotopeValues.hh"

using std::string;
using std::vector;
using std::map;


struct K37SublevelPopulations  // Spin 3/2 only!!
{
public:
//	K37SublevelPopulations();
//	K37SublevelPopulations(int the_sigma);  // completely polarized, so the_sigma=+/-1.
	K37SublevelPopulations(int the_sigma=1);  // completely polarized, so the_sigma=+/-1.
	~K37SublevelPopulations();
	
	double get_pop(string, int, int);
	void set_pop(string, int, int, double);
	void print_pops();
	void print_moments();
	
	//
	double get_P();
	double get_T();
	double get_Mz();
	double get_Mz2();
	double get_Mz3();
	
public: // functions that already lived in G4, and are included here for backward compatibility of some sort.
	double GetPolarization()  { return this->get_P(); }
	double GetAlignment()     { return this->get_T(); }
	
	void renormalize(bool verbose=false);
	
	void setup_sigma();
	void set_sigma_plus();
	void set_sigma_minus();
	int get_sigma();  // checks both physical- and recorded polarization directions.  Might break for M_z=0.
	
	void AdjustPolarization(double);  // not implemented yet!
	
private:
	void Setup_Pops_From_InputsMap();
	string atomic_filename;
	map<string, isotope_values * > theInputs;
	bool loadup_textfile(string paramfilename);
	double FindValue(const std::string &key_) const;
	double FindUncertainty(const std::string &key_) const;
	void print_isotope_values();  // for debugging.  prints from 'theInputs'.
	//
	

private:	
	bool sanity(int F, int M_F);   // Checks that F==1 || F==2, and |M_F| <= F.
	double get_scale(string, int, int);  // how to scale populations to get M_z, M_z2, M_z3.
	
private:
	vector<double> ground_F1;
	vector<double> ground_F2;
	vector<double> excited_F1;
	vector<double> excited_F2;
	

public:
void Setup_Pops_From_InputsMap( map<string, isotope_values * > theInputs);
	
	/*
	// failed functions:
public:
	void Setup_FromPolarizationOnly(double pol);  // use measured laser powers for OP:  <normal/repumper> to estimate NG2/NG1.  alignment comes out *almost* within 1 sigma experimental uncertainty.  
	void Setup_FromPolarizationAlignment(double pol, double ali);
private:
	void killall_pops();  // must be private!!
	void SetPops_Ns_NG1_NG2(double Ns, double NG1, double NG2, int sigma);
	*/
	
	//
private:
	double allowed_mismatch;
	double op_power_ratio;  // <P_op/P_re> = <NG2/NG2>, from laser power info with our toy model. 1.7557 +/- 0.2898 in 2014.
	bool is_sigma_plus;
	void swap_states();  // swap *only* the states, not the 'sigma' flag.  This MUST BE PRIVATE!
	
	
	// things that are in the G4 code version of this file, but that aren't actually implemented fully yet.  They're just here for reading...
	/*
public:
	// Setup Multipoles:  (these functions don't even work yet.)
	void Setup_FromPolarizationOnly(double pol);
	void Setup_FromPolarizationAlignment(double pol, double ali);
	void Setup_FromPolarizationAlignmentOctopole(double pol, double ali, double oct);  // this method is stupid in that it mixes conventions.
	void Setup_FromDipoleQuadrupoleOctopole(double dip, double quad, double oct);
	
	//
	void Initialize() { return; }  // currently does nothing.
	
private:
// these are private to prevent the user from calling them directly.
// Instead, messenger classes will use SetPolarizationAlignment(...)
// and SetPolarizationOnly(...).
// Actually, we'll make our instantiation of this whole class private.  These things can be public again.  this is probably a bad coding practice.
private:
	void SetPolarizationOnly(double pol);
	void SetAlignmentOnly(double ali);
	void SetOctopoleOnly(double oct);
	void Set_Mz_Only(double new_Mz);
	void Set_Mz2_Only(double new_Mz2);
	void Set_Mz3_Only(double new_Mz3);
public:	
// Utility classes...
	void MakeAlignmentFromPolarization();
	void MakeOctopoleFromPolarizationAlignment();
	*/
	
};


#endif
