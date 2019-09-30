// Author: Melissa Anholm - 2019

#ifndef K37SublevelPops_h
#define K37SublevelPops_h 1


#include <string>  
#include <sstream> 


using std::string;
using std::vector;


struct K37SublevelPopulations  // Spin 3/2 only!!
{
public:
	K37SublevelPopulations();
	K37SublevelPopulations(int the_sigma);  // completely polarized, so the_sigma=+/-1.
	~K37SublevelPopulations();
	
	double get_pop(string, int, int);
	void set_pop(string, int, int, double);
	void print_pops();
	
	//
	double get_P();
	double get_T();
	double get_Mz();
	double get_Mz2();
	double get_Mz3();
	
public: // functions that already lived in G4, and are included here for backward compatibility of some sort.
	double GetPolarization()  { return this->get_P(); }
	double GetAlignment()     { return this->get_T(); }
	
	void renormalize();
	void set_sigma_plus();
	void set_sigma_minus();
	int get_sigma();
	
private:
	bool is_sigma_plus;
	void swap_states();  // swap *only* the states, not the 'sigma' flag.  This MUST BE PRIVATE!
	
	bool sanity(int F, int M_F);
	double get_scale(string, int, int);  // how to scale populations to get M_z, M_z2, M_z3.
	
private:
	vector<double> ground_F1;
	vector<double> ground_F2;
	vector<double> excited_F1;
	vector<double> excited_F2;
	
	
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
