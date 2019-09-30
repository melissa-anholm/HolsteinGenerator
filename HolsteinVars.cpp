#include <string>  // for SS
#include <sstream> // for SS
#include <vector>  // for SS

#include "SplitString.hh"
#include "K37SublevelPopulations.cc"
#include "IsotopeValues.hh"

//#include <CLHEP/Units/SystemOfUnits.h>
#include <G4SystemOfUnits.hh>
	
using std::cout;
using std::endl;
using std::string;
using std::vector;

const double pi = std::atan(1.0)*4.0;

void check_if_NDEBUG()
{
#if defined(NDEBUG) 
	std::cout << "NDEBUG is defined!" << std::endl;
#else
	std::cout << "NDEBUG is not defined." << std::endl;
#endif
	return;
}

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

/*
// namespace SS is copypasta from "SplitString.hh".
namespace SS
{
	// below:  does this one *ever* get used?! ...yes, yes it does.
	std::vector<std::string> & split(const std::string &s, const char &delim, std::vector<std::string> &elems)
	{
		std::stringstream ss(s);
		std::string item;
		while (std::getline(ss, item, delim))
		{
			elems.push_back(item);
		}
		return elems;
	}
	
	// below:  this one is used.
	std::vector<std::string> split(const std::string &s, const char &delim)
	{
		bool verbose = false;
		std::vector<std::string> elems;
		elems = SS::split(s, delim, elems);
		if(verbose)
		{
			cout << "string:  " << s << endl;  
			int i_size = elems.size();
			for(int i=0; i<i_size; i++)
			{
				cout << "\telems[" << i << "] = " << elems[i] << endl;
			}
		}
		return elems;
	}
}
*/

/*
// class isotope_values is copypasta from Isotope.hh.  It's for loading up the inputs.
class isotope_values
{
public:
	isotope_values(const double &value_, const double &uncertainty_, const std::string &message_, const std::string &name_)
	:value(value_), uncertainty(uncertainty_), message(message_), name(name_)
	{};
	
	double GetValue()const       { return value; };
	double GetUncertainty()const { return uncertainty; };
	string GetMessage()const     { return message; };     // extra added by MJA.
	string GetName()const        { return name; };        // extra added by MJA.
	
	void SetUncertainty(const double &uncertainty_)
		{ uncertainty = uncertainty_; };
	void SetValue(const double &value_)
		{ value = value_; };
	void Print() 
	{
		std::cout << std::setw(30) << name;
		std::cout << std::setw(13) << value << std::setw(12) << uncertainty;
		std::cout << "    " << message << std::endl;
	};
	
private:
	double value;
	double uncertainty;
	std::string message;
	std::string name;  // added by MJA.  It probably isn't necessary.  the info is elsewhere.
};
*/

// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- //

/*
struct sublevel_populations  // Spin 3/2 only!!
{
public:
	sublevel_populations(int the_sigma);
	double get_pop(string, int, int);
	void set_pop(string, int, int, double);
	void print_pops();
	
	double get_Mz();
	double get_Mz2();
	double get_Mz3();
	
	double get_P();
	double get_T();
	
	void renormalize();
	void set_sigma_plus();
	void set_sigma_minus();
	int get_sigma();
	
private:
	bool is_sigma_plus;
	void swap_states();  // swap *only* the states, not the 'sigma' flag.  This MUST BE PRIVATE!
	
	bool sanity(int F, int M_F);
	double get_scale(string, int, int);  // how to scale populations to get M_z, M_z2, M_z3.
	std::vector<double> ground_F1;
	std::vector<double> ground_F2;
	std::vector<double> excited_F1;
	std::vector<double> excited_F2;
};
sublevel_populations::sublevel_populations(int the_sigma)
{
	excited_F1.push_back(0);
	excited_F1.push_back(0);
	excited_F1.push_back(0);
	
	excited_F2.push_back(0);
	excited_F2.push_back(0);
	excited_F2.push_back(0);
	excited_F2.push_back(0);
	excited_F2.push_back(0);
	
	ground_F1.push_back(0);
	ground_F1.push_back(0);
	ground_F1.push_back(0);
	
	ground_F2.push_back(0);
	ground_F2.push_back(0);
	ground_F2.push_back(0);
	ground_F2.push_back(0);
	ground_F2.push_back(1);
	
	if(the_sigma==1)       { is_sigma_plus = true; }
	else if(the_sigma==-1) { is_sigma_plus = false; }
	else                   { cout << "Bad sigma!!" << endl; assert(0); return; }
	
	// if it's sigma-, swap everything.
	if(!is_sigma_plus) { swap_states(); }
}
void sublevel_populations::swap_states()
{
	bool verbose = true;
	
	// clone the 4 vectors...
	vector<double> tmp_eF1(excited_F1);
	vector<double> tmp_eF2(excited_F2);
	vector<double> tmp_gF1(ground_F1);
	vector<double> tmp_gF2(ground_F2);
	
	// F=2:
	for(int i=0; i<5; i++)
	{
	//	cout << "F=2;  i=" << i << endl;
		excited_F2[i] = tmp_eF2[4-i];
		ground_F2[i]  = tmp_gF2[4-i];
	}
	// F=1
	for(int i=0; i<3; i++)
	{
	//	cout << "F=1;  i=" << i << endl;
		excited_F1[i] = tmp_eF1[2-i];
		ground_F1[i]  = tmp_gF1[2-i];
	}
	
	if(verbose)
	{
		cout << "printing pops from within swap_states():" << endl;
		print_pops();
	}
}
int sublevel_populations::get_sigma() // returns +/- 1 (or 0 if it's broken)
{
	bool physics_sigma_plus;
	if(get_Mz()>=0) { physics_sigma_plus=true; }
	else{ physics_sigma_plus=false; }
	
	if(is_sigma_plus && physics_sigma_plus)        { return 1.0; }
	else if(!is_sigma_plus && !physics_sigma_plus) { return -1.0; }
	else
	{
		cout << "Bad!  Physical polarization direction doesn't match up with the records!" << endl;
		cout << "record:  is_sigma_plus=" << is_sigma_plus << endl;
		cout << "physics: physics_sigma_plus=" << physics_sigma_plus << endl;
		cout << "     Mz: " << get_Mz() << endl;
		assert(0); // do I need a hard kill?  The code could go on, but it probbly shouldn't.  
		return 0;
	}
}

void sublevel_populations::set_sigma_plus()
{
	bool verbose = true;
	if(get_sigma()>0) 
	{ 
		if(verbose)
		{
			cout << "Called set_sigma_plus(), but there's nothing to be done." << endl;
		}
	} // don't do anything!
	else if(get_sigma()<0)
	{
		if(verbose)
		{
			cout << "Called set_sigma_plus().  Swapping states." << endl;
		}
		is_sigma_plus=true;
		swap_states();
	}
	else
	{
		cout << "Welp.  That's bad." << endl;
		assert(0);
	}
}
void sublevel_populations::set_sigma_minus()
{
	bool verbose = true;
	if(get_sigma()>0) 
	{
		if(verbose)
		{
			cout << "Called set_sigma_minus().  Swapping states." << endl;
		}
		is_sigma_plus=false;
		swap_states();
	}
	else if(get_sigma()<0) 
	{ 
		if(verbose)
		{
			cout << "Called set_sigma_minus(), but there's nothing to be done." << endl;
		}
	} // don't do anything!	
	else
	{
		cout << "Welp.  That's bad." << endl;
		assert(0);
	}
}

double sublevel_populations::get_Mz()
{
	bool verbose=false;
	renormalize();
	
	double running_Mz = 0.0;
	
	double the_pop   = 0.0;
	double the_scale = 0.0;

	//
	int F=2;
	for(int M_F=-2; M_F<=2; M_F++)
	{
		the_pop     = get_pop("excited", F, M_F);
		the_scale   = get_scale("M_z", F, M_F);
		running_Mz += the_pop*the_scale;
		//
		the_pop     = get_pop("ground", F, M_F);
		running_Mz += the_pop*the_scale;
		if(verbose)
		{
			cout << "F=" << F << "\tM_F=" << M_F << endl;
			cout << "\tpop = " << the_pop << "\tscale=" << the_scale << endl;
		}
	}
	F=1;
	for(int M_F=-F; M_F<=F; M_F++)
	{
		the_pop     = get_pop("excited", F, M_F);
		the_scale   = get_scale("M_z", F, M_F);
		running_Mz += the_pop*the_scale;
		//
		the_pop     = get_pop("ground", F, M_F);
		running_Mz += the_pop*the_scale;
	}
	
	return running_Mz;
}
double sublevel_populations::get_Mz2()
{
	renormalize();
	double running_Mz2 = 0.0;
	
	double the_pop   = 0.0;
	double the_scale = 0.0;

	//
	int F=2;
	for(int M_F=-2; M_F<=2; M_F++)
	{
		the_pop     = get_pop("excited", F, M_F);
		the_scale   = get_scale("M_z2", F, M_F);
		running_Mz2 += the_pop*the_scale;
		//
		the_pop     = get_pop("ground", F, M_F);
		running_Mz2 += the_pop*the_scale;
	}
	F=1;
	for(int M_F=-F; M_F<=F; M_F++)
	{
		the_pop     = get_pop("excited", F, M_F);
		the_scale   = get_scale("M_z2", F, M_F);
		running_Mz2 += the_pop*the_scale;
		//
		the_pop     = get_pop("ground", F, M_F);
		running_Mz2 += the_pop*the_scale;
	}
	
	return running_Mz2;
}
double sublevel_populations::get_Mz3()
{
	renormalize();
	double running_Mz3 = 0.0;
	
	double the_pop   = 0.0;
	double the_scale = 0.0;

	//
	int F=2;
	for(int M_F=-2; M_F<=2; M_F++)
	{
		the_pop     = get_pop("excited", F, M_F);
		the_scale   = get_scale("M_z3", F, M_F);
		running_Mz3 += the_pop*the_scale;
		//
		the_pop     = get_pop("ground", F, M_F);
		running_Mz3 += the_pop*the_scale;
	}
	F=1;
	for(int M_F=-F; M_F<=F; M_F++)
	{
		the_pop     = get_pop("excited", F, M_F);
		the_scale   = get_scale("M_z3", F, M_F);
		running_Mz3 += the_pop*the_scale;
		//
		the_pop     = get_pop("ground", F, M_F);
		running_Mz3 += the_pop*the_scale;
	}
	
	return running_Mz3;
}
double sublevel_populations::get_P()
{
	// spin 3/2 only
	double P = get_Mz()/(3.0/2.0);
	return P;
}
double sublevel_populations::get_T()
{
	// spin 3/2 only
	double T = 5.0/4.0 - get_Mz2();
	return T;
}

void sublevel_populations::set_pop(string level, int F, int M_F, double the_pop)
{
	// sanity checks first.
	if( level.compare(string("ground"))!=0 && level.compare(string("excited"))!=0 )
	{
		cout << "Level:  " << level << " is not a thing.  Try again." << endl;
		assert(0);
		return;
	}
	if( !sanity(F, M_F) ) { assert(0); return; } 
	
	// now actually set it.
	// Now actually do the thing.
//	std::vector<double> use_this_vector;
	int offset;
	if(F==1) { offset=1; }
	if(F==2) { offset=2; }
	
	if( !level.compare(string("ground")) ) // level is ground level
	{
		if(F==1)
		{
		//	use_this_vector = ground_F1;
			ground_F1[M_F+offset] = the_pop;
		}
		else if(F==2)
		{
		//	use_this_vector = ground_F2;
			ground_F2[M_F+offset] = the_pop;
		}
	}
	else if( !level.compare(string("excited")) ) // level is excited level
	{
		if(F==1)
		{
		//	use_this_vector = excited_F1;
			excited_F1[M_F+offset] = the_pop;
		}
		else if(F==2)
		{
		//	use_this_vector = excited_F2;
			excited_F2[M_F+offset] = the_pop;
		}
	}
	
//	use_this_vector[M_F+offset] = the_pop;
	return;
}
double sublevel_populations::get_pop(string level, int F, int M_F)
{
//	cout << "Getting population." << endl;
	// sanity checks first.
	if( !level.compare(string("ground")) && !level.compare(string("excited")) )
	{
		cout << "Level:  " << level << " is not a thing.  Try again." << endl;
		assert(0);
		return -1.0;
	}
	if( !sanity(F, M_F) ) { assert(0); return -1.0; }
	
	// Now actually do the thing.
	std::vector<double> use_this_vector;
	int offset;
	double the_pop;
	if( !level.compare(string("ground")) ) // level is ground level
	{
		if(F==1)
		{
			use_this_vector = ground_F1;
		}
		else if(F==2)
		{
			use_this_vector = ground_F2;
		}
	}
	else if( !level.compare(string("excited")) ) // level is excited level
	{
		if(F==1)
		{
			use_this_vector = excited_F1;
		}
		else if(F==2)
		{
			use_this_vector = excited_F2;
		}
	}
	if(F==1) { offset=1; }
	if(F==2) { offset=2; }
	
	the_pop = use_this_vector[M_F+offset];
	return the_pop;
}

bool sublevel_populations::sanity(int F, int M_F)
{
	if( F!=1 && F!=2)
	{
		cout << "F = " << F << " isn't a thing.  Try again."  << endl;
		return false;
	}
	if( M_F>F || M_F<-F )
	{
		cout << "M_F = " << M_F << " (F = " << F << ") is not a thing.  Try again." << endl;
		return false;
	}
	return true;
}
double sublevel_populations::get_scale(string the_parameter, int F, int M_F)
{
	if( the_parameter.compare(string("M_z")) != 0
	 && the_parameter.compare(string("M_z2")) != 0
	 && the_parameter.compare(string("M_z3")) != 0 )
	{
		cout << "The function isn't made for this!  Input should be \"M_z\", \"M_z2\", or \"M_z3\"" << endl;
		assert(0);
		return -1;
	}
	if( !sanity(F, M_F) ) { assert(0); return -1; }
	
	//
	double the_scale = 0.0;
	if( the_parameter.compare(string("M_z")) == 0 )
	{
		if(F==2)
		{
			if(M_F==-2) { the_scale = -3.0/2.0; }
			if(M_F==-1) { the_scale = -3.0/4.0; }
			if(M_F==0)  { the_scale =  0.0;     }
			if(M_F==1)  { the_scale =  3.0/4.0; }
			if(M_F==2)  { the_scale =  3.0/2.0; }
		}
		if(F==1)
		{
			if(M_F==-1) { the_scale = -5.0/4.0; }
			if(M_F==0)  { the_scale =  0.0;     }
			if(M_F==1)  { the_scale =  5.0/4.0; }
		}
	}
	if( the_parameter.compare(string("M_z2")) == 0 )
	{
		if(F==2)
		{
			if(M_F==-2) { the_scale = 9.0/4.0; }
			if(M_F==-1) { the_scale = 3.0/4.0; }
			if(M_F==0)  { the_scale = 1.0/4.0; }
			if(M_F==1)  { the_scale = 3.0/4.0; }
			if(M_F==2)  { the_scale = 9.0/4.0; }
		}
		if(F==1)
		{
			if(M_F==-1) { the_scale = 7.0/4.0; }
			if(M_F==0)  { the_scale = 1.0/4.0; }
			if(M_F==1)  { the_scale = 7.0/4.0; }
		}
	}
	if( the_parameter.compare(string("M_z3")) == 0 )
	{
		if(F==2)
		{
			if(M_F==-2) { the_scale = -27.0/8.0; }
			if(M_F==-1) { the_scale = -15.0/16.0;}
			if(M_F==0)  { the_scale =  0.0;      }
			if(M_F==1)  { the_scale =  15.0/16.0;}
			if(M_F==2)  { the_scale =  27.0/8.0; }
		}
		if(F==1)
		{
			if(M_F==-1) { the_scale = -41.0/16.0;}
			if(M_F==0)  { the_scale =  0.0;      }
			if(M_F==1)  { the_scale =  41.0/16.0;}
		}
	}
	return the_scale;
}
void sublevel_populations::renormalize()
{
	double running_sum = 0.0;
	for(int i=0; i<(int)excited_F1.size(); i++)
	{
		running_sum += excited_F1[i];
	}
	for(int i=0; i<(int)excited_F2.size(); i++)
	{
		running_sum += excited_F2[i];
	}
	for(int i=0; i<(int)ground_F1.size(); i++)
	{
		running_sum += ground_F1[i];
	}
	for(int i=0; i<(int)ground_F2.size(); i++)
	{
		running_sum += ground_F2[i];
	}
	
	// sanity check:
	if(running_sum==0)
	{
		cout << "That's bad.  Total population is zero." << endl;
		assert(0);
		return;
	}
	if(running_sum<0)
	{
		cout << "That's bad.  Total population is negative." << endl;
		assert(0);
		return;
	}
	
	// run through them all again.
	for(int i=0; i<(int)excited_F1.size(); i++)
	{
		excited_F1[i] /= running_sum;
	}
	for(int i=0; i<(int)excited_F2.size(); i++)
	{
		excited_F2[i] /= running_sum;
	}
	for(int i=0; i<(int)ground_F1.size(); i++)
	{
		ground_F1[i] /= running_sum;
	}
	for(int i=0; i<(int)ground_F2.size(); i++)
	{
		ground_F2[i] /= running_sum;
	}
	return;
}
void sublevel_populations::print_pops()
{
	int the_width=12;
	int the_precision=the_width-4;
	
	cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
	cout << "Excited State Populations for F=2:  " << endl;
	cout << std::setw(the_width) << std::left << "M_F=-2";
	cout << std::setw(the_width) << std::left << "M_F=-1";
	cout << std::setw(the_width) << std::left << "M_F=0";
	cout << std::setw(the_width) << std::left << "M_F=+1";
	cout << std::setw(the_width) << std::left << "M_F=+2" << endl;
	
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_pop("excited", 2,-2);
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_pop("excited", 2,-1);
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_pop("excited", 2,0);
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_pop("excited", 2,1);
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_pop("excited", 2,2) << endl;

	cout << "Excited State Populations for F=1:  " << endl;
	cout << std::setw(the_width) << " ";
	cout << std::setw(the_width) << std::left << "M_F=-1";
	cout << std::setw(the_width) << std::left << "M_F=0";
	cout << std::setw(the_width) << std::left << "M_F=+1" << endl;
	
	cout << std::setw(the_width) << " ";
	cout << std::setprecision(the_precision) << std::setw(the_width) << get_pop("excited", 1,-1);
	cout << std::setprecision(the_precision) << std::setw(the_width) << get_pop("excited", 1,0);
	cout << std::setprecision(the_precision) << std::setw(the_width) << get_pop("excited", 1,1) << endl;

	cout << endl;
	cout << "Ground State Populations for F=2:  " << endl;
	cout << std::setw(the_width) << std::left << "M_F=-2";
	cout << std::setw(the_width) << std::left << "M_F=-1";
	cout << std::setw(the_width) << std::left << "M_F=0";
	cout << std::setw(the_width) << std::left << "M_F=+1";
	cout << std::setw(the_width) << std::left << "M_F=+2" << endl;
	
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_pop("ground", 2,-2);
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_pop("ground", 2,-1);
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_pop("ground", 2,0);
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_pop("ground", 2,1);
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_pop("ground", 2,2);
	cout << endl;

	cout << "Ground State Populations for F=1:  " << endl;
	cout << std::setw(the_width) << " ";
	cout << std::setw(the_width) << std::left << "M_F=-1";
	cout << std::setw(the_width) << std::left << "M_F=0";
	cout << std::setw(the_width) << std::left << "M_F=+1" << endl;
	
	cout << std::setw(the_width) << " ";
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_pop("ground", 1,-1);
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_pop("ground", 1,0);
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_pop("ground", 1,1) << endl;
	cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
}
*/

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
	string paramfilename;
	std::string isotopeName;      // gets a value loaded into it, but probably doesn't ever get used.
	std::map<std::string, isotope_values * > theInputs;
	bool loadup_textfile();
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

HolsteinVars::HolsteinVars():
	I_spin(1.5), u(3.0/2.0), v(3.0/2.0), 
	Z_parent(19), N_parent(18), A_nucleons(37), Z_daughter(18), N_daughter(19),
	T_isospin(0), T3_parent(0), T3_daughter(0), M_F_isospin(0),
	mu_parent(0), sigma_mu_parent(0), mu_daughter(0), sigma_mu_daughter(0),
	quad_parent(0), sigma_quad_parent(0), quad_daughter(0), sigma_quad_daughter(0),
	is_twopercent(false), 
//	hbarc_eV_nm(0), 
	amu_to_mev(1), sigma_amu_to_mev(0), // uh ... these are terrible default values.
//	E0(0), M(0), Delta(0), 
//	m_e(0), sigma_m(0),            // m_electron
	deltaC(0), sigma_deltaC(0),  // isospin symmetry breaking.  somehow.
	g_V(0), g_A(0), g_P(0), g_M(0), 
	g_II(0), g_S(0), // the so-called "second-class" currents.
	M_F(0), M_GT(0), 
	M_r2(0), M_sr2(0), M_1y(0), M_sL(0), M_srp(0), M_rdotp(0), M_2y(0), 
	M_L(0), // only used for th. weak magnetism 'b'?
	M_Q(0), M_rp(0), // exp. measurement to go from 'g' to 'f'.
	a1(0), a2(0), c1(0), c2(0), 
	b(0), d(0), e(0), f(0), g(0), h(0), 
	j2(0), j3(0), 
	sigma_M1(0), sigma_M2(0),
	paramfilename(string("K_37_INPUT.txt"))
{
	loadup_textfile();                    // reads the text file into 'theInputs'.
	initialize_physics_parameters();      //
	initialize_spinfuncs(I_spin, I_spin); // sets u and v, and epsilon_uv and etc.
	
	return;
}


void HolsteinVars::initialize_physics_parameters()
{
	bool verbose = false;
	
	hbarc_eV_nm      = FindValue("HBARC");
	amu_to_mev       = FindValue("AMU_TO_MEV");
	sigma_amu_to_mev = FindUncertainty("AMU_TO_MEV");
	
	// kludge in the speed of light:
	speed_of_light = 299792458.0;// *meter/second;  // meters per second.
	
	// I_spin = u = v = 3/2 is only appropriate for the 98% branch.  Have to do something else for the 2% branch.
	I_spin     = FindValue("SPIN");
	u          = I_spin;
	v          = I_spin; 
	
	m_e        = ( FindValue("MASS_OF_ELECTRON") )*MeV;        // in MeV/c^2
//	sigma_m    = ( FindUncertainty("MASS_OF_ELECTRON") )*MeV;  // 
	
	M          = ( FindValue("AVERAGE_MASS_PARENT_DAUGHTER")*amu_to_mev + 0.5*m_e/MeV )*MeV;       // average mass. units propagated.
	
	Z_parent   = FindValue("NUMBER_OF_PROTONS");
	N_parent   = FindValue("NUMBER_OF_NEUTRONS");
	A_nucleons = Z_parent + N_parent;
	Z_daughter = FindValue("NUMBER_OF_PROTONS_DAUGHTER");  // 
	N_daughter = A_nucleons - Z_daughter;
	
	T_isospin   = std::abs((N_parent-Z_parent)/2.0);  // this is only true for isobaric analog decays.
	T3_parent   = (Z_parent-N_parent)/2.0;            // isospin sign convention:  protons are isospin-up.
	T3_daughter = (Z_daughter-N_daughter)/2.0; 
	
	mu_parent           = FindValue("PARENT_MAGNETIC_MOMENT");
	sigma_mu_parent     = FindUncertainty("PARENT_MAGNETIC_MOMENT");
	mu_daughter         = FindValue("DAUGHTER_MAGNETIC_MOMENT");
	sigma_mu_daughter   = FindUncertainty("DAUGHTER_MAGNETIC_MOMENT");
	
	quad_parent         = FindValue("PARENT_QUADRUPOLE_MOMENT");
	sigma_quad_parent   = FindUncertainty("PARENT_QUADRUPOLE_MOMENT");
	quad_daughter       = FindValue("DAUGHTER_QUADRUPOLE_MOMENT");
	sigma_quad_daughter = FindUncertainty("DAUGHTER_QUADRUPOLE_MOMENT");
	
	// Just kludge in some values, right here right now.
	sigma_M1 = 0.10e-6 * amu_to_mev * MeV;
	sigma_M2 = 0.22e-6 * amu_to_mev * MeV;

	// some calculated stuff...
	Delta = ( FindValue("MASS_OF_PARENT")*amu_to_mev - FindValue("MASS_OF_DAUGHTER")*amu_to_mev - m_e/MeV )*MeV;
	E0    = ( (Delta/MeV)*(1.0+(m_e/MeV)*(m_e/MeV)/(2.0*(M/MeV)*(Delta/MeV)) )/(1.0+(Delta/MeV)/(2.0*M)) )*MeV;
	
	double recoil_correction       = ( 1.0 + (Delta/MeV)/(2.0*(M/MeV)) ); // this *should* come out unitless...
	
	// deltaC from severijns 2008 "isospin symmetry breaking correction".  
	// it seems to go into M_F and M_GT.  only those.  what the hell is it?!  
	deltaC       = FindValue("DELTA_C");
	sigma_deltaC = FindUncertainty("DELTA_C");
	
	// Vector/Axial Coupling constants.
	g_V        = FindValue("GV");
	g_A        = FindValue("GA");
//	g_A = 0.922;  // HUGE KLUDGE!
	
	// Other misc. coupling constants...
	g_II       = FindValue("GII");        // value is zero.
	g_S        = FindValue("GS");         // it's zero.  
	g_P        = FindValue("GP");         // it's zero.
	g_M        = FindValue("GM");         // used only(?) for theoretical calculation of b_weakmagnetism.
	
	
	// matrix elements...
	M_F         = FindValue("M_F")  / sqrt(1.0-deltaC);  // corrected for *something*...  copied from g4.  makes it slightly bigger.
	M_GT        = FindValue("M_GT") / sqrt(1.0-deltaC);  // corrected for *something*...  copied from g4.  isospin mixing, apparently.  It shouldn't get a correction.
//	M_GT *= -1.0;  // KLUDGE KLUDGE KLUDGE
	
	// some other matrix elements...
	M_r2    = FindValue("M_R2")   / (hbarc_eV_nm*hbarc_eV_nm);  
	M_sr2   = FindValue("M_SR2")  / (hbarc_eV_nm*hbarc_eV_nm);  
	M_1y    = FindValue("M_1Y")   / (hbarc_eV_nm*hbarc_eV_nm);  
	M_sL    = FindValue("M_SL");                    // doesn't need to be divided??  ... it's zero anyway.
	M_srp   = FindValue("M_SRP")  /  hbarc_eV_nm;         // it's zero, and I'm dividing it anyway.
	M_rdotp = FindValue("M_RDOTP")/ (hbarc_eV_nm*hbarc_eV_nm);  // it's zero, but it's still getting divided by hbarc_eV_nm^2.
	M_2y    = FindValue("M_2Y")   / (hbarc_eV_nm*hbarc_eV_nm);  
	M_3y    = FindValue("M_3Y")   / (hbarc_eV_nm*hbarc_eV_nm);  
//	M_L     = FindValue("M_L");                     // 'b' th. calc.
	// M_Q, M_RP:  don't load values from file, set them up later from experimental 'g'.
	
	
	// calculate a1.  value copied from g4, error calculated by me.
//	a1  = g_V*M_F/recoil_correction;  // copied from g4.
	a1  = g_V*(M_F - (Delta/MeV)*(Delta/MeV)/6.0*M_r2 + (Delta/MeV)/3.0*M_rdotp) / ( 1.0 + (Delta/MeV)/(2.0*(M/MeV)) );
	a2  = g_V*(M_r2/6.0)*((M/MeV)*(M/MeV)) / recoil_correction;  // copied from g4.  matches Dan, even though it's huge.
	
	c1  = g_A*(M_GT - (Delta/MeV)*(Delta/MeV)/6.0*M_sr2 + M_1y*2.0*(Delta/MeV)*(Delta/MeV)/(6.0*sqrt(10.)) + A_nucleons*(Delta/MeV)/(2.0*M)*M_sL + (Delta/MeV)/2.0*M_srp ) / recoil_correction;
	c2  = g_A*((M_sr2/6.0) + (M_1y/(6.0*sqrt(10.))))*( ((M/MeV)*(M/MeV))/recoil_correction);
	
	M_F_isospin = sqrt( (T_isospin + T3_parent)*(T_isospin - T3_parent + 1.0) );  // use for b, g.
	
	// Holstein (21)-(22):
	b       = A_nucleons*M_F_isospin * sqrt((I_spin+1.0)/I_spin) * (mu_parent-mu_daughter);
	
	// experimental 'b' makes the PDF less negative.
	if(verbose)
	{
		cout << "M_F_isospin = " << M_F_isospin << endl;
		cout << "* experimental b = " << b << endl;
	}

	double d_I  = g_A*(-1.0*M_GT + (Delta/MeV)*(Delta/MeV)*M_sr2/6.0 + M_1y*(Delta/MeV)*((M/MeV)+(Delta/MeV)/6.0)/sqrt(10.) + (A_nucleons*M_sL) + ((M/MeV)*M_srp) ) / recoil_correction;
	double d_II = -1.0*g_II*A_nucleons*M_GT;
	d = d_I + d_II;
	
	// Eq. (22).  assumes CVC.
	e       = 0;
	
	// Use experimental quadrupole moments for 'g', rather than calculated M_Q.
	double jterm = sqrt( (I_spin + 1.0)*(2.0*I_spin + 3.0)/( I_spin*(2.0*I_spin - 1.0) ));
	g   = (2.0/3.0) * (M/MeV)*(M/MeV) * M_F_isospin * jterm*(quad_daughter-quad_parent) / (hbarc_eV_nm*hbarc_eV_nm);
	// Set M_Q by experimental 'g' value:
	M_Q = g / ( -4.0/3.0 * (M/MeV)*(M/MeV) * g_V );
	// Set M_rp by M_Q:
	M_rp = (E0/MeV)*M_Q/sqrt(6);
	
	// for f, use new value of M_rp
	f   = 2.0*g_V*(M/MeV)*M_rp;
	
	// h ...
	h   = -1.0*(g_A*2.0*((M/MeV)*(M/MeV))*M_1y/sqrt(10.) + g_P*(A_nucleons*A_nucleons)*M_GT) / recoil_correction;
	
	// j2 ..
	j2   = -g_A*2.0*((M/MeV)*(M/MeV))*M_2y/3.0;
	
	// j3 ?
	j3   = -g_A*2.0*((M/MeV)*(M/MeV))*M_3y/3.0;
	
	if(verbose)
	{
//		cout << "* experimental g = " << g << endl;
//		cout << "* experimental M_Q = " << M_Q << endl;
//		cout << "* experimental M_rp = " << M_rp << endl;
		cout << endl;
		cout << "--" << endl;
		cout << "T3_parent = " << T3_parent << " \t";
		cout << "T3_daughter = " << T3_daughter << endl;
		cout << "hbarc_eV_nm = " << hbarc_eV_nm << "\t";
		cout << "amu_to_mev = " << amu_to_mev << endl;
		cout << endl;
		cout << "recoil_correction = " << recoil_correction << endl;
		cout << "deltaC = " << deltaC << " +/- " << sigma_deltaC << endl;
		cout << endl;
	//	cout << "g_V = " << g_V << endl;
	//	cout << "M_F = " << M_F << endl;
		/*
		cout << "a1        = " << a1 << endl;
		cout << "by formula: " << g_V*(M_F - (Delta/MeV)*(Delta/MeV)/6.0*M_r2 + (Delta/MeV)/3.0*M_rdotp) / ( 1.0 + (Delta/MeV)/(2.0*(M/MeV)) ) << endl;
		cout << "prev. a1  = " << g_V*M_F/recoil_correction << endl;
		cout << "a2 = " << a2 << endl; 
		*/
	//	cout << "g_A = " << g_A << endl;
	//	cout << "M_GT = " << M_GT << endl;
		/*
		cout << "c1 = " << c1 << endl; 
		cout << "c2 = " << c2 << endl; 
		cout << "prev. c1 = " << g_A*M_GT/recoil_correction << endl;
		cout << endl;
		*/
		cout << "mu_parent = " << mu_parent << endl;
		cout << "mu_daughter = " << mu_daughter << endl;
		/*
		cout << "b = " << b << endl;
		cout << "d = " << d << endl;
		cout << "e = " << e << endl;
		cout << "f = " << f << endl;
		cout << "g = " << g << endl;
	//	cout << "M_1y = " << M_1y << endl;
	//	cout << "g_P = " << g_P << endl;
		cout << "h = " << h << endl;
		cout << "j2 = " << j2 << endl;  
		cout << "j3 = " << j3 << endl;  
		*/
		cout << "- - -" << endl;
		cout << "(M/MeV) = " << (M/MeV)  << endl;
		cout << "(Delta/MeV) = " << (Delta/MeV)  << endl;
		cout << "E0 = " << E0 << endl;
		cout << "E0/MeV = " << E0/MeV << endl;
		cout << "-- * --" << endl;
		
		/*
		cout << "M_F     = " << M_F << endl;
		cout << "M_GT    = " << M_GT << endl;
		cout << "M_r2    = " << M_r2 << endl;
		cout << "M_sr2   = " << M_sr2 << endl;
		cout << "M_1y    = " << M_1y << endl;
		cout << "M_sL    = " << M_sL << endl;
		cout << "M_srp   = " << M_srp << endl;
		cout << "M_rdotp = " << M_rdotp << endl;
		cout << "M_2y    = " << M_2y << endl;
		cout << "M_3y    = " << M_3y << endl;
		*/
		
		print_matrixelements();
		print_couplingconstants();
		print_calculatedJTW();
		print_holsteinalphabet();
	}
}

void HolsteinVars::print_vars()
{
	cout << "--" << endl;
	print_couplingconstants();
	cout << endl;
	print_matrixelements();
	cout << endl;
	print_holsteinalphabet();
	cout << endl;
	print_calculatedJTW();
	cout << endl;
	
	// print sublevel populations.
	// print spin functions?
	cout << "--" << endl;
}

void HolsteinVars::print_couplingconstants()
{
	cout << "g_V = " << g_V << endl;
	cout << "g_A = " << g_A << endl;
	cout << "g_P = " << g_P << endl;

	cout << "g_II= " << g_II << endl;
	cout << "g_S = " << g_S << endl;
	cout << "g_P = " << g_P << endl;
	cout << "g_M = " << g_M << endl;

}
void HolsteinVars::print_matrixelements()
{
	cout << "M_F     = " << M_F << endl;
	cout << "M_GT    = " << M_GT << endl;
	cout << "M_r2    = " << M_r2 << endl;
	cout << "M_sr2   = " << M_sr2 << endl;
	cout << "M_1y    = " << M_1y << endl;
	cout << "M_sL    = " << M_sL << endl;
	cout << "M_srp   = " << M_srp << endl;
	cout << "M_rdotp = " << M_rdotp << endl;
	cout << "M_2y    = " << M_2y << endl;
	cout << "M_3y    = " << M_3y << endl;
	//
	cout << "M_Q     = " << M_Q  << "\t(experimental, from g)" << endl;
	cout << "M_rp    = " << M_rp << "\t(experimental, from M_Q)" << endl;
}

void HolsteinVars::print_calculatedJTW()
{
	double rho = (g_A*M_GT)/(g_V*M_F);
//	double Abeta = -2.0*rho*(sqrt(3.0/5.0) - abs(rho)/5.0) / (1.0+rho*rho);
	
	cout << "rho   = (g_A*M_GT)/(g_V*M_F)" << endl;
	cout << "      = " << (g_A*M_GT)/(g_V*M_F) << endl;
//	cout << "Abeta = -2*rho*(sqrt(3/5) - rho/5) / (1+rho^2) [from proposal]" << endl;
//	cout << "      = " << -2.0*rho*(sqrt(3.0/5.0) - rho/5.0) / (1.0+rho*rho) << endl;
	cout << "Abeta = (2/5*rho^2 - 2*sqrt(3/5)*|rho|) / (rho^2+1) [from calculation]" << endl;
	cout << "      = " << (2.0/5.0*rho*rho - 2.0*sqrt(3.0/5.0)*abs(rho) ) / (rho*rho+1.0) << endl;
	cout << std::setprecision(8);
	cout << "(E0/MeV) = " << (E0/MeV) << endl;
}

void HolsteinVars::print_holsteinalphabet()
{
	// debug this stuff...
	cout << "a1 = " << a1 << endl;
	cout << "a2 = " << a2 << endl; 
	cout << "c1 = " << c1 << endl; 
	cout << "c2 = " << c2 << endl; 
	//
	cout << "b  = " << b << "\t(experimental magnetic moment)" << endl;
	cout << "d  = " << d << endl;
	cout << "e  = " << e << endl;
	cout << "f  = " << f << endl;
	cout << "g  = " << g << "\t(experimental quadrupole)" << endl;
	cout << "h  = " << h << endl;
	cout << "j2 = " << j2 << endl;  
	cout << "j3 = " << j3 << endl;  
}

bool HolsteinVars::initialize_spinfuncs(double u_, double v_)
{
	u = u_;
	v = v_;
	// I'm not going to bother to code these up for any values *other* than u=v=3/2.
	// ie, the 98% branch of 37K decay.
	// Maybe more later.
	assert( u == 3.0/2.0 && v == 3.0/2.0 && !is_twopercent ); // hard kill.
	if( u == 3.0/2.0 && v == 3.0/2.0 && !is_twopercent)
	{
		delta_uv   =  1.0;
		gamma_uv   =  1.0;
		lambda_uv  = -sqrt(2.0)/5.0;          // 2.0*sqrt(3.0);
		theta_uv   =  1.0;
		kappa_uv   =  1.0/(2.0*sqrt(2));      // 3.0/sqrt(2.0);
		epsilon_uv = -1.0/( 2.0*sqrt(5.0));   // -1.0/(2.0*sqrt(35.0));
		rho_uv     = -41.0/40.0;
		sigma_uv   = -41.0/(4.0*sqrt(35.0));  // 41.0/(4.0*sqrt(35.0));
		phi_uv     =  0.0;                    // 1.0/32.0 * sqrt(3.0/5.0);
	}
	else
	{
		cout << "ERROR:  Could not initialize the Holstein spin functions!!" << endl;
		cout << "        Everything will be broken now." << endl;
		//
		delta_uv   = 0;
		gamma_uv   = 0;
		lambda_uv  = 0;
		theta_uv   = 0;
		kappa_uv   = 0;
		epsilon_uv = 0;
		rho_uv     = 0;
		sigma_uv   = 0;
		phi_uv     = 0;
		
		return false;
	}
	return true;
}

bool HolsteinVars::loadup_textfile()
// puts the contents of the text file into 'theInputs'.
{
	bool verbose = false;
	
	std::fstream inputfile(paramfilename.c_str(), std::fstream::in);
	if(verbose)
	{
		if(inputfile) { cout << "paramfilename = " << paramfilename << " -- file opened." << endl; }
		else          { cout << "Could not open file:  " << paramfilename << endl; }
	}
	if(inputfile)
	{
		std::string line;
		std::vector<std::string> parsed;
		int lineNumber = 1;
		while(getline( inputfile, line ))
		{
			// don't even bother if the line starts with a "#".
			if( line.find_first_of("#") == 0 )
			{
				if(verbose) { cout << "*Line " << lineNumber << " -- skipping" << endl; }
				continue;
			}
			parsed = SS::split(line, ':');
			switch(parsed.size())
			{
				case 1:
				{
					cout << "Problem with input line: " << lineNumber << endl;
					break;
				}
				case 2:
				{
					isotopeName = parsed[0];
					if(verbose)
					{
						cout << "Case 2:" << endl;
						cout << "parsed[0] = " << parsed[0] << endl;
					}
					break;
				}
				case 3:
				{
					if (verbose)
					{
						cout << "before erasing, parsed[2] = " << parsed[2] << endl;
					}
					parsed[2].erase(parsed[2].begin(), parsed[2].begin() + parsed[2].find_first_of("#"));
					theInputs[parsed[0]]= new isotope_values(std::stod(parsed[1]), 0, parsed[2], parsed[0]);
					// in Case 3, parsed[0] is the name, parsed[1] is the value, 
					//  	0 is the (implied) uncertainty, and parsed[2] is the comment.
					if(verbose)
					{
						cout << "Case 3:" << "\t";
						cout << "parsed[0] = " << parsed[0] << endl;
						cout << "\tparsed[1] = " << parsed[1] << "\tparsed[2] = " << parsed[2];
						cout << endl;
					}
					break;
				}
				case 4:
				{
					parsed[3].erase(parsed[3].begin(), parsed[3].begin() + parsed[3].find_first_of("#"));
					theInputs[parsed[0]]= new isotope_values( std::stod(parsed[1]), std::stod(parsed[2]), parsed[3], parsed[0]);
					// In case 4, parsed[0] is the name, parsed[1] is the value, 
					//  	parsed[2] is the uncertainty, and parsed[3] is the comment.
					if(verbose)
					{
						cout << "Case 4:" << "\t";
						cout << "parsed[0] = " << parsed[0] << endl;
						cout << "\tparsed[1] = " << parsed[1] << "\tparsed[2] = " << parsed[2] << "\tparsed[3] = " << parsed[3];
						cout << endl;
					}
					break;
				}
				default:
				{
					cout << "Problem with input line: " << lineNumber << endl;
					break;
				}
			}
			parsed.clear();
			++lineNumber;
		}
	}
	else
	{
		std::cerr << "File: " << paramfilename << " could not be opened." << std::endl;
		assert(0);
	}
	inputfile.close();	
	// check if we have all of the parameters we need?
	return true;
}
double HolsteinVars::FindValue(const std::string &key_) const
{
	auto findIt = theInputs.find(key_);
	if(findIt != theInputs.end()) { return findIt->second->GetValue(); }
	else 
	{ 
		cout << "* Couldn't find entry for " << key_ << " (value)" << endl;
		assert(0);
		return 0; 
	}
}
double HolsteinVars::FindUncertainty(const std::string &key_) const
{
	auto findIt = theInputs.find(key_);
	
	if(findIt != theInputs.end()) { return findIt->second->GetUncertainty(); }
	else 
	{
		cout << "* Couldn't find entry for " << key_ << " (uncertainty)" << endl;
		assert(0);  // kill.
		return 0;
	}
}
void HolsteinVars::print_isotope_values()
// for debugging.  Only for values loaded from the file.
{
	for (auto themap : theInputs) 
		{ themap.second -> Print(); }
	return;
}


