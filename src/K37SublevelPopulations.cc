// Author: Melissa Anholm - 2019

#include <iostream>  // cout, endl
#include <iomanip>   // setw

#undef NDEBUG
#include<assert.h>

#include "G4UIsession.hh"  // G4cout, G4endl;

#include "K37SublevelPopulations.hh"

using std::cout;
using std::endl;

// ------------------------------------------------------------- //
K37SublevelPopulations::K37SublevelPopulations()
{
	K37SublevelPopulations(1);
}

K37SublevelPopulations::K37SublevelPopulations(int the_sigma)
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
K37SublevelPopulations::~K37SublevelPopulations()
{
	G4cout << "Deleting the sublevel populations. " << G4endl;
}

// ------------------------------------------------------------- //
void K37SublevelPopulations::swap_states()
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
int K37SublevelPopulations::get_sigma() // returns +/- 1 (or 0 if it's broken)
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

void K37SublevelPopulations::set_sigma_plus()
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
void K37SublevelPopulations::set_sigma_minus()
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

double K37SublevelPopulations::get_Mz()
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
double K37SublevelPopulations::get_Mz2()
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
double K37SublevelPopulations::get_Mz3()
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
double K37SublevelPopulations::get_P()
{
	// spin 3/2 only
	double P = get_Mz()/(3.0/2.0);
	return P;
}
double K37SublevelPopulations::get_T()
{
	// spin 3/2 only
	double T = 5.0/4.0 - get_Mz2();
	return T;
}

void K37SublevelPopulations::set_pop(string level, int F, int M_F, double the_pop)
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
double K37SublevelPopulations::get_pop(string level, int F, int M_F)
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

bool K37SublevelPopulations::sanity(int F, int M_F)
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
double K37SublevelPopulations::get_scale(string the_parameter, int F, int M_F)
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
void K37SublevelPopulations::renormalize()
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
void K37SublevelPopulations::print_pops()
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
