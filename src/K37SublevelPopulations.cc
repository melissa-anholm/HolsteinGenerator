// Author: Melissa Anholm - 2019

#include <iostream>  // cout, endl
#include <iomanip>   // setw

#undef NDEBUG
#include<assert.h>

#include "G4UIsession.hh"  // G4cout, G4endl;

#include "K37SublevelPopulations.hh"
#include "K37Config.hh"

using std::cout;
using std::endl;

// ------------------------------------------------------------- //
/*
K37SublevelPopulations::K37SublevelPopulations()
{
	cout << "Called K37SublevelPopulations() [without parameters!]." << endl;
	
	K37SublevelPopulations(1);  // this thing.  this thing is what fucks it up.
	cout << "K37SublevelPopulations() is done creating K37SublevelPopulations(1)." << endl;
}
*/

K37SublevelPopulations::K37SublevelPopulations(int the_sigma) :  // fully polarized up or down.
	allowed_mismatch(1.0e-15),  // 1e-15 works well for Wisely.
	op_power_ratio(1.7557)//,     // 1.7557 +/- 0.2898 for 2014, from laser power measurements.  
//	atomic_filename(string("K_37_POPULATIONS_INPUT.txt"))
{
//	cout << "Called K37SublevelPopulations(" << the_sigma << ")." << endl;
	
	G4String configPath = CONFIGURATION_DIRECTORY;
	atomic_filename = configPath + "K_37_POPULATIONS_INPUT.txt";

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
	
	if(the_sigma==1)       { is_sigma_plus = true;  }
	else if(the_sigma==-1) { is_sigma_plus = false; }
	else                   { cout << "Bad sigma!!" << endl; assert(0); return; }
	
	// if it's sigma-, swap everything.
	if(!is_sigma_plus) { swap_states(); }
	
	
	// OK, so that's a perfectly good default way to start it.  
	// But actually, let's load population data from the text file instead. 
	// I'll leave the above thing alone because I don't know whether my std::vectors 
	// will be the right size otherwise.
//	cout << "We'll setup the pops from the inputs map..." << endl;
	
	theInputs = SS::loadup_textfile(atomic_filename);     // reads the atomic text file into 'theInputs'.
	Setup_Pops_From_InputsMap();
//	cout << "pops have been set up from the inputs map." << endl;
}

void K37SublevelPopulations::Setup_Pops_From_InputsMap( /*map<string, isotope_values * > theInputs*/)
{
// This will be a slow, computationally inefficient way to set up the populations.  
// It will also be pretty inelegant.  I don't care.
	
//	theInputs = SS::loadup_textfile(atomic_filename);     // reads the atomic text file into 'theInputs'.
	
	bool tmp_sigma = FindValue("IS_SIGMA_PLUS");
	if(tmp_sigma) { is_sigma_plus=true;  }
	else          { is_sigma_plus=false; }
	
//	set_pop(string level, int F, int M_F, double the_pop) // level=="ground" || level=="excited"
//	cout << "ok, we're going to go set a bunch of pops now." << endl;
	
	set_pop("ground", 2,  2, FindValue("GROUND_F2_Mplus2") );
	set_pop("ground", 2,  1, FindValue("GROUND_F2_Mplus1") );
	set_pop("ground", 2,  0, FindValue("GROUND_F2_M0") );
	set_pop("ground", 2, -1, FindValue("GROUND_F2_Mminus1") );
	set_pop("ground", 2, -2, FindValue("GROUND_F2_Mminus2") );
	
	set_pop("ground", 1,  1, FindValue("GROUND_F1_Mplus1") );
	set_pop("ground", 1,  0, FindValue("GROUND_F1_M0") );
	set_pop("ground", 1, -1, FindValue("GROUND_F1_Mminus1") );


	set_pop("excited", 2,  2, FindValue("EXCITED_F2_Mplus2") );
	set_pop("excited", 2,  1, FindValue("EXCITED_F2_Mplus1") );
	set_pop("excited", 2,  0, FindValue("EXCITED_F2_M0") );
	set_pop("excited", 2, -1, FindValue("EXCITED_F2_Mminus1") );
	set_pop("excited", 2, -2, FindValue("EXCITED_F2_Mminus2") );

	set_pop("excited", 1,  1, FindValue("EXCITED_F1_Mplus1") );
	set_pop("excited", 1,  0, FindValue("EXCITED_F1_M0") );
	set_pop("excited", 1, -1, FindValue("EXCITED_F1_Mminus1") );
	
	// ok, now what?  check if it's normalized?
//	cout << "We'll renormalize." << endl;
	renormalize();
	
	// check if it has the right sigma?
}

K37SublevelPopulations::~K37SublevelPopulations()
{
//	cout << "why has this even been called??" << endl;
	G4cout << "Deleting the sublevel populations.  Apparently." << G4endl;
}

// ------------------------------------------------------------- //
void K37SublevelPopulations::swap_states() // PRIVATE!  swap populations, but not the sigma flag.
{
	bool verbose = false;
	
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
	if(get_Mz()>=0) { physics_sigma_plus=true;  }
	else            { physics_sigma_plus=false; }
	
	if(is_sigma_plus && physics_sigma_plus)        { return  1.0; }
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
void K37SublevelPopulations::setup_sigma() // sets is_sigma_plus from polarization value.
{
	double pol = get_P();
	if(pol >= 0) { is_sigma_plus = true;  }
	else         { is_sigma_plus = false; }
}

void K37SublevelPopulations::set_sigma_plus() // calls functions to swap populations, *and* adjusts the is_sigma_plus flag.
{
	bool verbose = true;
	if(get_sigma()>0) // sigma +
	{ 
		if(verbose)
		{
			cout << "Called set_sigma_plus(), but there's nothing to be done." << endl;
		}
	} // don't do anything!
	else if(get_sigma()<0) // sigma -
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
		return;
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
//	cout << "Called get_Mz()." << endl;
	bool verbose=false;
//	bool verbose = true;
//	renormalize();
	
	double running_Mz = 0.0;
	double the_pop   = 0.0;
	double the_scale = 0.0;

	//
	int F=2;
	for(int M_F=-2; M_F<=2; M_F++)
	{
		the_pop     = get_pop("excited", F, M_F);
	//	cout << "get_Mz() hasn't died horribly at finding the_pop." << endl;
	//	cout << "the_pop = " << the_pop << endl;
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
//	bool verbose = false;
//	renormalize();
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
		//
	//	if(verbose)
	//	{
	//		G4cout << "Ground:   F=" << F << ", M_F=" << M_F << ":  pop=" << get_pop("ground", F, M_F)  << ", scale=" << get_scale("M_z2", F, M_F);
	//			G4cout << ";\tsubtotal=" << get_pop("ground", F, M_F)*get_scale("M_z2", F, M_F) << G4endl;
	//		G4cout << "Excited:  F=" << F << ", M_F=" << M_F << ":  pop=" << get_pop("excited", F, M_F) << ", scale=" << get_scale("M_z2", F, M_F);
	//			G4cout << ";\tsubtotal=" << get_pop("excited", F, M_F)*get_scale("M_z2", F, M_F) << G4endl;
	//		G4cout << "running_Mz2 = " << running_Mz2 << G4endl;
	//	}
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
		//
	//	if(verbose)
	//	{
	//		G4cout << "Ground:   F=" << F << ", M_F=" << M_F << ":  pop=" << get_pop("ground", F, M_F)  << ", scale=" << get_scale("M_z2", F, M_F);// << G4endl;
	//			G4cout << ";\tsubtotal=" << get_pop("ground", F, M_F)*get_scale("M_z2", F, M_F) << G4endl;
	//		G4cout << "Excited:  F=" << F << ", M_F=" << M_F << ":  pop=" << get_pop("excited", F, M_F) << ", scale=" << get_scale("M_z2", F, M_F);// << G4endl;
	//			G4cout << ";\tsubtotal=" << get_pop("excited", F, M_F)*get_scale("M_z2", F, M_F) << G4endl;
	//		G4cout << "running_Mz2 = " << running_Mz2 << G4endl;
	//	}
	}
	
	return running_Mz2;
}
double K37SublevelPopulations::get_Mz3()
{
//	renormalize();
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
	double Mz2 = get_Mz2();
//	G4cout << "Total Mz2 = " << Mz2 << G4endl;
	
	double T = 5.0/4.0 - Mz2;
	return T;
}

void K37SublevelPopulations::set_pop(string level, int F, int M_F, double the_pop) // level=="ground" || level=="excited"
{
	// sanity checks first.
	if( level.compare(string("ground"))!=0 && level.compare(string("excited"))!=0 )
	{
		cout << "Level:  " << level << " is not a thing.  Try again." << endl;
		assert(0);
		return;
	}
	if( !sanity(F, M_F) ) { assert(0); return; } 
	// If we got here, it passed all sanity checks.
	
	// now actually set it.
	// Now actually do the thing.
	int offset;
	if(F==1) { offset=1; }
	if(F==2) { offset=2; }
	
	if( !level.compare(string("ground")) ) // level is ground level
	{
		if(F==1)
		{
			ground_F1[M_F+offset] = the_pop;
		}
		else if(F==2)
		{
			ground_F2[M_F+offset] = the_pop;
		}
	}
	else if( !level.compare(string("excited")) ) // level is excited level
	{
		if(F==1)
		{
			excited_F1[M_F+offset] = the_pop;
		}
		else if(F==2)
		{
			excited_F2[M_F+offset] = the_pop;
		}
	}
	return;
}
double K37SublevelPopulations::get_pop(string level, int F, int M_F)
{
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

bool K37SublevelPopulations::sanity(int F, int M_F) // F==1 || F==2;  |M_F| <= F.
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
void K37SublevelPopulations::renormalize(bool verbose)
{
//	verbose = true;
	bool do_the_thing = false;
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
	if( abs(running_sum - 1.0) > allowed_mismatch ) 
	{
		do_the_thing = true; // it's false unless we do this here.
	}
	
	if(verbose)
	{
	//	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_P();
		cout << std::scientific << std::setprecision(4);
		cout << "-" << endl;
		cout << "renormalizeation:  " << endl;
		cout << "running_sum       = " << running_sum << endl;
		cout << "running_sum - 1.0 = " << running_sum - 1.0 << endl;
		cout << "allowed_mismatch  = " << allowed_mismatch << endl;
		cout << "do_the_thing = " << do_the_thing << endl;
		cout << "-" << endl;
	}
	
	// run through them all again.
	if(do_the_thing)
	{
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
	}
	if(verbose)
	{
		cout << "We're done renormalizing." << endl;
	}
	return;
}
void K37SublevelPopulations::print_pops()
{
	int the_width=12;
	int the_precision=the_width-4;
	
	cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
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
	cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
}
void K37SublevelPopulations::print_moments()
{
	int the_width=12;
	int the_precision=6;
	
//	int the_width_3 = 10;
	int the_width_2 = 6;
	int the_precision_2 = 4;
	
	cout << "-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --" << endl;
	cout << "P = ";
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_P();
	cout << endl;
	cout << "T = ";
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_T();	
	cout << endl;
	//
	cout << "-- " << endl;
	//
	cout << "Mz  = ";
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_Mz();	
	cout << "( max = " << std::fixed << std::setw(the_width_2) << std::setprecision(the_precision_2) << 3.0/2.0 << " )";
	cout << endl;
	
	cout << "Mz2 = ";
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_Mz2();	
	cout << "( max = " << std::fixed << std::setw(the_width_2) << std::setprecision(the_precision_2) << 9.0/4.0 << " )";
	cout << endl;
	
	cout << "Mz3 = ";
	cout << std::fixed << std::setw(the_width) << std::setprecision(the_precision) << get_Mz3();	
	cout << "( max = " << std::fixed << std::setw(the_width_2) << std::setprecision(the_precision_2) << 27.0/8.0 << " )";
	cout << endl;
	//
	cout << "-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --" << endl;
	return;
}
// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // 

void K37SublevelPopulations::AdjustPolarization(double new_Pol)
{
	int verbose=2;
	if(verbose)
	{
		cout << "Adjusting polarization from P_old=" << get_P() << " to P_new=" << new_Pol << endl;
	}
	double abs_old_Pol = abs(get_P());
	if( abs_old_Pol+allowed_mismatch >= 1.0) // we've lost sublevel information.  reload states from inputs map.
	{
		K37SublevelPopulations::Setup_Pops_From_InputsMap(); // also sets the sigma.
	//	if(verbose>1)
	//	{
	//		cout << "we've set up the populations from the saved inputs.  probably.  what even are the saved inputs?!" << endl;
	//		print_pops();
	//	}
		abs_old_Pol = abs(get_P());  // we set this at the beginning, but we're going to use it for calculations later, and we want the reset value.
	}
	if(abs(new_Pol) > 1.0)
	{
		cout << "Can't set polarization to " << new_Pol << endl;
		cout << "Polarization will remain unchanged." << endl;
		return;
	}
	
	// check normalization first?
	renormalize();
	
	
	// check that new_pol and old_pol even have the same sigma??
	bool new_is_sigma_plus;
	if(new_Pol >= 0) { new_is_sigma_plus = true;  }
	else             { new_is_sigma_plus = false; }
//	cout << "new_is_sigma_plus:  " << new_is_sigma_plus << endl;
//	cout << "old is_sigma_plus:  " << is_sigma_plus << endl;
	if(new_is_sigma_plus != is_sigma_plus) // if they're not the same sigma, we'll just swap the atomic states around now.
	{
	//	cout << "new sigma is different than (new) old sigma." << endl;
		swap_states();
		setup_sigma();  // sets is_sigma_plus to its new physically correct value.
	}
//	else if(new_is_sigma_plus == is_sigma_plus)
//	{
//		cout << "new sigma is same as old sigma." << endl;
//	}
//	else
//	{
//		cout << "sigma is broken." << endl;
//	}
	if(new_is_sigma_plus != is_sigma_plus)
	{
		cout << "Well fuck, you just broke mathematical logic." << endl;
		assert(0);
		return;
	}
	
	// OK, the sigma is the same now, in the correct "new" direction.  
	double P_stretched_old, P_other_old;
	
	double pop_stretched_old, pop_other_old;
	if(is_sigma_plus) { pop_stretched_old = get_pop("ground", 2, 2);  }
	else              { pop_stretched_old = get_pop("ground", 2, -2); }
	
	// for this math, we assume everything is sigma plus, then put the signs back in at the end.
	P_stretched_old = pop_stretched_old;         // sigma +
	P_other_old     = abs_old_Pol - P_stretched_old; // all positive because of how we've defined other things.  sigma+.
	//
	pop_other_old = 1.0 - pop_stretched_old;
	
	double pop_stretched_new, pop_other_new;
	pop_stretched_new = (pop_other_old*abs(new_Pol) - P_other_old) / (pop_other_old - P_other_old);
	pop_other_new = 1.0 - pop_stretched_new;
	if( pop_stretched_new < 0 || pop_other_new < 0)
	{
		cout << "pop_stretched_new = " << pop_stretched_new << endl;
		cout << "pop_other_new     = " << pop_other_new << endl;
		cout << "*** Cannot set polarization to this value using the assumptions we've made about the sublevel population distribution." << endl;
		cout << "*** Polarization remains unchanged." << endl;
		return;
	}
	
	double rescale_other = pop_other_new / pop_other_old;
	if(rescale_other < 0)  // allow it to be equal to zero, for P=+/-1.  We don't divide by this.
	{ 
		cout << "This is bad.  Very bad." << endl; 
		assert(0);
		return;
	}
	
//	bool verbose=false;
//	if(verbose)
//	{
//		cout << "abs_old_Pol:  " << abs_old_Pol << ";\tnew_Pol:  " << new_Pol << endl;
//		cout << "rescale_other   = " << rescale_other << endl;
//		cout << "1/rescale_other = " << 1.0/rescale_other << endl;
//		cout << "pop_stretched_old = " << pop_stretched_old << ";\tpop_other_old = " << pop_other_old << endl;
//		cout << "pop_stretched_old + pop_other_old = " << pop_stretched_old + pop_other_old << endl;
//		cout << "pop_stretched_new = " << pop_stretched_new << ";\tpop_other_new = " << pop_other_new << endl;
//		cout << "pop_stretched_new + pop_other_new = " << pop_stretched_new + pop_other_new << endl;
//	}
	
	
	// Now we rescale:
	int F, MF;
	double new_tmp_pop;
	
	F=2; 
	for(MF=2; MF>=-2; MF--)
	{
		new_tmp_pop = rescale_other*get_pop("excited", F, MF);
		set_pop("excited", F, MF, new_tmp_pop );
	}
	//
	F=2; MF=1;
	new_tmp_pop = rescale_other*get_pop("ground", F, MF);
	set_pop("ground", F, MF, new_tmp_pop );
	F=2; MF=0;
	new_tmp_pop = rescale_other*get_pop("ground", F, MF);
	set_pop("ground", F, MF, new_tmp_pop );
	F=2; MF=-1;
	new_tmp_pop = rescale_other*get_pop("ground", F, MF);
	set_pop("ground", F, MF, new_tmp_pop );
	
	if(is_sigma_plus)
	{
		F=2; MF=2;
	//	new_tmp_pop = rescale_other*get_pop("ground", F, MF);
		new_tmp_pop = pop_stretched_new;
		set_pop("ground", F, MF, new_tmp_pop );
		
		F=2; MF=-2;  // sigma+.  treat this state like the others.
		new_tmp_pop = rescale_other*get_pop("ground", F, MF);
		set_pop("ground", F, MF, new_tmp_pop );
	}
	else
	{
		F=2; MF=2;  // sigma-.  treat this state like the others.
		new_tmp_pop = rescale_other*get_pop("ground", F, MF);
		set_pop("ground", F, MF, new_tmp_pop );
		
		F=2; MF=-2;
	//	new_tmp_pop = rescale_other*get_pop("ground", F, MF);
		new_tmp_pop = pop_stretched_new;
		set_pop("ground", F, MF, new_tmp_pop );
	}
	//
	F=1; 
	for(MF=1; MF>=-1; MF--)
	{
		new_tmp_pop = rescale_other*get_pop("excited", F, MF);
		set_pop("excited", F, MF, new_tmp_pop );
	}
	F=1; 
	for(MF=1; MF>=-1; MF--)
	{
		new_tmp_pop = rescale_other*get_pop("ground", F, MF);
		set_pop("ground", F, MF, new_tmp_pop );
	}

//	if(verbose)
//	{
//		cout << "Populations are rescaled according to polarization, but haven't been normalized again yet..." << endl;
//		print_pops();
//		print_moments();
//		cout << "....Renormalizing now." << endl;
//	}
	
	// Still normalized?
	renormalize(false);

}


// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // 
// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // 
// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // 
// The stuff below is just all copied from HolsteinVars.  It's a very inelegant solution.  
// Originally, Spencer probably implemented this code in a much nicer way.  Sorry, Spencer.
// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // 
/*
bool K37SublevelPopulations::loadup_textfile(string paramfilename)
// puts the contents of the text file into 'theInputs'.
{
	bool verbose = false;
	string isotopeName;
	
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
*/
double K37SublevelPopulations::FindValue(const std::string &key_) const
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
double K37SublevelPopulations::FindUncertainty(const std::string &key_) const
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
void K37SublevelPopulations::print_isotope_values()
// for debugging.  Only for values loaded from the file.
{
	for (auto themap : theInputs) 
		{ themap.second -> Print(); }
	return;
}




/*
// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // 
// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // 
void K37SublevelPopulations::killall_pops() // when this gets done, there is no population. it's unphysical.  must be private.
{
	excited_F1.clear();
	excited_F2.clear();
	ground_F1.clear();
	ground_F2.clear();
	//
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
	ground_F2.push_back(0);
}
void K37SublevelPopulations::SetPops_Ns_NG1_NG2(double Ns, double NG1, double NG2, int sigma)
{
	double the_sum = Ns+NG1+NG2;
	if( abs(the_sum-1.0) > allowed_mismatch )
	{
		cout << "Ns+NG1+NG2 = " << the_sum << ";  Ns+NG1+NG2-1 = " << the_sum - 1.0 << ";\tRenormalizing now." << endl;
		Ns  = Ns /the_sum;
		NG1 = NG1/the_sum;
		NG2 = NG2/the_sum;
	}
	//
	killall_pops();
	if(sigma>0)
	{
		//                F, M_F,pop
		set_pop("ground", 2, 2,  Ns);
		set_pop("ground", 2, 1,  NG2);
		set_pop("ground", 1, 1,  NG1);
	}
	else
	{
		//                F,M_F,  pop
		set_pop("ground", 2, -2,  Ns);
		set_pop("ground", 2, -1,  NG2);
		set_pop("ground", 1, -1,  NG1);
	}
}
void K37SublevelPopulations::Setup_FromPolarizationOnly(double pol)
{
	bool verbose = true;
	if( pol>1.0 || pol <-1.0 )
	{
		cout << "Bad.  Can't have polarization of " << pol << endl;
		assert(0);
		return;
	}
	//
	int sigma = 0;
	if(pol>0) { sigma= 1;            }
	else      { sigma=-1; pol*=-1.0; }
	//
	//	double the_ratio = 1.7557;  // <NG2/NG2>, from laser power info with our toy model. 1.7557 +/- 0.2898.
	
	//
	double Ns  = ((1.0+op_power_ratio)*pol - (5.0/6.0 + op_power_ratio/2.0)) / (1.0/6.0+op_power_ratio/2.0);
	double NG1 = (pol - Ns) / ( 5.0/6.0 + op_power_ratio/2.0);
	double NG2 = op_power_ratio*NG1;
	
//	Don't bother with this, we do it in SetPops_Ns_NG1_NG2(...) anyway.
//	// Normalize 'em now, just in case?
//	double the_sum = Ns+NG1+NG2;
//	if(the_sum != 1.0)
//	if( abs(the_sum-1.0) > allowed_mismatch )
//	{
//		cout << "Ns+NG1+NG2 = " << the_sum << ";  Ns+NG1+NG2-1 = " << the_sum - 1.0 << ";\tRenormalizing now." << endl;
//		Ns  = Ns /the_sum;
//		NG1 = NG1/the_sum;
//		NG2 = NG2/the_sum;
//	}
	
	SetPops_Ns_NG1_NG2(Ns, NG1, NG2, sigma);
	
	if(verbose)
	{
		cout << "Input polarization was: " << pol*double(sigma) << endl;
		cout << "Ns  = " << Ns << endl;
		cout << "NG1 = " << NG1 << endl;
		cout << "NG2 = " << NG2 << endl;
		cout << "New polarization is:    " << get_P() << endl;
		cout << "New alignment is:       " << get_T() << endl;
		cout << "4P+T = " << 4.0*get_P() + get_T() << ";  (expect 3)" << endl;
		cout << "--" << endl;
	}
	//
	return;
}
void K37SublevelPopulations::Setup_FromPolarizationAlignment(double pol, double ali)
{
	bool verbose = true;
	if( pol>1.0 || pol<-1.0 )
	{
		cout << "Bad.  Can't have polarization of " << pol << endl;
		assert(0);
		return;
	}
	if( ali>1.0 || ali<-1.0 )
	{
		cout << "Bad.  Can't have alignment of " << ali << endl;
		assert(0);
		return;
	}
	//
	int sigma = 0;
	if(pol>0) { sigma= 1;            }
	else      { sigma=-1; pol*=-1.0; }
	//
	double NG1 = 18.0*pol - 10.0 + 8.0*ali;
	double NG2 = 16.0/3.0 - 6.0*pol - 8.0/3.0*ali;
	double Ns  = 17.0/3.0 - 12.0*pol - 16.0/3.0*ali;
	
	SetPops_Ns_NG1_NG2(Ns, NG1, NG2, sigma);
	
	if(verbose)
	{
		cout << "Input polarization was: " << pol*double(sigma) << endl;
		cout << "Ns  = " << Ns << endl;
		cout << "NG1 = " << NG1 << endl;
		cout << "NG2 = " << NG2 << endl;
		cout << "New polarization is:    " << get_P() << endl;
		cout << "New alignment is:       " << get_T() << endl;
		cout << "--" << endl;
	}
	//

	//
	return;
}
// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // 
*/

/*
// ------------------------------------------------------------- //
void K37SublevelPopulations::Setup_FromPolarizationOnly(double pol)
{
	SetPolarizationOnly(pol);
	MakeAlignmentFromPolarization();
	MakeOctopoleFromPolarizationAlignment();
	return;
}
void K37SublevelPopulations::Setup_FromPolarizationAlignment(double pol, double ali)
{
	SetPolarizationOnly(pol);
	SetAlignmentOnly(ali);
	MakeOctopoleFromPolarizationAlignment();
	return;
}
void K37SublevelPopulations::Setup_FromPolarizationAlignmentOctopole(double pol, double ali, double oct)  // this method is stupid in that it mixes conventions.  Probably just avoid using it.
{
	SetPolarizationOnly(pol);
	SetAlignmentOnly(ali);
	SetOctopoleOnly(oct);
	return;
}
void K37SublevelPopulations::Setup_FromDipoleQuadrupoleOctopole(double dip, double quad, double oct)
{
//	G4cout << "WARNING:  It is not yet supported to set up multipole moments directly." << G4endl;
//	Set_Mz_Only();
//	Set_Mz2_Only();
	Set_Mz3_Only(oct);	
	return;
}


void K37SublevelPopulations::MakeAlignmentFromPolarization()  // doesn't work yet!!
{
	G4cout << "*** WARNING:  Setting up the alignment parameter from only polarization information is not yet supported." << G4endl;
	// T = 2.0-3.0*P ;
	SetAlignmentOnly(1.0);
	return;
}
void K37SublevelPopulations::MakeOctopoleFromPolarizationAlignment()
{
	G4cout << "*** WARNING:  Setting up the octopole moment from only polarization and alignment information is not yet supported." << G4endl;
//	SetOctopoleOnly(1.0);  what even *is* a "default" octopole moment??
	return;
}


// ------------------------------------------------------------- //
void K37SublevelPopulations::SetPolarizationOnly(double pol) // Must be private!
{
	G4cout << "*** WARNING:  Setting up polarization directly is not yet supported." << G4endl;
//	
//	if (fabs(pol) <= 1.0) 
//	{
//		polarization_ = pol;
//		G4cout << "Polarization set to " << pol << G4endl;
//	} 
//	else 
//	{
//		G4cout << "WARNING: Polarization " << pol
//			   << " not in allowed range.  No changes made." << G4endl;
//	}
//	
}
void K37SublevelPopulations::SetAlignmentOnly(double ali) // Must be private!
{
	G4cout << "*** WARNING:  Setting up alignment directly is not yet supported." << G4endl;
	
//	if (ali <= 1.0 && ali >= 0.0) 
//	{
//		alignment_ = ali;
//		G4cout << "Alignment set to " << ali << G4endl;
//	} 
//	else if (ali < 0.0 && ali >= -1.0) 
//	{
//		ali = fabs(ali);
//		G4cout << "WARNING: Alignment must be positive (for both polarization states)";
//		G4cout << G4endl;
//		G4cout << "Note that this sign convention is opposite that in JTW 1957 and that";
//		G4cout << G4endl;
//		G4cout << "this change is fully documented in the \"JTWEvent.cc\" source file";
//		G4cout << G4endl;
//		G4cout << "*********************************************************" << G4endl;
//		G4cout << "   Setting alignment to T = +" << ali << G4endl;
//		G4cout << "*********************************************************" << G4endl;
//		alignment_ = ali;
//	} 
//	else 
//	{
//		G4cout << "WARNING: Alignment " << ali
//			   << " not in allowed range.  No changes made." << G4endl;
//	}
}
void K37SublevelPopulations::SetOctopoleOnly(double oct)  // Must be private!
{
	G4cout << "*** WARNING:  Octopole moment is not yet supported." << G4endl;
	return;
}

void K37SublevelPopulations::Set_Mz_Only(double new_Mz)
{
	G4cout << "*** WARNING:  Setting Mz directly is not yet supported." << G4endl;
	return;
}
void K37SublevelPopulations::Set_Mz2_Only(double new_Mz2)
{
	G4cout << "*** WARNING:  Setting Mz2 directly is not yet supported." << G4endl;
	return;
}
void K37SublevelPopulations::Set_Mz3_Only(double new_Mz3)
{
	G4cout << "*** WARNING:  Setting Mz3 directly is not yet supported." << G4endl;
	return;
}
*/
