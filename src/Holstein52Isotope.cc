// Author: Melissa Anholm - 2019

#include <string>  // for SS
#include <sstream> // for SS
#include <vector>  // for SS

#undef NDEBUG
#include<assert.h>

#include "Holstein52Isotope.hh"  // formerly HolsteinVars
#include "K37SublevelPopulations.hh"
#include "IsotopeValues.hh"
#include "K37Config.hh"
//#include "SplitString.hh"
#include "globals.hh"


using std::cout;
using std::endl;
using std::string;
using std::vector;

//const double pi = std::atan(1.0)*4.0;

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
// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- //

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
	sigma_M1(0), sigma_M2(0)
//	nuclear_filename(string("K_37_INPUT.txt"))
//	atomic_filename(string("K_37_POPULATIONS_INPUT.txt"))
{
	G4cout << "Creating a new HolsteinVars()." << G4endl;
	G4String configPath = CONFIGURATION_DIRECTORY;
	G4String nuclear_filename = configPath + "K_37_INPUT.txt";
	theInputs = SS::loadup_textfile(nuclear_filename);    // reads the nuclear text file into 'theInputs'.
//	loadup_textfile(atomic_filename);     // reads the atomic text file into 'theInputs'.
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
//	cout << std::setprecision(8);
	cout << std::setprecision(50);
	cout << "(E0/MeV)  = " << (E0/MeV) << endl;
	cout << "(m_e/MeV) = " << (m_e/MeV) <<  endl;
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

/*
//bool HolsteinVars::loadup_textfile(string paramfilename)
map<string, isotope_values * > HolsteinVars::loadup_textfile(string paramfilename)
map<string, isotope_values * > SS::loadup_textfile(string paramfilename)
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
	//return true;
	return theInputs;
}
*/

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


