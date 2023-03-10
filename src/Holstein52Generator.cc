// Author: Melissa Anholm - 2019

#include <cmath>
#include <iostream>  // cout, endl

#undef NDEBUG
#include<assert.h>

#include<Randomize.hh>  // also includes:  CLHEP/Random/RandFlat.h, and #define's G4RandFlat   CLHEP::RandFlat.


#include "Holstein52Generator.hh"
#include "Holstein52Isotope.hh"  // formerly HolsteinVars
//	//#include "IsotopeValues.hh"
//	//#include "SplitString.hh"
#include "K37SublevelPopulations.hh"
#include "K37FermiFunction.hh"

using std::cout;
using std::endl;
using std::abs;
using std::sqrt;  // do I need this??  who knows, but I can't take the chance.

//  Holstein52Generator class relies on other classes: 
//  	HolsteinVars
//  		isotope_values
//  		SS (splitstring)
//  	K37SublevelPopulations


// I think this is the equivalent class of K37EventGenerator.



// ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- // ---- //
void Holstein52Generator::SetAcceptanceMode(string the_mode)
{
	cout << "SetAcceptanceMode(...) doesn't work yet." << endl;
	return;
}

/*
Holstein52Generator::Holstein52Generator() // broken?
{
	cout << "Initializing Holstein52Generator without input parameters.  " << endl;
	cout << "Creating a new set of HolsteinVars, and a new set of K37SublevelPopulations.  " << endl;
	cout << "I think this is broken and will cause segfaults when we try to look for the HolsteinVars or the K37SublevelPopulations." << endl;
	cout << "Hard kill it now." << endl;
	assert(0);
	
	//
//	
//	HolsteinVars * HV           = new HolsteinVars();
//	K37SublevelPopulations * pops = new K37SublevelPopulations(1);
//	cout << "things are created." << endl;
//	Holstein52Generator( (HolsteinVars*)HV, (K37SublevelPopulations*)pops );
//	
	Params   = new HolsteinVars();
	the_pops = new K37SublevelPopulations(1);
	Holstein52Generator( (HolsteinVars*)Params, (K37SublevelPopulations*)the_pops );
}
*/

//Holstein52Generator::Holstein52Generator(HolsteinVars * HV, K37AtomicSetup * atomic_setup):
//Holstein52Generator::Holstein52Generator(HolsteinVars * HV, K37AtomicSetup * atomic_setup, K37FermiFunction * FF):
Holstein52Generator::Holstein52Generator(HolsteinVars * HV, K37AtomicSetup * atomic_setup):
	// nuclear parameters:
	delta_M1(0), delta_M2(0), delta_deltaC(0), 
	delta_mu_parent(0), delta_mu_daughter(0), delta_quad_parent(0), delta_quad_daughter(0), 
	// run parameters:
	use_roc(true),
	use_cone(false), cone_costheta(0),
	prob_max(0.2), 
	do_coulomb(true)
{
	// Nuclear:
	Params = HV;
	
	M_F_isospin         = Params->get_M_F_isospin();    // exact
	I_spin              = Params->get_I_spin();         // exact
	A_nucleons          = Params->get_A_nucleons();     // exact
	hbarc_eV_nm         = Params->get_hbarc_eV_nm();    // exact.  unitless. ... it is *not* fucking unitless.  ... it needs to be unitless for now though, because who fucking knows what units the matrix elements are supposed to have.
	Z_daughter          = Params->get_Z_daughter();     // exact
	speed_of_light      = Params->get_speed_of_light(); // exact.  double.  no G4units.
//	alpha_finestructure = Params->get_finestructure();
	alpha               = Params->get_alpha();
	
	E0    = Params->get_E0();    // uncertain (propagation) MeV.  comes out in MeV from HolsteinVars.
	M     = Params->get_M();     // uncertain (propagation) MeV.  comes out in MeV from HolsteinVars.
	Delta = Params->get_Delta(); // uncertain (propagation) MeV.  comes out in MeV from HolsteinVars.  units propagated.
	m_e   = Params->get_m_e();   // uncertain (direct).     MeV.  comes out in MeV from HolsteinVars.
	
	deltaC = Params->get_deltaC(); // uncertain (direct)
	
	quad_parent   = Params->get_quad_parent();   // uncertain (direct) 
	quad_daughter = Params->get_quad_daughter(); // uncertain (direct) 
	mu_parent     = Params->get_mu_parent();     // uncertain (direct) 
	mu_daughter   = Params->get_mu_daughter();   // uncertain (direct) 
	
	// dont' set up g_V and g_A, because I'll want to be able to update their values from within G4..
//	g_V  = Params->get_g_Vector();  // exact
//	g_A  = Params->get_g_Axial();  // exact-ish??
	
	g_II = Params->get_g_II(); // exact
	g_S  = Params->get_g_S();  // exact
	g_P  = Params->get_g_P();  // exact
	g_M  = Params->get_g_M();  // exact
	
	M_F     = Params->get_M_F();     // uncertain (propagation, deltaC) 
	M_GT    = Params->get_M_GT();    // uncertain (propagation, deltaC) 
	M_r2    = Params->get_M_r2();    // exact
	M_sr2   = Params->get_M_sr2();   // exact
	M_sL    = Params->get_M_sL();    // exact
	M_srp   = Params->get_M_srp();   // exact
	M_rdotp = Params->get_M_rdotp(); // exact
	M_1y    = Params->get_M_1y();    // exact
	M_2y    = Params->get_M_2y();    // exact
	M_3y    = Params->get_M_3y();    // exact
	M_Q     = Params->get_M_Q();     // uncertain (propagation) -- I think this doesn't get used directly...
	M_rp    = Params->get_M_rp();    // uncertain (propagation)
	// M_L = Params->M_L;  // exact
	
	a1 = Params->get_a1(); // uncertain (propagation, M1 M2, M_F)          // note:  a1 gets re-calculated in randomize_nuclear, so don't have to worry about propagating changes to g_V right now.  This is just a placeholder value.  
	a2 = Params->get_a2(); // uncertain (propagation, M1 M2)
	c1 = Params->get_c1(); // uncertain (propagation, M1 M2, M_GT, (g_A) )
	c2 = Params->get_c2(); // uncertain (propagation, M1 M2, (g_A) )
	
	b  = Params->get_b();  // uncertain (propagation, mu_parent mu_daughter)
	d  = Params->get_d();  // uncertain (propagation, M1 M2, M_GT, (g_A) )
	e  = Params->get_e();  // exact by isospin(?)+CVC.
	g  = Params->get_g();  // uncertain (propagation, M1 M2, quad_parent, quad_daughter)
	f  = Params->get_f();  // uncertain (propagation, M1 M2, quad_parent, quad_daughter)
	h  = Params->get_h();  // uncertain (propagation, M1 M2, M_GT, (g_A) )
	
	j2 = Params->get_j2(); // uncertain (propagation, M1 M2, (g_A) )
	j3 = Params->get_j3(); // uncertain (propagation, M1 M2, (g_A) )
	
	
	jtw_a1 = (Params->get_g_Vector())*Params->FindValue("M_F");  // hasn't had isospin correction turned on yet.
	jtw_c1 = (Params->get_g_Axial())*Params->FindValue("M_GT"); // hasn't had isospin correction turned on yet.
	jtw_xi    = 2.0*(jtw_a1*jtw_a1 + jtw_c1*jtw_c1);
	jtw_Abeta = 2.0*(jtw_c1*jtw_c1*(1.0/(3.0/2.0 + 1.0)) + 2.0*jtw_a1*jtw_c1*sqrt((3.0/2.0)/(3.0/2.0 + 1.0)) ) / (2.0*(jtw_a1*jtw_a1 + jtw_c1*jtw_c1));
	jtw_rho   = jtw_c1/jtw_a1;
	
	// b_kludge, but I'll need this later anyway:
//	alpha_finestructure = 
	jtw_gamma = sqrt(1.0-alpha*alpha* Z_daughter*Z_daughter);
	
	// Fermi Function:
	the_FF = new K37FermiFunction();
	X_coulomb = HV->get_X(); 
	Y_coulomb = HV->get_Y();
//	alpha     = HV->get_alpha();

	
	// Atomic:
	the_atomic_setup = atomic_setup;
	
	Mz  = the_atomic_setup->GetPops()->get_Mz();   // 1.5
	Mz2 = the_atomic_setup->GetPops()->get_Mz2();  // 2.25
	Mz3 = the_atomic_setup->GetPops()->get_Mz3();  // 3.375
	initialize_lambdafuncs();  // must have Mz, Mz2, Mz3 set up already.
	
	// Decay-specific Parameters:
	initial_momentum = G4ThreeVector(0.0*MeV,0.0*MeV,0.0*MeV);  // 
	
	//
//	cout << "in Holstein52Generator::Holstein52Generator(...), we have g_A = " << (Params->get_g_Axial()) << endl;
}


void Holstein52Generator::randomize_nuclear(bool doit)  // needs Params.
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
		delta_M1 = ( G4RandGauss::shoot(0, Params->get_sigma_M1()/MeV) )*MeV;  // MeV
		delta_M2 = ( G4RandGauss::shoot(0, Params->get_sigma_M2()/MeV) )*MeV;  // MeV
	}
	
	// masses and energies.  These are all in MeV.
	G4double delta_M = (delta_M1 + delta_M2)/2.0;  // both things are already in G4units of MeV.
	G4double delta_Delta = (delta_M1 - delta_M2);  // both things are already in G4units of MeV.
	
	M     = Params->get_M()     + delta_M;      // both things in G4units of MeV already.
	Delta = Params->get_Delta() + delta_Delta;  // both things in G4units of MeV already.
	m_e   = Params->get_m_e();     // don't even bother with this uncertainty right now.  comes out in G4units from HolsteinVars.
	E0    = (Delta/MeV)*(1.0+(m_e/MeV)*(m_e/MeV)/(2.0*(M/MeV)*(Delta/MeV)) )/(1.0+(Delta/MeV)/(2.0*(M/MeV))) * MeV;  // error is propagated already;  M1, M2 --> M, Delta --> E0.
	rc    = 1.0 + (Delta/MeV)/(2.0*(M/MeV));  // unitless.
	
//	// delta_E0 isn't really needed.  It's more vestigial.
//	double delta_E0 = (0.5 + ( pow(M-0.5*Delta, 2)-m_e*m_e )/( 2.0*pow(M+0.5*Delta,2) ))*delta_M1 + (-M+0.5*Delta)/(M+0.5*Delta)*delta_M2;
	
	// isospin correction (goes into M_F and M_GT):
	double delta_deltaC = 0;
	if(doit)
	{
		delta_deltaC = G4RandGauss::shoot(0, Params->get_sigma_deltaC() );
	}
	deltaC = Params->get_deltaC() + delta_deltaC;
	M_F  = Params->get_M_F()  * sqrt(1.0-Params->get_deltaC() ) / sqrt(1.0-deltaC);  // uncorrect the prev. correction, they re-correct with error.
	M_GT = Params->get_M_GT() * sqrt(1.0-Params->get_deltaC() ) / sqrt(1.0-deltaC);  // this probably doesn't capture all the physics of M_GT and deltaC ...
	/// !!!  MELISSA LOOK HERE!  
	/// !!!  I THINK THESE VALUES FOR M_F AND M_GT HAVE THE WRONG ISOSPIN SYMMETRY ERROR PROPAGATION, 
	/// !!!  AND NOT JUST BECAUSE OF SUBTLETIES!
	
	
	// magnetic moments:
	if(doit)
	{
		delta_mu_parent   = G4RandGauss::shoot(0, Params->get_sigma_mu_parent());
		delta_mu_daughter = G4RandGauss::shoot(0, Params->get_sigma_mu_daughter());
	}
	mu_parent         = Params->get_mu_parent()   + delta_mu_parent;   
	mu_daughter       = Params->get_mu_daughter() + delta_mu_daughter;
	
	// quadrupole moments:
	if(doit)
	{
		delta_quad_parent   = G4RandGauss::shoot(0, Params->get_sigma_quad_parent());    // MeV
		delta_quad_daughter = G4RandGauss::shoot(0, Params->get_sigma_quad_daughter());  // MeV
	}
	quad_parent   = Params->get_quad_parent()   + delta_quad_parent;
	quad_daughter = Params->get_quad_daughter() + delta_quad_daughter;
	
	// Use experimental quadrupole moments for 'g'-->M_Q, rather than calculated M_Q.
	double jterm = sqrt( (I_spin + 1.0)*(2.0*I_spin + 3.0)/( I_spin*(2.0*I_spin - 1.0) ));
	M_Q = (2.0/3.0) * M_F_isospin * jterm* (quad_daughter-quad_parent) / (hbarc_eV_nm*hbarc_eV_nm) / ( -4.0/3.0 * (Params->get_g_Vector()) );
	// Set M_rp by M_Q:
	M_rp = (E0/MeV)*M_Q/sqrt(6);
	
	// Calculate everything else:
	a1  = (Params->get_g_Vector())*(M_F - (Delta/MeV)*(Delta/MeV)/6.0*M_r2 + (Delta/MeV)/3.0*M_rdotp) / rc;
	a2  = (Params->get_g_Vector())*((M/MeV)*(M/MeV))*(M_r2/6.0) / rc;  
	c1  = (Params->get_g_Axial())*(M_GT - (Delta/MeV)*(Delta/MeV)/6.0*M_sr2 + M_1y*2.0*(Delta/MeV)*(Delta/MeV)/(6.0*sqrt(10.)) + A_nucleons*(Delta/MeV)/(2.0*(M/MeV))*M_sL + (Delta/MeV)/2.0*M_srp ) / rc;
	c2  = (Params->get_g_Axial())*((M/MeV)*(M/MeV))*( M_sr2/6.0 + M_1y/(6.0*sqrt(10.)) ) / rc;
	
	b   = A_nucleons*M_F_isospin * sqrt((I_spin+1.0)/I_spin) * (mu_parent-mu_daughter);
	
	double d_I  = (Params->get_g_Axial())*(-1.0*M_GT + (Delta/MeV)*(Delta/MeV)*M_sr2/6.0 + M_1y*(Delta/MeV)*((M/MeV)+(Delta/MeV)/6.0)/sqrt(10.) + (A_nucleons*M_sL) + ((M/MeV)*M_srp) ) / rc;
	double d_II = -1.0*g_II*A_nucleons*M_GT;
	d   = d_I + d_II;
	
	e   = 0; // Eq. (22).  assumes CVC.
	g   = (2.0/3.0) * (M/MeV)*(M/MeV) * M_F_isospin * jterm*(quad_daughter-quad_parent) / (hbarc_eV_nm*hbarc_eV_nm);
	f   = 2.0*(Params->get_g_Vector())*(M/MeV)*M_rp; // for f, use new value of M_rp
	h   = -1.0*((Params->get_g_Axial())*2.0*((M/MeV)*(M/MeV))*M_1y/sqrt(10.) + g_P*(A_nucleons*A_nucleons)*M_GT) / rc;
	
	j2  = -(Params->get_g_Axial())*2.0*((M/MeV)*(M/MeV))*M_2y/3.0;
	j3  = -(Params->get_g_Axial())*2.0*((M/MeV)*(M/MeV))*M_3y/3.0;
	
	// must still update the jtw parameters, because I might have changed g_A or M_GT.
	jtw_a1 = (Params->get_g_Vector())*M_F *sqrt(1.0-deltaC);  // turn isospin correction back off for this.  it's already included in M_F.
	jtw_c1 = (Params->get_g_Axial())*M_GT*sqrt(1.0-deltaC);  // turn isospin correction back off for this.  it's already included in M_GT.
	
	jtw_xi    = 2.0*(jtw_a1*jtw_a1 + jtw_c1*jtw_c1);
	jtw_Abeta = 2.0*(jtw_c1*jtw_c1*(1.0/(3.0/2.0 + 1.0)) + 2.0*jtw_a1*jtw_c1*sqrt((3.0/2.0)/(3.0/2.0 + 1.0)) ) / (2.0*(jtw_a1*jtw_a1 + jtw_c1*jtw_c1));
	jtw_rho   = jtw_c1/jtw_a1;
	
//	cout << "in Holstein52Generator::randomize_nuclear(...), we have g_A = " << (Params->get_g_Axial()) << endl;

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
		
	//	jtw_a1 = (Params->get_g_Vector())*M_F *sqrt(1.0-deltaC);  // turn isospin correction back off for this.  it's already included in M_F.
	//	jtw_c1 = (Params->get_g_Axial())*M_GT*sqrt(1.0-deltaC);  // turn isospin correction back off for this.  it's already included in M_GT.
		
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
	/*
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
	*/
	return;
}

void Holstein52Generator::initialize_lambdafuncs()
{
	get_lambda0();  // 1.0
	get_lambda1();
	get_lambda2();
	get_lambda3();
	get_lambda4();  // 0.0
}
double Holstein52Generator::get_lambda0()
{
	Lambda0=1; 
	return Lambda0;
}
double Holstein52Generator::get_lambda1()
{
	Lambda1 = Mz/I_spin; // Polarization.
	return Lambda1;
}
double Holstein52Generator::get_lambda2()
{
	double T = 5.0/4.0 - Mz2;
	Lambda2 = T * (2.0*I_spin - 1.0)/(I_spin+1.0);
	return Lambda2;
}
double Holstein52Generator::get_lambda3()
{
	double numerator =   5.0*Mz3;
	double denominator = I_spin*(3.0*I_spin*I_spin + 3.0*I_spin - 1.0);
	Lambda3 = Mz/I_spin - numerator/denominator;

	return Lambda3;
}
double Holstein52Generator::get_lambda4()
{
	Lambda4=0; 
	return Lambda4; 
}

double Holstein52Generator::F_0(double E) // a1 a2 c1 c2 (M/MeV) (E0/MeV) (m_e/MeV) b d h
{
	bool verbose=false;
	double t1, t2, t3, t4, t5;
	double bsm_term = 0;
	t1 = a1*a1 + 2.0*a1*a2/(3.0*(M/MeV)*(M/MeV))*(     (m_e/MeV)*(m_e/MeV) +  4.0*E*(E0/MeV) + 2.0*(m_e/MeV)*(m_e/MeV)*(E0/MeV)/E -  4.0*E*E);
	t2 = c1*c1 + 2.0*c1*c2/(9.0*(M/MeV)*(M/MeV))*(11.0*(m_e/MeV)*(m_e/MeV) + 20.0*E*(E0/MeV) - 2.0*(m_e/MeV)*(m_e/MeV)*(E0/MeV)/E - 20.0*E*E);
	t3 = -2.0*(E0/MeV)/(3.0*(M/MeV))*(c1*c1 + c1*d - c1*b);
	t4 =  2.0*E /(3.0*(M/MeV))*(3.0*a1*a1 + 5.0*c1*c1 -2.0*c1*b);
	t5 = -1.0*(m_e/MeV)*(m_e/MeV)/(3.0*(M/MeV)*E)*(-3.0*a1*e + 2.0*c1*c1 + c1*d - 2.0*c1*b - c1*h*((E0/MeV)-E)/(2.0*(M/MeV)));
	
	double g_S = Params -> get_g_Scalar();
	double g_T = Params -> get_g_Tensor();
	
	bsm_term = M_F*M_F*g_S*g_S + M_GT*M_GT*g_T*g_T;
	
	double F0 = t1 + t2 + t3 + t4 + t5 + bsm_term;
	if(verbose) 
	{ 
		cout << "F0 = " << F0 << endl;
		cout << "a1*a1 + c1*c1 = " << a1*a1 + c1*c1 << endl;   // JTW xi*factor.
		cout << "F_0(E) bsm_term = " << bsm_term << endl;
	}
	return F0;
}
double Holstein52Generator::F_1(double E) // a1 c1 c2 a2 (E0/MeV) (m_e/MeV) (M/MeV) b d j2 g f 
{
	bool verbose=false;
	
	double f1, f2, f3;
	double t1a, t1b, t2a, t2b, t3a, t3b;
	double bsm_term = 0;
	//
	double u          = Params->get_u(); 
	double delta_uv   = Params->get_delta_uv();
	double gamma_uv   = Params->get_gamma_uv();
	double lambda_uv  = Params->get_lambda_uv(); 
	
	/*
	f1 = Params->delta_uv*sqrt(Params->u/(Params->u+1)) * 2.0;
	t1a = a1*c1 - ((E0/MeV)/(3.0*(M/MeV)))*(a1*c1 + a1*d - a1*b) + (E/(3.0*(M/MeV)))*(7.0*a1*c1 - a1*b + a1*d);  
	t1b = (a1*c2 + c1*a2)*(4.0*E*((E0/MeV)-E)+3.0*(m_e/MeV)*(m_e/MeV))/(3.0*(M/MeV)*(M/MeV));
	
	f2 = Params->gamma_uv/(Params->u+1.0);
	t2a = c1*c1 + 2.0*c1*c2*(8.0*E*((E0/MeV)-E) + 3.0*(m_e/MeV)*(m_e/MeV))/(3.0*(M/MeV)*(M/MeV));
	t2b = -2.0*(E0/MeV)/(3.0*(M/MeV))*(c1*c1 + c1*d - c1*b) + E/(3.0*(M/MeV))*(11.0*c1*c1 - c1*d - 5.0*c1*b);
	
	f3 = Params->lambda_uv/(Params->u+1.0);
	t3a = -5.0*E/(M/MeV)*c1*f + sqrt(3.0/2.0)*c1*g*((E0/MeV)*(E0/MeV) - 11.0*(E0/MeV)*E + 6.0*(m_e/MeV)*(m_e/MeV) + 4.0*E*E)/(6.0*(M/MeV)*(M/MeV));
	t3b = -3.0*c1*j2*(8.0*E*E - 5.0*(E0/MeV)*E - 3.0*(m_e/MeV)*(m_e/MeV))/(6.0*(M/MeV)*(M/MeV));
	*/
	
	f1 = delta_uv*sqrt( u/(u+1)) * 2.0;
	t1a = a1*c1 - ((E0/MeV)/(3.0*(M/MeV)))*(a1*c1 + a1*d - a1*b) + (E/(3.0*(M/MeV)))*(7.0*a1*c1 - a1*b + a1*d);  
	t1b = (a1*c2 + c1*a2)*(4.0*E*((E0/MeV)-E)+3.0*(m_e/MeV)*(m_e/MeV))/(3.0*(M/MeV)*(M/MeV));
	
	f2 = gamma_uv/(u+1.0);
	t2a = c1*c1 + 2.0*c1*c2*(8.0*E*((E0/MeV)-E) + 3.0*(m_e/MeV)*(m_e/MeV))/(3.0*(M/MeV)*(M/MeV));
	t2b = -2.0*(E0/MeV)/(3.0*(M/MeV))*(c1*c1 + c1*d - c1*b) + E/(3.0*(M/MeV))*(11.0*c1*c1 - c1*d - 5.0*c1*b);
	
	f3 = lambda_uv/(u+1.0);
	t3a = -5.0*E/(M/MeV)*c1*f + sqrt(3.0/2.0)*c1*g*((E0/MeV)*(E0/MeV) - 11.0*(E0/MeV)*E + 6.0*(m_e/MeV)*(m_e/MeV) + 4.0*E*E)/(6.0*(M/MeV)*(M/MeV));
	t3b = -3.0*c1*j2*(8.0*E*E - 5.0*(E0/MeV)*E - 3.0*(m_e/MeV)*(m_e/MeV))/(6.0*(M/MeV)*(M/MeV));
	
	double g_S = Params -> get_g_Scalar();
	double g_T = Params -> get_g_Tensor();
	bsm_term = -1.0*f1*M_F*M_GT*g_S*g_T + f2*M_GT*M_GT*g_T*g_T;
	
	double F1 = f1*(t1a + t1b) + f2*(t2a + t2b) + f3*(t3a + t3b) + bsm_term;
	if(verbose) 
	{ 
		cout << "F1 = " << F1 << endl;
		cout << "c1*c1*(1.0/(3.0/2.0 + 1.0)) + 2*a1*c1*sqrt( (3.0/2.0) /(3.0/2.0 + 1.0)) = " << c1*c1*(1.0/(3.0/2.0 + 1.0)) + 2.0*a1*c1*sqrt( (3.0/2.0) /(3.0/2.0 + 1.0)) << endl; // JTW xi*factor*A_beta.
		cout << "F_1(E) bsm_term = " << bsm_term << endl;
	} 
	return F1;
}
double Holstein52Generator::F_2(double E) // c1 (E0/MeV) (M/MeV) a1 d b g j2 f j3 c2 (no dep. on a2 ??)
{
	bool verbose=false;
	
	// did I double-check this?
	double f2;
	double t1, t2, t3, t4;
	//
	double u          = Params->get_u(); 
	double theta_uv   = Params->get_theta_uv(); 
	double delta_uv   = Params->get_delta_uv();
	double kappa_uv   = Params->get_kappa_uv(); 
	double epsilon_uv = Params->get_epsilon_uv();

	t1 = theta_uv*E/(2.0*(M/MeV))*(c1*c1 + 8.0*c1*c2*((E0/MeV)-E)/(3.0*(M/MeV)) - c1*d - c1*b);
	f2 = -1.0*delta_uv*E/(M/MeV)*sqrt( (u*(u+1.0))/((2.0*u-1.0)*(2.0*u+3.0)) );
	t2 = sqrt(3.0/2.0)*a1*f + a1*g*(E+2.0*(E0/MeV))/(4.0*(M/MeV)) - sqrt(3.0/2.0)*a1*j2*((E0/MeV)-E)/(2.0*(M/MeV));
	t3 = kappa_uv*E/(2.0*(M/MeV))*(-3.0*c1*f -sqrt(3.0/2.0)*c1*g*((E0/MeV)-E)/(M/MeV) + 3.0*c1*j2*((E0/MeV)-2.0*E)/(2.0*(M/MeV)));
	t4 = epsilon_uv*c1*j3*(21.0*E*E)/(8.0*(M/MeV)*(M/MeV));

	/*
	t1 = Params->theta_uv*E/(2.0*(M/MeV))*(c1*c1 + 8.0*c1*c2*((E0/MeV)-E)/(3.0*(M/MeV)) - c1*d - c1*b);
	f2 = -1.0*Params->delta_uv*E/(M/MeV)*sqrt( (Params->u*(Params->u+1.0))/((2.0*Params->u-1.0)*(2.0*Params->u+3.0)) );
	t2 = sqrt(3.0/2.0)*a1*f + a1*g*(E+2.0*(E0/MeV))/(4.0*(M/MeV)) - sqrt(3.0/2.0)*a1*j2*((E0/MeV)-E)/(2.0*(M/MeV));
	t3 = Params->kappa_uv*E/(2.0*(M/MeV))*(-3.0*c1*f -sqrt(3.0/2.0)*c1*g*((E0/MeV)-E)/(M/MeV) + 3.0*c1*j2*((E0/MeV)-2.0*E)/(2.0*(M/MeV)));
	t4 = Params->epsilon_uv*c1*j3*(21.0*E*E)/(8.0*(M/MeV)*(M/MeV));
	*/
	
	double F2 = t1 + f2*t2 + t3 + t4;
	if(verbose) { cout << "F2 = " << F2 << endl; }
	return F2;
}
double Holstein52Generator::F_3(double E) // a1 j3 (M/MeV) c1 g j2 j3
{
	bool verbose=false;
	
	// did I double-check this?
	double f1;
	double t1, t2, t3;
	
	double u        = Params->get_u(); 
	double delta_uv = Params->get_delta_uv();
	double rho_uv   = Params->get_rho_uv();
	double sigma_uv = Params->get_sigma_uv(); 
	
	f1 = -1.0*delta_uv*(3.0*u*u + 3.0*u-1.0)*sqrt(u)/sqrt( (u-1.0)*(u+1.0)*(u+2.0)*(2.0*u-1.0)*(2.0*u+3.0) );
	t1 = a1*j3*sqrt(15.0)*E*E/(4.0*(M/MeV)*(M/MeV));
	t2 = rho_uv/(u+1.0)*(5.0*E*E)/(4.0*(M/MeV)*(M/MeV))*(sqrt(3)*c1*g + sqrt(2)*c1*j2);
	t3 = -1.0*sigma_uv/(u+1.0)*(5.0*E*E)/(2.0*(M/MeV)*(M/MeV))*c1*j3;
	
	/*
	f1 = -1.0*Params->delta_uv*(3.0*Params->u*Params->u + 3.0*Params->u-1.0)*sqrt(Params->u)/sqrt( (Params->u-1.0)*(Params->u+1.0)*(Params->u+2.0)*(2.0*Params->u-1.0)*(2.0*Params->u+3.0) );
	t1 = a1*j3*sqrt(15.0)*E*E/(4.0*(M/MeV)*(M/MeV));
	t2 = Params->rho_uv/(Params->u+1.0)*(5.0*E*E)/(4.0*(M/MeV)*(M/MeV))*(sqrt(3)*c1*g + sqrt(2)*c1*j2);
	t3 = -1.0*Params->sigma_uv/(Params->u+1.0)*(5.0*E*E)/(2.0*(M/MeV)*(M/MeV))*c1*j3;
	*/
	
	double F3 = f1*t1 + t2 + t3;
	if(verbose) { cout << "F3 = " << F3 << endl; }
	return F3;
}

double Holstein52Generator::FF_dF1_Euvs(double E)  
// CC to the F_i(E, J, J', 0) functions of (B9).  E in units of MeV.  
// I've used 'a1' in place of 'a' (same thing for c) for these expressions --
// so the coulomb corrections don't get recoil-order corrections.
{
//	return 0.0;
	double prefactor = +(8.0*alpha*Z_daughter)/(3.0*M_PI);
	double a_term = a1*a1 * (4.0*E*(X_coulomb+Y_coulomb) + E0*X_coulomb + (m_e*m_e/E)*(X_coulomb+2.0*Y_coulomb) );
	double c_term = c1*c1 * (E*(16.0/3.0*X_coulomb + 4.0*Y_coulomb) - 1.0/3.0*E0*X_coulomb + (m_e*m_e/E)*(X_coulomb+2.0*Y_coulomb) );
	
	double dF = prefactor*(a_term + c_term);
	
//	return sign*(-8.)*alpha*Z_d/3./M_PI* 
//		( sq(a1_tot)*(4.*E*(X+Y) + E0*X + me2*(X+2.*Y)/E ) +
//		  sq(c1_tot)*(4.*E*(4.*X/3.+Y) - E0*X/3. + me2*(X+2.*Y)/E )  );
	//
	return dF;
}
double Holstein52Generator::FF_dF4_Euvs(double E)
{
//	return 0.0;
	double delta_uv   = Params->get_delta_uv();
	double gamma_uv   = Params->get_gamma_uv();
	double u          = 3.0/2.0;
	
	double prefactor = +(8.0*alpha*Z_daughter)/(3.0*M_PI) * E*(5.0*X_coulomb+4.0*Y_coulomb);
	double ac_term = delta_uv * sqrt( u/(u+1.0) ) * 2.0*a1*c1;
	double c_term  = (+1.0)*gamma_uv/(u+1.0)*c1*c1;
	
	double dF = prefactor*(ac_term + c_term);
	
//	double tmp1 =0;
//	double tmp2 =0;
//	double tmp3 =0;
//	tmp1 = sign*(-8.)*alpha*Z_d*E*(5.*X+4.*Y)/3./M_PI;
//	tmp2 = sqrt(J/(J+1.))*2.*a1_tot*c1_tot;
//	tmp3 = sign*gamma_uv()*sq(c1_tot)/(J+1.);
//	
//	return tmp1*(tmp2-tmp3);
	
	return dF;
}
double Holstein52Generator::FF_dF7_Euvs(double E)
{
//	return 0.0;
	double delta_uv   = Params->get_delta_uv();
	double gamma_uv   = Params->get_gamma_uv();
	double u          = 3.0/2.0;
	
	double prefactor = +(8.0*alpha*Z_daughter)/(3.0*M_PI) * (E0-E)*X_coulomb;
	double ac_term = delta_uv * sqrt( u/(u+1.0) ) * 2.0*a1*c1;
	double c_term  = (-1.0)*gamma_uv/(u+1.0)*c1*c1;
	
	double dF = prefactor*(ac_term + c_term);
	
//	double tmp1 =0;
//	double tmp2 =0;
//	double tmp3 =0;
//	tmp1 = sign*(-8.)*alpha*Z_d/3./M_PI*(E0-E)*X;
//	tmp2 = sqrt(J/(J+1.))*2.*a1_tot*c1_tot;
//	tmp3 = sign*gamma_uv()*sq(c1_tot)/(J+1.);
//	
//	return tmp1*(tmp2+tmp3);

	return dF;
}

double Holstein52Generator::dF0_coulomb(double E) // CC to the F_i(E) functions of (B10).  E in units of MeV.
{
	return FF_dF1_Euvs(E);
}
double Holstein52Generator::dF1_coulomb(double E)
{
	return ( FF_dF4_Euvs(E) + 1.0/3.0 * FF_dF7_Euvs(E) );
}


G4double Holstein52Generator::pbeta(G4double E)  // *only* in get_probability so far...
{
	G4double pbetac = sqrt( pow(E/MeV,2) - pow((m_e/MeV),2) );  // MeV for everything...
	return pbetac;  // comes out in G4units of MeV too.
}
G4double Holstein52Generator::get_Ebeta(G4double pbeta)  // not used???
{
	G4double E = sqrt( pow(pbeta/MeV,2) + pow((m_e/MeV), 2) );
	return E;
}
double Holstein52Generator::get_v_from_p(G4double pbeta) // [v] = unitless m/s, pbeta in G4units of MeV.
{
	double beta = (pbeta/MeV)* sqrt( (1.0)/((m_e/MeV)*(m_e/MeV) + (pbeta/MeV)*(pbeta/MeV)) );
	double v = beta*speed_of_light;  // meters/second 
	return v;
}
G4double Holstein52Generator::get_p_from_v(double vbeta)  // ...
{
	double gamma = 1.0 / sqrt(1.0 - vbeta*vbeta/(speed_of_light*speed_of_light)); // unitless....
	G4double pbeta = gamma*(m_e/MeV)*vbeta/speed_of_light * MeV;
	return pbeta;
}

double Holstein52Generator::FermiFunction( double E_kin_MeV )
{
//	bool verbose = false;
	FF_val = 1.0;
	if(!do_coulomb) { return FF_val; }
	
	// otherwise...
	FF_val = the_FF->getVFF(E_kin_MeV);
	
//	if(verbose)
//	{
//		G4cout << "T_beta = " << E_kin_MeV << " MeV;\t F(Z,E) = " << FF_val << G4endl;
//	}
	
	return FF_val;
}

double Holstein52Generator::get_probability(G4double E, double costheta) // want to use costheta_lab
{
	int verbose = 0;
	
	double Gv2       = 1.0; // kludge
	double cos2theta = 1.0; // kludge
	double prefactor = 2.0*Gv2*cos2theta/( pow(2.0*pi, 4) );
	
//	double the_scaling = prefactor*FermiFunction( Params->get_Z_daughter(), (E/MeV) ) * pow((E0/MeV)-(E/MeV), 2)*(E/MeV)*pbeta(E);
	// Call FF using KE rather than Etot.
//	double the_scaling = prefactor*FermiFunction( (E-m_e)/MeV ) * pow((E0/MeV)-(E/MeV), 2) * (E/MeV) * pbeta(E);
	double the_scaling = prefactor * pow((E0/MeV)-(E/MeV), 2) * (E/MeV) * pbeta(E);
	//Ebeta_kin_MeV = (Ebeta - m_e)/MeV;
	
	
	if(verbose>1)
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
	
	double t0 = (F_0(E/MeV)*FermiFunction( (E-m_e)/MeV ) + dF0_coulomb(E/MeV) );
	double t1 = (F_1(E/MeV)*FermiFunction( (E-m_e)/MeV ) + dF1_coulomb(E/MeV) ) * ( Lambda1 * costheta*(pbeta(E)/MeV)/(E/MeV) );
	double t2 = (F_2(E/MeV)*FermiFunction( (E-m_e)/MeV ) )                      * ( Lambda2 * (pbeta(E)/MeV)*(pbeta(E)/MeV) / ((E/MeV)*(E/MeV)) * (costheta*costheta - 1.0/3.0) );
	double t3 = (F_3(E/MeV)*FermiFunction( (E-m_e)/MeV ) )                      * ( Lambda3 * (pow(costheta*(pbeta(E)/MeV)/(E/MeV), 3) - 3.0/5.0*costheta*pow((pbeta(E)/MeV)/(E/MeV), 3)) );
	
	double g_V = Params -> get_g_Vector();
	double g_A = Params -> get_g_Axial();
	double g_S = Params -> get_g_Scalar();
	double g_T = Params -> get_g_Tensor();
	double b_term = (m_e/E) * FermiFunction( (E-m_e)/MeV ) * (-2.0*jtw_gamma) * ( M_F*M_F*g_V*g_S + M_GT*M_GT*g_A*g_T  );
	
//	double kludge_b_term = jtw_gamma*b_kludge*( (m_e/MeV)/(E/MeV) );
//	double the_prob      = the_scaling*(t0+t1+t2+t3);
	holstein_Abeta       = F_1(E/MeV)/t0;  // Do I still believe this, what with the coulomb corrections?
	holstein_probability = the_scaling*(t0+t1+t2+t3+b_term);
	
	if(verbose>0)
	{
	//	holstein_Abeta       = F_1(E/MeV)/t0;  // Do I still believe this, what with the coulomb corrections?  Also, does this var. do anything??
		
		cout << "Holstein52Generator::get_probability:  the_prob = " << holstein_probability << endl;
	//	cout << "\tthe_scaling = " << the_scaling << endl;
	//	cout << "\tholstein_Abeta = " << holstein_Abeta << endl;
	//	cout << "\tt0 = " << t0 << endl;
	//	cout << "\tt1 = " << t1 << endl;
	//	cout << "\tt2 = " << t2 << endl;
	//	cout << "\tt3 = " << t3 << endl;
		cout << "T_beta = " << (E-m_e)/MeV << " MeV;\t F(Z,E) = " << FF_val << G4endl;
		cout << "F_0 = " << F_0(E/MeV) << ";\tdF0 = " << dF0_coulomb(E/MeV) << endl;
		cout << "F_1 = " << F_1(E/MeV) << ";\tdF1 = " << dF1_coulomb(E/MeV) << endl;
	}
	return holstein_probability;
}

double Holstein52Generator::get_jtw_probability(G4double E, double costheta)
{
	double Gv2       = 1.0; // kludge
	double cos2theta = 1.0; // kludge
	double prefactor = 2.0*Gv2*cos2theta/( pow(2.0*pi, 4) );
//	double the_scaling = prefactor*FermiFunction(Params->get_Z_daughter(), (E/MeV))*pow((E0/MeV)-(E/MeV), 2)*(E/MeV)*pbeta(E);
	double the_scaling = prefactor * FermiFunction( (E-m_e)/MeV ) * pow((E0/MeV)-(E/MeV), 2)*(E/MeV)*pbeta(E);
	
	double fake_F0 = (jtw_a1*jtw_a1 + jtw_c1*jtw_c1);
	double fake_F1 = jtw_c1*jtw_c1*(1.0/(3.0/2.0 + 1.0)) + 2.0*jtw_a1*jtw_c1*sqrt((3.0/2.0)/(3.0/2.0 + 1.0));
	
//	double kludge_b_term = jtw_gamma*b_kludge*( (m_e/MeV)/(E/MeV) );  // have I definitely got the units right here??  ..if b_kludge==0, this term is zero too.
	
//	double the_term = fake_F0 + fake_F1 * Lambda1 * costheta*(pbeta(E)/MeV)/(E/MeV) + kludge_b_term;
	double the_term = fake_F0 + fake_F1 * Lambda1 * costheta*(pbeta(E)/MeV)/(E/MeV) + 0.0;
	double the_prob = the_scaling*the_term;
	
	jtw_probability = the_prob;
	return jtw_probability;
}

bool Holstein52Generator::shoot_decayevent( G4double the_monoenergy ) // uses param. runfast to decide whether to check *both* acceptances if one has already failed.
{
	bool verbose        = false;
//	jtw_acceptance      = false;
	holstein_acceptance = false;
//	det_acceptance      = false;
	
	this -> randomize_nuclear(false);
//	print_vars();  // get values of E0 and m_e that we're using.
	// Instead of this->randomize_atomic(...), just do the things from the function.
	// Later, for efficiency, I should avoid doing this every single event, but for now it's fine.
	Mz  = the_atomic_setup->GetPops()->get_Mz();
	Mz2 = the_atomic_setup->GetPops()->get_Mz2();
	Mz3 = the_atomic_setup->GetPops()->get_Mz3();
	initialize_lambdafuncs();  // this depends on the sublevel populations, which should be set up above.
	//
	did_the_printing=0;
	while( !holstein_acceptance )
	{
		this -> randomize_direction(the_monoenergy); // pick Ebeta, costheta (using the cone, if we're doing that.)  If the monoenergy is the default (-10), we don't use it.
		this -> check_holstein_acceptance();
		
		if(did_the_printing==0)
		{
		//	cout << "{" << costheta_lab << ", " << phi << "}, " << endl;
		//	cout << "{" << costheta_lab << ", " << test_prob << "}, " << endl;  // this one is good...
		//	cout << "wtf?" << endl;
		//	did_the_printing=1;
		}
		did_the_printing++;
	}
	
	if(verbose)
	{
	//	cout << "pdf_acceptance  = " << int(pdf_acceptance) << endl;
	//	cout << "holstein_acceptance = " << int(holstein_acceptance) << endl;
	//	cout << "jtw_acceptance = " << int(jtw_acceptance) << endl;
	//	cout << "det_acceptance  = " << int(det_acceptance) << endl;
	//	cout << "Event accepted? = " << int(det_acceptance && (holstein_acceptance || jtw_acceptance)) << endl;
	}
	return holstein_acceptance;
}

void Holstein52Generator::randomize_direction( G4double the_monoenergy ) // saves initial_momentum and initial_velocity in class.  Also Ebeta.  Checks the cone, and is optimized for that.
{
	bool verbose=false;
	
	// pick E
	if(the_monoenergy == -10)
	{
	//	cout << "m_e/MeV = " << m_e/MeV << endl;
	//	cout << "E0/MeV = " << E0/MeV << endl;
		Ebeta = G4RandFlat::shoot(m_e/MeV, E0/MeV)*MeV;  // m_e and E0 have to be in the same units.  also, in units of MeV.
	}
	else
	{
	//	Ebeta = the_monoenergy;
		Ebeta_kin_MeV = the_monoenergy/MeV;
		Ebeta = Ebeta_kin_MeV*MeV + m_e;
	
	//	MeV*Ebeta_kin_MeV + m_e= Ebeta;
	//	Ebeta = 
	}
	
	// pick costheta
	costheta_lab = G4RandFlat::shoot(-1.0, 1.0);
	while( use_cone && (abs(costheta_lab) < cone_costheta) ) // the max. is on theta.  cos(theta_max) is the min.  this is a really slow way to do this.
	{
	//	cout << "shooting somemore." << endl;
		costheta_lab = G4RandFlat::shoot(-1.0, 1.0);
	}
//	cout << "Picked cos(theta) = " << costheta_lab << endl;

	// pick phi
//	double phi = G4RandFlat::shoot(0.0, 2.0*pi);
	phi = G4RandFlat::shoot(0.0, 2.0*pi);  // I've only made phi get saved because I want to print it out from another function for debugging.  it's not really needed.
	
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
		cout << "(Ebeta/MeV)                                               = " << (Ebeta/MeV) << endl;
		cout << "pbeta(Ebeta)                                              = " << pbeta(Ebeta) << endl;
		cout << "get_v_from_p( pbeta(Ebeta) )                              = " << get_v_from_p( pbeta(Ebeta) ) << endl;
		cout << "    initial_velocity.mag()                                = " << initial_velocity.mag() << endl;
		cout << "    vx/c=" << initial_velocity.x()/speed_of_light << ";\tvy/c=" << initial_velocity.y()/speed_of_light << ";\tvz/c=" << initial_velocity.z()/speed_of_light << endl;
		cout << "    sqrt(vx^2 + vy^2 + vz^2)                              = ";
		cout << sqrt( initial_velocity.x()*initial_velocity.x() 
					+ initial_velocity.y()*initial_velocity.y() 
					+ initial_velocity.z()*initial_velocity.z() ) << endl;

		cout << "get_p_from_v( get_v_from_p( pbeta(Ebeta) ) )/MeV          = " << get_p_from_v( get_v_from_p( pbeta(Ebeta)) )/MeV << endl;
		cout << "    initial_momentum.mag()/MeV                            = " << initial_momentum.mag()/MeV << endl;
		cout << "    px/MeV=" << initial_momentum.x()/MeV << ";\tpy/MeV=" << initial_momentum.y()/MeV << ";\tpz/MeV=" << initial_momentum.z()/MeV  << endl;
		cout << "    sqrt(px^2 + py^2 + pz^2)                              = ";
		cout << sqrt( initial_momentum.x()/MeV*initial_momentum.x()/MeV 
					+ initial_momentum.y()/MeV*initial_momentum.y()/MeV 
					+ initial_momentum.z()/MeV*initial_momentum.z()/MeV ) << endl;
		cout << "get_Ebeta( get_p_from_v( get_v_from_p( pbeta(Ebeta) ) ) )/MeV = " << get_Ebeta( get_p_from_v( get_v_from_p( pbeta(Ebeta)) ) )/MeV << endl;
	}
	
//	if(did_the_printing>=3)
//	{
//	//	cout << "{" << costheta_lab << ", ";// << phi << "}, " << endl;
//	//	cout << costheta_lab << "}, " << endl;
//		cout << "{" << costheta_lab << ", " << phi << "}, " << endl;
//	//	cout << "{" << costheta_lab << ", " << Ebeta/MeV << "}, " << endl;
//	}

//	cout << "{" << costheta_lab << ", " << phi << "}, " << endl;
//	cout << "{" << costheta_lab << ", " << Ebeta/MeV << "}, " << endl;
	
	return;
}

bool Holstein52Generator::check_holstein_acceptance()  // uses initial_momentum
{
	bool verbose = false;
	
//	pdf_acceptance      = false;
//	jtw_acceptance      = false;
	holstein_acceptance = false;
	
	if(initial_momentum.mag() == 0)
	{
		cout << "No, this probably hasn't been initialized.  Rejected!"  << endl;
		assert(0);
	}
//	double prob_max = 0.2;  
//	double test_prob = G4RandFlat::shoot(0.0, prob_max);
	test_prob = G4RandFlat::shoot(0.0, prob_max);
	
	holstein_probability = get_probability(Ebeta, costheta_lab ); // this is redundant, but whatevs.
	
//	if(verbose)
//	{
//		cout << "E = " << Ebeta << ";  costheta=" << costheta_lab << ";  holstein_prob = " << holstein_probability << ";  test_prob = " << test_prob << endl;
//	//	cout << "holstein_prob = " << holstein_probability << /*";  jtw_prob = " << jtw_probability << */";  max = " << prob_max << ";  test_prob = " << test_prob << endl;
//	}
	//
	if(holstein_probability > prob_max)
	{
		cout << "Bad!  We've broken our acceptance/rejection method!" << endl;
		cout << "You should increase prob_max within the code.  " << endl;
		cout << "prob_max = " << prob_max << endl;
		assert(0); // hard kill.
		return true;
	}
	
	// Check Holstein Acceptance:
	if( test_prob <= holstein_probability)
	{
		holstein_acceptance = true;
	//	if(verbose) { cout << "\tEvent accepted by Holstein PDF!" << endl; }
	}
	else
	{
		holstein_acceptance = false;
	//	if(verbose) { cout << "\tEvent rejected by Holstein PDF." << endl; }
	}
	
	return holstein_acceptance;
}

bool Holstein52Generator::check_detector_acceptance(G4ThreeVector initial_thing)  // probably only used in TheHolstein.
{
	initial_position = initial_thing;
	
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

void Holstein52Generator::print_results()
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
