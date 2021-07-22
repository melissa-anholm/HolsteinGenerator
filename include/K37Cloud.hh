// Authors: Spencer Behling, Benjamin Fenker, Melissa Anholm - 2012

#ifndef K37Cloud_h
#define K37Cloud_h 1

#include "G4ParticleDefinition.hh" // calls:  #include <CLHEP/Units/SystemOfUnits.h>
#include "G4SystemOfUnits.hh" // calls:  #include <CLHEP/Units/SystemOfUnits.h> and then lets us use unit names for free.
#include "G4ThreeVector.hh"
#include "globals.hh"


class K37Cloud 
{
public:
	K37Cloud();
	~K37Cloud();

public:
	void PrintCloud();
	
	// For the whole cloud:
	G4ThreeVector GetInitialCloudPosition() { return initial_position;           }
	G4ThreeVector GetFinalCloudPosition()   { return final_position;             }
	
	G4ThreeVector GetInitialCloudSize()     { return initial_cloud_size;         }
	G4ThreeVector GetFinalCloudSize()       { return final_cloud_size;           }
	
	G4double GetFreeExpansionTime()         { return expansion_before_polarized; }
	G4double GetOP_CycleTime()              { return cycleTime;                  }

//	G4ThreeVector GetCloudCenter()      { return initial_position;           }
//	G4ThreeVector GetTemperature()      { return temperature;                }
//	G4ThreeVector GetSailVelocity()     { return sail_velocity;              }
//	G4ThreeVector GetVelocitySigma()    { return velocity_sigma;             }  // must be already set up with temperature?
//	void SetCloudCenter(G4ThreeVector center);
	
	void SetInitialCloudPosition(G4ThreeVector center);
	void SetFinalCloudPosition(G4ThreeVector center);
	
	void SetInitialCloudSize(G4ThreeVector size);
	void SetFinalCloudSize(G4ThreeVector size);
	
	void SetInitialCloudSize(G4double size) { SetInitialCloudSize(G4ThreeVector(size, size, size)); };
	void SetFinalCloudSize(G4double size)   { SetFinalCloudSize(G4ThreeVector(size, size, size)); };
	
	void SetFreeExpansionTime(G4double time);
	void SetOP_CycleTime(G4double time);
	
	
	void SetDecayPosition(G4ThreeVector center)
	{
		SetInitialCloudPosition(center);
		SetFinalCloudPosition(center);
		SetInitialCloudSize(0.0);
		SetFinalCloudSize(0.0);
	}
	
	void SetGenH1(bool doit) 
	{
		gen_h1 = doit;
		if(gen_h1) 
		{ 
			gen_h2=false;
			gen_h3=false;
			gen_h4=false;
			gen_rmcp=false;
			gen_ring=false;
			//
			SetInitialCloudPosition( G4ThreeVector( (45.0+1.0)*mm, 97.0*mm, (45.0+1.0)*mm ) );
			SetFinalCloudPosition(   G4ThreeVector( (45.0+1.0)*mm, 97.0*mm, (45.0+1.0)*mm ) );
			SetInitialCloudSize(     G4ThreeVector(1.0*mm, 0.0, 1.0*mm) );
			SetFinalCloudSize(       G4ThreeVector(1.0*mm, 0.0, 1.0*mm) );
			
			double max_opposite     = 47.0;
			double max_hypotenuse   = sqrt( 47.0*47.0 + (97.0)*(97.0) );
			//
			double min_opposite     = 45.0;
			double min_hypotenuse   = sqrt( 45.0*45.0 + (97.0)*(97.0) );
			
			max_theta = asin(max_opposite / max_hypotenuse);
			min_theta = asin(min_opposite / min_hypotenuse);
			
			if(gen_hoopedge)
			{
				SetInitialCloudPosition( G4ThreeVector( 45.0*mm, 97.5*mm, 45.0*mm ) );
				SetFinalCloudPosition(   G4ThreeVector( 45.0*mm, 97.5*mm, 45.0*mm ) );
				SetInitialCloudSize(     G4ThreeVector(0, 0.5*mm, 0) );
				SetFinalCloudSize(       G4ThreeVector(0, 0.5*mm, 0) );
				
				double min_opposite     = 45.0;
				double min_hypotenuse   = sqrt( 45.0*45.0 + (97.0)*(97.0) );
				min_theta = asin(min_opposite / min_hypotenuse);
				max_theta = min_theta;
			}
			adjacent = 0;  // shouldn't get used.
		}
	}
	void SetGenH2(bool doit) 
	{
		gen_h2 = doit;
		if(gen_h2) 
		{ 
			gen_h1=false;
			gen_h3=false;
			gen_h4=false;
			gen_rmcp=false;
			gen_ring=false;
			//
			SetInitialCloudPosition( G4ThreeVector( (40.0+1.5)*mm, 75.0*mm, (40.0+1.5)*mm ) );
			SetFinalCloudPosition(   G4ThreeVector( (40.0+1.5)*mm, 75.0*mm, (40.0+1.5)*mm ) );
			SetInitialCloudSize(     G4ThreeVector(1.5*mm, 0.0, 1.5*mm) );
			SetFinalCloudSize(       G4ThreeVector(1.5*mm, 0.0, 1.5*mm) );
			
			double max_opposite     = 43.0;
			double max_hypotenuse   = sqrt( 43.0*43.0 + (75.0)*(75.0) );
			//
			double min_opposite     = 40.0;
			double min_hypotenuse   = sqrt( 40.0*40.0 + (75.0)*(75.0) );
			
			max_theta = asin(max_opposite / max_hypotenuse);
			min_theta = asin(min_opposite / min_hypotenuse);
			
			if(gen_hoopedge)
			{
				SetInitialCloudPosition( G4ThreeVector( 40.0*mm, 75.5*mm, 40.0*mm ) );
				SetFinalCloudPosition(   G4ThreeVector( 40.0*mm, 75.5*mm, 40.0*mm ) );
				SetInitialCloudSize(     G4ThreeVector(0, 0.5*mm, 0) );
				SetFinalCloudSize(       G4ThreeVector(0, 0.5*mm, 0) );
				
				double min_opposite     = 40.0;
				double min_hypotenuse   = sqrt( 40.0*40.0 + (75.0)*(75.0) );
				min_theta = asin(min_opposite / min_hypotenuse);
				max_theta = min_theta;
			}
			adjacent = 0;  // shouldn't get used.
		}
	}
	void SetGenH3(bool doit) 
	{
		gen_h3 = doit;
		if(gen_h3) 
		{ 
			gen_h1=false;
			gen_h2=false;
			gen_h4=false;
			gen_rmcp=false;
			gen_ring=false;
			
			SetInitialCloudPosition( G4ThreeVector(0., 57.0*mm, (37.9+0.75)*mm) );
			SetFinalCloudPosition(   G4ThreeVector(0., 57.0*mm, (37.9+0.75)*mm) );
			SetInitialCloudSize(     G4ThreeVector(45.0*mm, 0.0, 0.75*mm) );
			SetFinalCloudSize(       G4ThreeVector(45.0*mm, 0.0, 0.75*mm) );
			
			double max_opposite     = 45.0;
			double max_hypotenuse   = sqrt( 57.0*57.0 + (37.9+0.75)*(37.9+0.75) + 45.0*45.0 );
			//
			adjacent                = sqrt( 57.0*57.0 + (37.9+0.75)*(37.9+0.75) );
			max_theta = asin(max_opposite / max_hypotenuse);
			min_theta = 0;  // shouldn't get used.
			//
			if(gen_hoopedge)
			{
				SetInitialCloudPosition( G4ThreeVector(0., (57.0+0.5)*mm, (37.9)*mm) );
				SetFinalCloudPosition(   G4ThreeVector(0., (57.0+0.5)*mm, (37.9)*mm) );
				SetInitialCloudSize(     G4ThreeVector(45.0*mm, 0.5*mm, 0.0*mm) );
				SetFinalCloudSize(       G4ThreeVector(45.0*mm, 0.5*mm, 0.0*mm) );
				
				double max_opposite     = 45.0;
				double max_hypotenuse   = sqrt( 57.5*57.5 + (37.9)*(37.9) + 45.0*45.0 );
				adjacent                = sqrt( 57.5*57.5 + (37.9)*(37.9) );
				max_theta = asin(max_opposite / max_hypotenuse);
				min_theta = 0;  // shouldn't get used.
			}
		}
	};
	void SetGenH4(bool doit) 
	{
		gen_h4 = doit;
		if(gen_h4) 
		{ 
			gen_h1=false;
			gen_h2=false;
			gen_h3=false;
			gen_rmcp=false;
			gen_ring=false;
			
			SetInitialCloudPosition( G4ThreeVector(0., 29.0*mm, (37.9+0.75)*mm) );
			SetFinalCloudPosition(   G4ThreeVector(0., 29.0*mm, (37.9+0.75)*mm) );
			SetInitialCloudSize(     G4ThreeVector(45.0*mm, 0.0, 0.75*mm) );
			SetFinalCloudSize(       G4ThreeVector(45.0*mm, 0.0, 0.75*mm) );
			
			double max_opposite     = 45.0;
			double max_hypotenuse   = sqrt( 29.0*29.0 + (37.9+0.75)*(37.9+0.75) + 45.0*45.0 );
			adjacent                = sqrt( 29.0*29.0 + (37.9+0.75)*(37.9+0.75) );
			max_theta = asin(max_opposite / max_hypotenuse);  // this will always be positive.
			min_theta = 0;  // shouldn't get used.
			
		//	// max_cos_theta = 1.0;
		//	// pick a cos(theta), then hyp = adjacent/cos(theta) -->
		//	// ( adjacent/cos(theta) )^2 = adj^2 + x^2
		//	// x^2 = ( adjacent / cos(theta) )^2 - adjacent^2;
			
			//
			if(gen_hoopedge)
			{
				SetInitialCloudPosition( G4ThreeVector(0., (29.0+0.5)*mm, (37.9)*mm) );
				SetFinalCloudPosition(   G4ThreeVector(0., (29.0+0.5)*mm, (37.9)*mm) );
				SetInitialCloudSize(     G4ThreeVector(45.0*mm, 0.5*mm, 0.0*mm) );
				SetFinalCloudSize(       G4ThreeVector(45.0*mm, 0.5*mm, 0.0*mm) );
				
				double max_opposite     = 45.0;
				double max_hypotenuse   = sqrt( 29.5*29.5 + (37.9)*(37.9) + 45.0*45.0 );
				adjacent                = sqrt( 29.5*29.5 + (37.9)*(37.9) );
				max_theta = asin(max_opposite / max_hypotenuse);  // this will always be positive.
				min_theta = 0;  // shouldn't get used.
			}
		}
	};
	
	void SetGenRMCP(bool doit)
	{
		gen_rmcp = doit;
		if(gen_rmcp) 
		{ 
			gen_h1=false;
			gen_h2=false;
			gen_h3=false;
			gen_h4=false;
			gen_ring=false;
			
			// should be 101.48 instead of 104.4
			
			SetInitialCloudPosition( G4ThreeVector(0., 104.4*mm, 0.0) );
			SetFinalCloudPosition(   G4ThreeVector(0., 104.4*mm, 0.0) );
			SetInitialCloudSize(     G4ThreeVector(41.5*mm, 0.0, 41.5*mm) );
			SetFinalCloudSize(       G4ThreeVector(41.5*mm, 0.0, 41.5*mm) );
			
			double max_opposite     = 41.5;
			double adjacent         = 104.4;
		//	double max_hypotenuse   = sqrt( max_opposite*max_opposite + (104.4)*(104.4) );
			
		//	max_theta = asin(max_opposite / max_hypotenuse);  // this will always be positive.
			max_theta = atan(max_opposite / adjacent);
			min_theta = 0;  // shouldn't get used.
		}
	}
	void SetGenRing(bool doit)
	{
		gen_ring = doit;
		if(gen_ring) 
		{ 
			gen_h1=false;
			gen_h2=false;
			gen_h3=false;
			gen_h4=false;
			gen_rmcp=false;
			
			// should be 99.48 instead of 102.4
			
			SetInitialCloudPosition( G4ThreeVector( (41.5+0.75)*mm, 102.4*mm, (41.5+0.75)*mm ) );
			SetFinalCloudPosition(   G4ThreeVector( (41.5+0.75)*mm, 102.4*mm, (41.5+0.75)*mm ) );
			SetInitialCloudSize(     G4ThreeVector(0.75*mm, 0.0, 0.75*mm) );
			SetFinalCloudSize(       G4ThreeVector(0.75*mm, 0.0, 0.75*mm) );
			
			//
			double max_opposite     = 43.0;
			double max_hypotenuse   = sqrt( 43.0*43.0 + (102.4)*(102.4) );
			//
			double min_opposite     = 41.5;
			double min_hypotenuse   = sqrt( 41.5*41.5 + (102.4)*(102.4) );
			
			max_theta = asin(max_opposite / max_hypotenuse);
			min_theta = asin(min_opposite / min_hypotenuse);
			
			if(gen_hoopedge)
			{
				// should be 100.48 instead of 103.4
				
				SetInitialCloudPosition( G4ThreeVector( 41.5*mm, 103.4*mm, 41.5*mm ) );
				SetFinalCloudPosition(   G4ThreeVector( 41.5*mm, 103.4*mm, 41.5*mm ) );
				SetInitialCloudSize(     G4ThreeVector(0, 1.0*mm, 0) );
				SetFinalCloudSize(       G4ThreeVector(0, 1.0*mm, 0) );
				
				double min_opposite     = 41.5;
				double min_hypotenuse   = sqrt( 41.5*41.5 + (103.4)*(103.4) );
				min_theta = asin(min_opposite / min_hypotenuse);
				max_theta = min_theta;
			}
			adjacent = 0;  // shouldn't get used.
		}
	}
	
	void SetUseEdge(bool doit)
	{
		gen_hoopedge=doit;
		if(gen_h1)
		{
			SetGenH1(true);
		}
		if(gen_h2)
		{
			SetGenH2(true);
		}
		if(gen_h3)
		{
			SetGenH3(true);
		}
		if(gen_h4)
		{
			SetGenH4(true);
		}
		if(gen_rmcp)
		{
			// doesn't matter.  no edge.
		}
		if(gen_ring)
		{
			SetGenRing(true);
		}
	}
	
	bool is_h1()       { return gen_h1; };
	bool is_h2()       { return gen_h2; };
	bool is_h3()       { return gen_h3; };
	bool is_h4()       { return gen_h4; };
	bool is_rmcpgen()         { return gen_rmcp; };
	bool is_ceramicringgen()  { return gen_ring; }
	bool is_edgegen()         { return gen_hoopedge; };
	
	double get_adjacent()     { return adjacent; };
	double get_maxtheta()     { return max_theta; };
	double get_mintheta()     { return min_theta; };
	
private:
	bool gen_h1;
	bool gen_h2;
	bool gen_h3;
	bool gen_h4;
	bool gen_rmcp;
	bool gen_ring;
	bool gen_hoopedge;
	
	double max_theta;
	double min_theta;
	double adjacent; // NO UNITS!
//	double r_from_trap;
	
	// Given that I have an initial cloud size and initial cloud position, 
	// I might set up the cloud's evolution by Temperature and Sail Velocity, 
	// or equivalently I might use Final Size and Final Position.  
	// It would be a mistake to try to code both in.
	// Ben used temperature and sail velocity, but I don't really like it.
	// I might change it later.
	/*
	void SetTemperature(G4ThreeVector temp);
	void SetTemperature(G4double temp);
	void SetSailVelocity(G4ThreeVector vel);
	*/


	
private:
//	void set_up_sail_velocity();  // initial_position, final_position, and cycleTime must already be set up.
//	void set_up_temperature();    // must have initial_size, final_size, and cycleTime already set up.
	
	

private:
// For the whole cloud:
	G4ThreeVector initial_position;
	G4ThreeVector final_position; //
	
	G4ThreeVector initial_cloud_size;
	G4ThreeVector final_cloud_size;
	
//	G4ThreeVector temperature;
//	G4ThreeVector velocity_sigma;
	G4double cycleTime;
	G4double expansion_before_polarized;
	
//	G4ThreeVector sail_velocity;  // legacy.
	
public:
//	void SetupVelocitySigma(G4ThreeVector temperature);  // We might like this to be private, but our instantiation of this class is itself private, so wev.
//	G4double CalcSigma(G4double temperature);            // We might like this to be private, but our instantiation of this class is itself private, so wev.
//	void Initialize();                                   // We might like this to be private, but our instantiation of this class is itself private, so wev.
//private:	
//	G4bool initialize_complete_;
	
};

#endif
