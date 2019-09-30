#include <string>  
#include <sstream> 


using std::string;
using std::vector;


struct sublevel_populations  // Spin 3/2 only!!
{
public:
	sublevel_populations(int the_sigma);  // completely polarized, so the_sigma=+/-1.
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
