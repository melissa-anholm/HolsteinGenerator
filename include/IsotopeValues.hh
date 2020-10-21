// Authors: Spencer Behling, Benjamin Fenker, and Melissa Anholm  - 2013

#ifndef IsotopeValues_h
#define IsotopeValues_h 1

#include <string> 
#include <iostream>  // cout, endl
#include <iomanip>   // setw
#include <sstream>
#include <vector>
#include <map>
#include "fstream"   // 

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;


class isotope_values
{
public:
	isotope_values(const double &value_, const double &uncertainty_, const string &message_, const string &name_)
	:value(value_), uncertainty(uncertainty_), message(message_), name(name_)
	{};
	
	double GetValue()const       { return value;       };
	double GetUncertainty()const { return uncertainty; };
	string GetMessage()const     { return message;     };     // extra added by MJA.
	string GetName()const        { return name;        };     // extra added by MJA.
	
	void SetUncertainty(const double &uncertainty_)
		{ uncertainty = uncertainty_; };
	void SetValue(const double &value_)
		{ value = value_; };
	void Print() 
	{
		cout << std::setw(30) << name;
		cout << std::setw(13) << value << std::setw(12) << uncertainty;
		cout << "    " << message << endl;
	};
	
private:
	double value;
	double uncertainty;
	string message;
	string name;  // added by MJA.  It probably isn't necessary.  the info is elsewhere.
};


namespace SS
{
	// below:  does this one *ever* get used?! ...yes, yes it does.
	std::vector<std::string> & split(const std::string &s, const char &delim, std::vector<std::string> &elems);
	
	// below:  this one is used.
	std::vector<std::string> split(const std::string &s, const char &delim);
	
	map<string, isotope_values * > loadup_textfile(string paramfilename);
	
}


#endif
