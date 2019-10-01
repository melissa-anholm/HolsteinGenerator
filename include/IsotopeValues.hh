#include <string> 
#include <iostream>  // cout, endl
#include <iomanip>   // setw

using std::string;


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
		std::cout << std::setw(30) << name;
		std::cout << std::setw(13) << value << std::setw(12) << uncertainty;
		std::cout << "    " << message << std::endl;
	};
	
private:
	double value;
	double uncertainty;
	string message;
	string name;  // added by MJA.  It probably isn't necessary.  the info is elsewhere.
};
