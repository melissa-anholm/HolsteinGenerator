// Authors: Spencer Behling, Benjamin Fenker, and Melissa Anholm  - 2013

#ifndef SplitString_h
#define SplitString_h 1

#include <string>
#include <sstream>
#include <vector>

using std::cout;
using std::endl;

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
#endif
