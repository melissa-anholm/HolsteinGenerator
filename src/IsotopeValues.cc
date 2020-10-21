// Authors: Spencer Behling, Benjamin Fenker, and Melissa Anholm  - 2013

#undef NDEBUG
#include<assert.h>


#include "IsotopeValues.hh"

/*
#ifndef SplitString_h
#define SplitString_h 1

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include "fstream"   // 


#include "IsotopeValues.hh"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;
*/

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

	map<string, isotope_values * > loadup_textfile(string paramfilename)
	// puts the contents of the text file into 'theInputs'.  returns theInputs.
	{
		bool verbose = false;
		string isotopeName;
		map<string, isotope_values * > theInputs;
		
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


}




//#endif
