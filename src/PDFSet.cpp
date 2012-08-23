#include "PDFSet.h"
#include <math.h>

//Standard Includes
#include <sstream>
#include <fstream>
#include <stdexcept>

//LHAPDF Includes
#include "PDFGrid.h"

//Constants
static const double PI = 3.1415926535;

// Masses of quarks in GeV.  Used to calculate the number of quarks that are active at a given energy range Q2
static const double quarkMasses[] = {0.0017, 0.0041, 0.1, 1.29, 4.1, 172};

//TODO: maybe allow for general lambdaQCD, set to 0.4 for now
static double lambdaQCD = 0.4;

//YAML-CPP Includes
#include "yaml-cpp/yaml.h"

namespace LHAPDF {
	PDFSet::~PDFSet() {
		//Delete members
		for (std::map<Member_t, PDF*>::iterator member = members.begin(); member != members.end(); ++member) {
			delete (member->second);
		}
	}
	
	PDF& PDFSet::getMember( const Member_t id ) {
		std::map<Member_t, PDF*>::const_iterator member = members.find(id);
		
		if (member == members.end()) {
			//Load member
			std::stringstream mpath;
			mpath << path << "mbr_" << id << ".LHm";
			
			members[id] = PDFGrid::load( mpath.str(), *this );
			
			return *(members[id]);
		}
		else {
			//Member is already loaded
			return *(member->second);
		}
	}

	PDFSet* PDFSet::load( const std::string& path ) {
		//Create empty set
		PDFSet* set = new PDFSet( path );
		
		//Generate path to meta file
		std::stringstream info;
		info << path << "/meta.LHinfo";
							
		//Open meta file
		std::ifstream file( info.str().c_str(), std::ifstream::in );
		
		if( !file.good() ) {
			//Clean up loose memory
			delete set;
			
			//Give back human error
			std::stringstream error;
			error << "Failed to open: " << info.str();
			
			throw std::runtime_error( error.str() );
		}
		
		std::cout << "YAML start set" << std::endl;
		//Start parsing with YAML
		YAML::Parser parser( file );
		YAML::Node metadoc;

		//THIS CAN FAIL IF THERE ARE SYNTAX ERRORS
		parser.GetNextDocument( metadoc );
		
		//Read all key value pairs
		for( YAML::Iterator it = metadoc.begin(); it != metadoc.end(); ++it ) {
			std::string key;
			it.first() >> key;
			
			if (key == "SetFlavours") {
				std::cout << "FLAVOUR" << std::endl;
				
				//Treat flavours differently
				try {
				for (YAML::Iterator fit = it.second().begin(); fit != it.second().end(); ++fit) { 
					PID_t flavour;
					fit->Read( flavour );
					
					//ERROR: THIS MIGHT CAUSE SOME SORT OF OVERFLOW ERROR, when testing when negative
					// numbers were given as anti particles pid it seemed to do something with a 16 bit
					// number not a native int... (64bit), however it works currently and I don't feel
					// like breaking it again...
										
					set->flavours.push_back( flavour );
				} }
				catch( std::exception e )
				{
					std::cout << "EXCEPT: " <<  e.what() << std::endl;
					}
			}
			else if(key == "SetReferences") {
				//Treat references differently?
				//bla bla bla
			}
			else {
				//Simple key value pair
				std::string value;
				it.second() >> value;
								
				set->meta[key] = value;
			}
		}
				
		return set;
	}
	
	PDFSet* PDFSet::loadByName( const std::string& name ) {
		//Search default path for name
		
		std::stringstream path;
		//LHAPDFPATH/ << NAME
		path << "/usr/local/share/lhapdf/PDFsets/";
		path << name << "/";
		
		//Use load by path to load file
		return load( path.str() );
	}

	//ALPHAS Calculation
	double PDFSet::alphaS( double q2 ) const {
		// Find the number of quarks with mass thresholds less than given energy scale Q2
		bool isNf = false;
		int nf = 0;
		while ( isNf == false ) {
		if ( quarkMasses[nf] <= q2 && nf < 6) {
			nf++;
		} else {
			isNf = true;
		}
		}
	
		// Calculate beta-functions
		double beta0 = 11 - ((2 / 3.0) * nf);
		double beta1 = 51 - ((19 / 3.0) * nf);
		double beta2 = 2857 - ((5033 / 9.0) * nf) + ((325 / 27.0) * nf * nf);
	
		// Pre-calculate powers of these 
		double beta02 = beta0 * beta0;
		double beta12 = beta1 * beta1;
	
		// Pre-calculate powers of ln(Q2/lambdaQCD) 
		double x = q2 / lambdaQCD;
		double lnx = log(x);
		double lnlnx = log(lnx);
		
		// Pre-calculate powers of expansion term y, where y is 1/ln(Q2/lambdaQCD)
		double y = 1 / lnx;
		double y2 = y * y;
	
		// Calculate expansion terms Ay[a_0*y^0 - a_1*y^1 + By^2[a_20 + a_21 - a_22]] where a_ij are series coefficients
		double A = 4 * PI / beta0;
		double a_0 = 1;
		double a_1 = -2 * beta1 * lnlnx / beta02;
		double B = 4 * beta12 / (beta02 * beta02);
		double a_20 = (lnlnx - 0.5) * (lnlnx - 0.5);
		double a_21 = beta2 * beta0 / (8 * beta12);
		double a_22 = 5 / 4.0;
	
		// Combine terms
		double alphaS = A * y * (a_0 - (a_1*y) + (B * y2 * (a_20 + a_21 - a_22)));
		
		return alphaS;
	}
}
