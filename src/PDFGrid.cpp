//LHAPDF Includes
#include "PDFGrid.h"
#include "Interpolator.h"
#include "factory.h"

//Standard Includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

//YAML-CPP Includes
#include "yaml-cpp/yaml.h"

namespace LHAPDF {
	bool PDFGrid::hasInterpolator() const {
		if (interpolator != 0)
			return true;
		else 
			return false;
	}
	
	void PDFGrid::setInterpolator( Interpolator* i ) {
		if (interpolator != 0) {
			if (allocatedInterpolator)
				//Delete the allocated interpolator
				delete interpolator;
			else
				//Release association with PDFGrid
				interpolator->bind( *this );
		}
		
		//Set Interpolator
		interpolator = i;
		allocatedInterpolator = false;
		
		//Allow preprocessing of PDFGrid
		interpolator->bind( *this );
	}
	
	void PDFGrid::setDefaultInterpolator() {
		if (interpolator != 0) {
			if (allocatedInterpolator)
				//Delete the allocated interpolator
				delete interpolator;
			else
				//Release association with PDFGrid
				interpolator->unbind( *this );
		}
		
		std::cout << "Ejected old interpolator" << std::endl;
		
		//Use Interpolator factory
		std::string name = set.getMetadata("SetInterpolator");
		interpolator = createInterpolator( name );
		allocatedInterpolator = true;
		
		std::cout << "Loaded new interplator" << std::endl;
	}
	
	bool PDFGrid::hasExtrapolator() const {
		if (extrapolator != 0)
			return true;
		else 
			return false;
	}
	
	void PDFGrid::setExtrapolator( Extrapolator* e ) {
		if (extrapolator != 0) {
			if (allocatedExtrapolator)
				delete extrapolator;
			else
				extrapolator->unbind( *this );
		}
		allocatedExtrapolator = false;
		
		extrapolator = e;
		extrapolator->bind( *this );
	}
	
	void PDFGrid::setDefaultExtrapolator() {
		if (extrapolator != 0) {
			if (allocatedExtrapolator)
				delete extrapolator;
			else
				extrapolator->unbind( *this );
		}
	
		//Use Extrapolator factory
		extrapolator = createExtrapolator( set.getMetadata("SetExtrapolator") );
		allocatedExtrapolator = true;
	}
	
	double PDFGrid::xfxQ2( const PID_t id, const X_t x, const Q2_t q2 ) const
	{
		///Error checks
		//Physical Ranges
		if (!inPhysicalRangeQ2( x, q2 )) {
			std::stringstream error;
			error << "Unphysical range given: " << x << " " << q2;
			
			throw std::runtime_error( error.str() );
		}
		
		//Undefined PIDs
		if (flavors.find( id ) == flavors.end()) {
			std::stringstream error;
			error << "Undefined flavour requested: " << id;
			
			throw std::runtime_error( error.str() );
		}
		
		//OUTSIDE IS EXTRAPOLATION : INSIDE IS INTERPOLATION
		if (inRangeQ2( x, q2 )) {
			//Valid Interpolator
			if( !hasInterpolator()) {
				std::stringstream error;
				error << "Undefined interpolator.";
			
				throw std::runtime_error( error.str() );
			}	
			
			//Invoke interpolator to find value
			return interpolator->interpolateQ2( *this, id, x, q2 );
		}
		else {
			//Valid Extrapolator
			if (!hasExtrapolator()) {
				std::stringstream error;
				error << "Undefined extrapolator.";
			
				throw std::runtime_error( error.str() );
			}
		
			//Invoke extrapolator to find value
			return extrapolator->extrapolateQ2( *this, id, x, q2 );
		}
	}
	
	Idx_t bisectionSearch( const double x, const AxisKnots& knots ) {
		Idx_t low = 0, high = knots.size()-1;
		
		while ((high - low) > 1) {
			Idx_t middle = low + ( (high - low) >> 1 );
			
			if (x > knots[middle])
				low = middle;
			else
				high = middle;
		}
		
		return low;
	}
	
	void PDFGrid::lookupClosestLow( const X_t x, const Q2_t q2, Idx_t& xidx, Idx_t& q2idx ) const {    
		//DO ERROR CHECKS ON X AND Q?
		//TRUST INTERPOLATORS?
		
		//Bisection search to find closest low limit
		//POSSIBLE: Implement a hunting algorithm to have a guess
		xidx = bisectionSearch( x, xknots );
		q2idx = bisectionSearch( q2, q2knots );
	}
	
	Idx_t binarySearch( const double x, const AxisKnots& points )
	{
		Idx_t low = 0, high = points.size()-1;
		
		while (true) {
			Idx_t middle = low + ((high - low) >> 2);
			
			if (middle == low || middle == high)
				break;
			
			if (x > points[middle])
				low = middle;
			else
				high = middle;
		}
		
		return low;
	}
	
	void PDFGrid::lookupClosest( const X_t x, const Q2_t q2, Idx_t& xidx, Idx_t& q2idx ) const {
		//Binary Search Axes
		xidx = binarySearch( x, xknots );
		q2idx = binarySearch( q2, q2knots );
	}
	
	inline void tokenize( const std::string& line, AxisKnots& knots ) {
		size_t start = 0;
		size_t pos = 0;
		
		do {
			pos = line.find_first_of( " ", start );
			
			knots.push_back( atof( line.substr( start, (pos - start) ).c_str() ) );
			
			start = pos+1;
		} while (pos != std::string::npos);
	}
	
	PDFGrid* PDFGrid::load( const std::string& path, const PDFSet& set ) {
		//File stream to member
		std::cout << "Loading member " << path << std::endl;
		
		std::ifstream file( path.c_str() );
		//file.exceptions( std::ifstream::failbit | std::ifstream::badbit | std::ifstream::eofbit );
		
		//Initiate PDFGrid
		std::cout << "Initialize set" << std::endl;
		
		PDFGrid* grid = new PDFGrid( set );
		
		//Read header
		std::string line;
		std::stringstream header;
		header << "---\n";
		
		//Get next line
		while (true) {
			getline( file, line );
									
			if (line.empty()) {
				break;
			}
			else if (line.size()==3) {
				//FIND A BETTER DELIM METHOD!
				break;
			}
			
			header << line << "\n";
		}
						
		//Parse header
		YAML::Parser parser( header );
		YAML::Node headdoc;
		parser.GetNextDocument( headdoc );
		
		headdoc["MemberName"] >> grid->name;
		headdoc["MemberID"] >> grid->id;
		
		std::cout << "YAML done..." << std::endl;
		
		//Parse Xs
		for (YAML::Iterator xsit = headdoc["Xs"].begin(); xsit != headdoc["Xs"].end(); ++xsit) {
			double x;
			(*xsit) >> x;
			
			grid->xknots.push_back(x);
		}
		
		//Parse Q2s
		for (YAML::Iterator q2sit = headdoc["Q2s"].begin(); q2sit != headdoc["Q2s"].end(); ++q2sit) {
			double q2;
			(*q2sit) >> q2;
			
			grid->q2knots.push_back(q2);
		}
		
		//Check for valid file
		if ( !file.good() ) {
			std::cout << "File not good" << std::endl;
			throw -1;
		}
		
		std::cout << "Load data from " << path << std::endl;
		
		//Allocate flavor data
		std::vector<PID_t>::const_iterator piditer = set.getFlavours().begin();
		for( ; piditer != set.getFlavours().end(); ++piditer ) {
			//std::cout << *piditer << std::endl;
			grid->flavors[*piditer] = new double[grid->xknots.size()*grid->q2knots.size()];
		}
		
		double* d = grid->flavors[1];
		
		//Parse flavor data
		/*uint32_t lcount = 0;
		while (file.good()) {
			getline( file, line );
			
			if( !file.good() ) {
				break;
			}
			
			size_t start = 0;
			size_t pos;
			std::vector<PID_t>::const_iterator piditer = set.getFlavours().begin();
			for (; piditer != set.getFlavours().end(); ++piditer) {
				pos = line.find_first_of( " ", start );
				
				grid->flavors[*piditer][lcount] = atof( line.substr( start, (pos - start) ).c_str() );
				//std::cout << line.substr( start, (pos - start) ).c_str() << " : " << grid->flavors[*piditer][lcount] << std::endl;
				
				start = pos+1;
			}
			
			++lcount;
		}*/
		
		uint32_t cline = 0;
		for(; cline < grid->xknots.size()*grid->q2knots.size(); ++cline ) {
			//GET LINE			
			if( !file.good() ) {
				std::stringstream error;
				error << cline;
				
				std::cout << "needed: " << grid->xknots.size()*grid->q2knots.size() << std::endl;
				std::cout << "received: " << cline << std::endl;
				
				throw std::runtime_error( error.str() );
			}
			getline( file, line );
			
			//PARSE LINE
			size_t start = 0;
			size_t pos;
			std::vector<PID_t>::const_iterator piditer = set.getFlavours().begin();
			for (; piditer != set.getFlavours().end(); ++piditer) {
				pos = line.find_first_of( " ", start );
				
				grid->flavors[*piditer][cline] = atof( line.substr( start, (pos - start) ).c_str() );
				//std::cout << line.substr( start, (pos - start) ).c_str() << " : " << grid->flavors[*piditer][lcount] << std::endl;
				
				start = pos+1;
			}
		}
				
		//std::cout << "Load inter/extra/polators" << std::endl;
 			
		//Set defaults
		grid->setDefaultInterpolator();
		grid->setDefaultExtrapolator();
		
		//std::cout << "Loaded inter/extra/polators" << std::endl;
		
		return grid;
	}
}