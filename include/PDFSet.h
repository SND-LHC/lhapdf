//Include Guard
#pragma once

//Standard Includes
#include <string>

//LHAPDF Includes
#include "types.h"
#include "PDF.h"

namespace LHAPDF {
	//Type Declarations
	typedef unsigned int Standard_t;
	typedef unsigned int Order_t;
	
	//PDFSet Declaration
	/// The PDFSet is the most general container.
	
	/// A Set contains PDF 'members' and meta data that is either Set general or default meta data for
	/// members.
	class PDFSet {
	public:
		/// Deletes allocated memory when deconstructing.
		~PDFSet();
		
		/// Retrieves general metadata
		inline std::string getMetadata( const std::string& ) const;
		
		/// Specific meta data
		inline std::string getName() const;
		inline std::string getDescription() const;
		
		inline Order_t getOrder() const;
		inline double getAlphaS() const;
		
		inline const std::vector<PID_t>& getFlavours() const;
			
		/// Get selected member.
		
		/// This will load the member if it is not loaded yet or simply return a pointer to the loaded
		/// instance.
		PDF& getMember( const Member_t );
		
		/// Loads the set by path
		static PDFSet* load( const std::string& );
		
		/// Loads the set by name
		
		/// This function loads uses a default path to load a named set.
		static PDFSet* loadByName( const std::string& );

		//TEMP
		double alphaS( double q2 ) const;
		
	private:
		//Declaring unaccessible functions
		PDFSet( const std::string& path ) : path(path) {}
		PDFSet( const PDFSet& set );
		void operator=( const PDFSet& set );
		
		/// Holds path for member retrieving
		std::string path;
			
		/// Holds all loaded members
		std::map<Member_t, PDF*> members;
		
		/// Holds all meta data
		std::map<std::string, std::string> meta;
		
		/// Holds all flavours defined by this set
		std::vector<PID_t> flavours;
	};
	
	const std::vector<PID_t>& PDFSet::getFlavours() const {
		return flavours;
	}
	
	std::string PDFSet::getMetadata( const std::string& key ) const {
		std::map<std::string, std::string>::const_iterator data = meta.find( key );
		
		if( data != meta.end() )
			return data->second;
		else {
			//Metadata not found
			std::stringstream error;
			error << "Metadata for key: " << key << " not found.";
			
			throw std::runtime_error( error.str() );
		}
	}
}