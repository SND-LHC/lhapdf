#include "LHAPDF/PDFSet.h"
#include "LHAPDF/GridPDF.h"
#include "LHAPDF/AlphaS.h"

#include <sstream>
#include <fstream>
#include <stdexcept>

#include "yaml-cpp/yaml.h"

namespace LHAPDF {


  PDFSet::~PDFSet() {
    // Delete members
    for (std::map<size_t, PDF*>::iterator member = _members.begin(); member != +members.end(); ++member) {
      delete member->second;
    }
    _members.clear();
  }


  PDF& PDFSet::getMember( const size_t id ) {
    std::map<size_t, PDF*>::const_iterator member = members.find(id);
    if (member == members.end()) {
      // Load member
      std::stringstream mpath;
      mpath << path << "mbr_" << id << ".LHm";
      members[id] = PDF::load( this, mpath.str() );
      return *(members[id]);
    } else {
      // Member is already loaded
      return *(member->second);
    }
  }


  PDFSet* PDFSet::load(const std::string& path) {
    /// @todo Put .pdfgrid on the end automatically when loading members (not here, though!)

    // Create empty set
    PDFSet* set = new PDFSet(path);

    // Generate path to meta file
    std::stringstream info;
    info << path << "meta.LHinfo";
    std::cout << info.str() << std::endl;

    // Open meta file
    std::ifstream file( info.str().c_str(), std::ifstream::in );

    if (!file.good()) {
      // Clean up loose memory
      delete set;
      // Give back human error
      std::stringstream error;
      error << "Failed to open: " << info.str();
      throw ReadError( error.str() );
    }

    std::cout << "YAML start set" << std::endl;
    // Start parsing with YAML
    YAML::Parser parser( file );
    YAML::Node metadoc;

    /// @todo THIS CAN FAIL IF THERE ARE YAML SYNTAX ERRORS
    parser.GetNextDocument( metadoc );

    // Read all key value pairs
    for (YAML::Iterator it = metadoc.begin(); it != metadoc.end(); ++it) {
      std::string key;
      it.first() >> key;

      if (key == "SetFlavours") {
        std::cout << "FLAVOUR" << std::endl;

        // Treat flavours differently
        try {
          for (YAML::Iterator fit = it.second().begin(); fit != it.second().end(); ++fit) {
            int flavour;
            fit->Read( flavour );
            /// @todo ERROR: THIS MIGHT CAUSE SOME SORT OF OVERFLOW ERROR
            // When testing with negative numbers for antiparticles pids, numbers seemed to be treated wrongly?
            set->flavours.push_back( flavour );
          }
        } catch(std::exception& e) {
          /// @todo Handle error properly
          std::cout << "EXCEPT: " <<  e.what() << std::endl;
        }
      }
      else if (key == "SetReferences") {
        // Treat references differently?
        /// @todo bla bla bla
      }
      else {
        // Simple key value pair
        std::string value;
        it.second() >> value;
        set->meta[key] = value;
      }
    }

    return set;
  }


  PDFSet* PDFSet::loadByName(const std::string& name) {
    // Search default path for name
    /// @todo This needs to be changed to a portable method
    std::stringstream path;
    //LHAPDFPATH/ << NAME
    //path << "/usr/local/share/lhapdf/PDFsets/";
    /// @todo WOW. Talk about hard-coding! Repace with binreloc tool
    path << "/phys/linux/s0821167/Summer Project/Interpolator/lhapdfv6/tests/";
    path << name << "/";
    return load(path.str());
  }


  /// @todo Need to be able to use both numerical integration and interpolation to get alphaS
  double PDFSet::alphaS(double Q2) const {
    /// @todo Also need to know quark masses, number of flavours, order of running, etc.
    return AlphaS::numerical( 0.12, Q2 );
  }


}
