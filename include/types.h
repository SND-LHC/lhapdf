#pragma once

#include <stdint.h>
#include <cstdlib>
#include <vector>
#include <map>

namespace LHAPDF {
	///
	typedef int PID_t;
	
	///
	typedef double* PIDdata;
	
	///
	typedef std::map<PID_t,PIDdata> PIDmap;
	
	///
	typedef std::vector<double> AxisKnots;
	
	///
	typedef double X_t;
	
	///
	typedef double Q_t;
	
	///
	typedef double Q2_t;
	
	///
	typedef size_t Member_t;
	
	///
	typedef size_t Idx_t;
}