#ifndef SOLVER_H
#define SOLVER_H
#include "Environment.h"
#include <string>

namespace Update {

class Solver {
public:
	Solver(const std::string& _name = "") : name(_name), precision(0.0000000001), maxSteps(1000) { }
	virtual ~Solver() { }

	virtual bool solve(DiracOperator* , const reduced_dirac_vector_t& , reduced_dirac_vector_t&  , reduced_dirac_vector_t const* = 0) {
		if (isOutputProcess()) std::cout << "Solver not implemented by this class " << name << std::endl;
		return false;
	}

#ifdef ENABLE_MPI
	virtual bool solve(DiracOperator* dirac, const extended_dirac_vector_t& source, extended_dirac_vector_t& solution, extended_dirac_vector_t const* initial_guess = 0) {
		reduced_dirac_vector_t red_source = source;
		reduced_dirac_vector_t red_solution;
		bool res = false;
		if (initial_guess != 0) {
			reduced_dirac_vector_t red_initial_guess = *initial_guess;
			res = this->solve(dirac, red_source, red_solution, &red_initial_guess);
		}
		else {
			res = this->solve(dirac, red_source, red_solution, 0);
		}
		solution = red_solution;
		return res;
	}
#endif

	void setPrecision(real_t _epsilon) {
		precision = _epsilon;
	}
	
	real_t getPrecision() const {
		return precision;
	}

	void setMaximumSteps(unsigned int _maxSteps) {
		maxSteps = _maxSteps;
	}
	
	unsigned int getMaximumSteps() const {
		return maxSteps;
	}

	real_t getLastError() const {
		return lastError;
	}
	
	unsigned int getLastSteps() const {
		return lastSteps;
	}

protected:
	std::string name;

	real_t precision;
	unsigned int maxSteps;
	real_t lastError;
	unsigned int lastSteps;

};

}


#endif

