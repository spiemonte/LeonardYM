#ifndef MMMRMULTISHIFTSOLVER_H_
#define MMMRMULTISHIFTSOLVER_H_
#include "MultishiftSolver.h"

namespace Update {

class MMMRMultishiftSolver : public MultishiftSolver {
public:
	MMMRMultishiftSolver(real_t _epsilon = 0.00000001, unsigned int _maxSteps = 3000);
	virtual ~MMMRMultishiftSolver();

	virtual bool solve(DiracOperator* dirac, const extended_dirac_vector_t& source, std::vector<extended_dirac_vector_t>& solutions, const std::vector<real_t>& shifts);

private:
#ifndef ALIGNED_OPT
	//Temporary vectors
	reduced_dirac_vector_t r;
	reduced_dirac_vector_t p;
#else
	//Temporary vectors
	reduced_soa_dirac_vector_t r;
	reduced_soa_dirac_vector_t p;
#endif

	//Omega is the overrelaxation parameter
	real_t omega;
};

} /* namespace Update */
#endif /* MMMRMULTISHIFTSOLVER_H_ */
