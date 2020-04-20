#ifndef OVERLAPCHIRALROTATION_H_
#define OVERLAPCHIRALROTATION_H_

#include "LatticeSweep.h"
#include "fermion_measurements/StochasticEstimator.h"
#include "dirac_operators/DiracOperator.h"
#include "inverters/Solver.h"

namespace Update {

class OverlapChiralRotation : public Update::LatticeSweep, Update::StochasticEstimator {
public:
	OverlapChiralRotation();
	OverlapChiralRotation(const OverlapChiralRotation& toCopy);
	~OverlapChiralRotation();

	virtual void execute(environment_t& environment);

	static void registerParameters(po::options_description& desc);

	void chiral_rotation_psibar(DiracOperator* dirac, const reduced_dirac_vector_t& input, reduced_dirac_vector_t& output, unsigned int steps, const std::complex<real_t>& theta) const;
	void chiral_rotation_psi(DiracOperator* dirac, const reduced_dirac_vector_t& input, reduced_dirac_vector_t& output, unsigned int steps, const std::complex<real_t>& theta) const;

private:
	reduced_dirac_vector_t randomNoise;
	reduced_dirac_vector_t tmp;
	reduced_dirac_vector_t tmp1;
	reduced_dirac_vector_t tmp2;
	reduced_dirac_vector_t tmp3;
	reduced_dirac_vector_t inverse;
	reduced_dirac_vector_t barpsi_rotated;
	reduced_dirac_vector_t psi_rotated;
	DiracOperator* diracOperator;
	Solver* inverter;
};

} /* namespace Update */

#endif /* OverlapChiralRotation_H_ */
