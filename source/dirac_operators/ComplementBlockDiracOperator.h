#ifndef COMPLEMENTBLOCKDIRACOPERATOR_H_
#define COMPLEMENTBLOCKDIRACOPERATOR_H_
#include "BlockDiracWilsonOperator.h"
#include "BlockDiracOperator.h"
#include "inverters/BiConjugateGradient.h"
#include "inverters/GMRESR.h"

namespace Update {

class ComplementBlockDiracOperator : public BlockDiracOperator {
public:
	ComplementBlockDiracOperator(DiracOperator* _diracOperator, BlockDiracOperator* _redBlockDiracOperator, BlockDiracOperator* _blackBlockDiracOperator);
	ComplementBlockDiracOperator(const ComplementBlockDiracOperator& toCopy);
	//ComplementBlockDiracOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa = 0.);
	~ComplementBlockDiracOperator();

	virtual void multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input);

	virtual void multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const std::complex<real_t>& alpha);

	virtual FermionForce* getForce() const;
	
	virtual void setKappa(real_t _kappa);
	
	virtual void setLattice(const extended_fermion_lattice_t& _lattice);

	void setPrecision(double precision);
	void setMaximumSteps(int steps);
private:
	DiracOperator* diracOperator;
	BlockDiracOperator* redBlockDiracOperator;
	BlockDiracOperator* blackBlockDiracOperator;
	//BiConjugateGradient* biConjugateGradient;
	GMRESR* gmresr;
	reduced_dirac_vector_t tmp1;
	reduced_dirac_vector_t tmp2;
};

} /* namespace Update */
#endif /* RANDOMDIRACWILSONOPERATOR_H_ */
