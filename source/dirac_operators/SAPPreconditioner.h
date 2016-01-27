/*
 * ComplementBlockDiracWilsonOperator.h
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#ifndef SAPPRECONDITIONER_H_
#define SAPPRECONDITIONER_H_
#include "ComplementBlockDiracOperator.h"

namespace Update {

class SAPPreconditioner : public DiracOperator {
public:
	SAPPreconditioner(DiracOperator* _diracOperator, ComplementBlockDiracOperator* _K);
	//ComplementBlockDiracOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa = 0.);
	~SAPPreconditioner();

	virtual void multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input);

	virtual void multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const std::complex<real_t>& alpha);

	virtual FermionForce* getForce() const;
	
	virtual void setKappa(real_t _kappa);
	
	virtual void setLattice(const extended_fermion_lattice_t& _lattice);
	
	void setSteps(int _steps);

	void setPrecision(double precision);
private:
	DiracOperator* diracOperator;
	ComplementBlockDiracOperator* K;
	reduced_dirac_vector_t tmp1;
	reduced_dirac_vector_t tmp2;
	reduced_dirac_vector_t tmp3;
	
	int steps;
	real_t precision;
};

} /* namespace Update */
#endif /* RANDOMDIRACWILSONOPERATOR_H_ */
