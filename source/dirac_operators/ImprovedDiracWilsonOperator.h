/*
 * ImprovedDiracWilsonOperator.h
 *
 *  Created on: May 4, 2012
 *      Author: spiem_01
 */

#ifndef IMPROVEDDIRACWILSONOPERATOR_H_
#define IMPROVEDDIRACWILSONOPERATOR_H_

#include "DiracOperator.h"

namespace Update {

class ImprovedDiracWilsonOperator : public DiracOperator {
public:
	ImprovedDiracWilsonOperator();
	ImprovedDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, double _kappa = 0., double _csw = 1., bool _gamma5 = true);
	~ImprovedDiracWilsonOperator();

	/**
	 * This routine multiplies the Dirac operator to input and stores the result in output
	 * @param output
	 * @param input
	 */
	virtual void multiply(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & input);
	/**
	 * This routine multiplies the Dirac operator to vector1 and stores the result in output adding alpha*vector2
	 * @param output
	 * @param vector1
	 * @param vector2
	 * @param alpha
	 */
	virtual void multiplyAdd(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & vector1, const reduced_dirac_vector_t & vector2, const complex& alpha);

	virtual void setLattice(const extended_fermion_lattice_t& _lattice);

	virtual FermionForce* getForce() const;

	real_t getCSW() const;
	void setCSW(real_t _csw);

protected:
	//The clover term
	real_t csw;
	//The field strength
	//FermionicGroup (* F)[6];//TODO
	reduced_field_strength_lattice_t F;

	void updateFieldStrength(const extended_fermion_lattice_t& _lattice);
};

} /* namespace Update */
#endif /* IMPROVEDDIRACWILSONOPERATOR_H_ */
