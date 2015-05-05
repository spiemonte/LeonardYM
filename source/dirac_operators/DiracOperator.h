/*
 * DiracOperator.h
 *
 *  Created on: Mar 28, 2012
 *      Author: spiem_01
 */

#ifndef DIRACOPERATOR_H_
#define DIRACOPERATOR_H_

#include "Environment.h"
#include "hmc_forces/FermionForce.h"

#include <string>

namespace Update {

class DiracOperator {
public:
	DiracOperator();
	DiracOperator(const extended_fermion_lattice_t& _lattice, double _kappa = 0., bool _gamma5 = true);
    virtual ~DiracOperator();
    static DiracOperator* getInstance(const std::string& name, unsigned int power, const StorageParameters& parameters);
	
#ifdef ENABLE_MPI
	void multiply(standard_dirac_vector_t& output, const standard_dirac_vector_t& input) {
		reduced_dirac_vector_t _output, _input = input;
		this->multiply(_output,_input);
		output = _output;
	}

	void multiply(extended_dirac_vector_t& output, const extended_dirac_vector_t& input) {
		reduced_dirac_vector_t _output, _input = input;
		this->multiply(_output,_input);
		output = _output;
	}
#endif
	
    /**
	 * This routine multiplies the Dirac operator to input and stores the result in output
	 * @param output
	 * @param input
	 */
    virtual void multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) = 0;
    /**
	 * This routine multiplies the Dirac operator to vector1 and stores the result in output adding alpha*vector2
	 * @param output
	 * @param vector1
	 * @param vector2
	 * @param alpha
	 */
    virtual void multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t & vector2, const complex& alpha) =0;

    virtual FermionForce* getForce() const = 0;

    /**
     * This function sets the kappa parameter of the operator
     * @param _kappa
     */
    virtual void setKappa(real_t _kappa);

    /**
     * This function returns back the kappa parameter of the operator
     * @return kappa
     */
    real_t getKappa() const;

    /**
     * This function returns back the pointer to the lattice
     * @return the pointer to the lattice
     */
    const reduced_fermion_lattice_t *getLattice() const;

    /**
     * This function set the lattice
     * @param lattice
     */
    virtual void setLattice(const extended_fermion_lattice_t& _lattice);

	void setGamma5(bool _gamma5);
	bool getGamma5() const;

protected:
    reduced_fermion_lattice_t lattice;

    real_t kappa;

	bool gamma5;
};

} /* namespace Update */
#endif /* DIRACOPERATOR_H_ */
