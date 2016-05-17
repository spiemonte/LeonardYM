/*
 * ImprovedDiracWilsonOperator.h
 *
 *  Created on: May 4, 2012
 *      Author: spiem_01
 */

#ifndef EVENODDIMPROVEDDIRACWILSONOPERATOR_H_
#define EVENODDIMPROVEDDIRACWILSONOPERATOR_H_

#include "ImprovedDiracWilsonOperator.h"
#include "MatrixTypedef.h"

typedef Eigen::Matrix< std::complex<Update::real_t> , 4*Update::diracVectorLength, 4*Update::diracVectorLength > clover_matrix_t;
enum Part {EVEN = 0, ODD};

namespace Update {

class EvenOddImprovedDiracWilsonOperator : public ImprovedDiracWilsonOperator {
public:
	EvenOddImprovedDiracWilsonOperator();
	EvenOddImprovedDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, double _kappa = 0., double _csw = 1., bool _gamma5 = true);
	~EvenOddImprovedDiracWilsonOperator();

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

public:
	void multiplyOddOdd(reduced_dirac_vector_t & output, Part part);
	void multiplyOddOddMinusIdentity(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & input, Part part);
	void multiplyEvenEvenInverse(reduced_dirac_vector_t & output);
	void multiplyEvenOdd(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & input, Part part);

private:
	//The field strength
	clover_matrix_t *cloverMatrixInverse;//TODO

	void calculateInverseEvenEven();

	
};

} /* namespace Update */
#endif /* IMPROVEDDIRACWILSONOPERATOR_H_ */
