#ifndef MULTIGRIDOPERATOR_H
#define MULTIGRIDOPERATOR_H
#include "MultiGridVector.h"
#include "dirac_operators/DiracOperator.h"

namespace Update {

class MultiGridOperator {
public:
	MultiGridOperator();

	void multiply(multigrid_vector_t& output, const multigrid_vector_t& input);

	//The add operation is applied only to the inner Dirac operator
	void multiplyAdd(multigrid_vector_t& output, const multigrid_vector_t& input, const complex& alpha);

	//To set the Dirac operator
	void setDiracOperator(DiracOperator* _dirac) {
		dirac = _dirac;
	}

	//Add a vector to the vector space used for block deflation
	void addVector(const reduced_dirac_vector_t& vector) {
		vectorspace.push_back(&vector);
	}

	//Remove all the vectors in the vector space
	void clearVectorSpace() {
		vectorspace.clear();
	}

	matrix_t asMatrix();

protected:
	std::vector<reduced_dirac_vector_t const*> vectorspace;
	DiracOperator* dirac;

	reduced_dirac_vector_t tmp_input;
	reduced_dirac_vector_t tmp_output;
};

}

#endif

