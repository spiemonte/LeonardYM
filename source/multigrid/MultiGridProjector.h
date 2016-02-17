#ifndef MULTIGRIDPROJECTOR_H
#define MULTIGRIDPROJECTOR_H
#include "MultiGridVector.h"

namespace Update {

class MultiGridProjector {
public:
	MultiGridProjector();

	void apply(multigrid_vector_t& output, const reduced_dirac_vector_t& input);

	void apply(reduced_dirac_vector_t& output, const multigrid_vector_t& input);

	//Add a vector to the vector space used for block deflation
	void addVector(const reduced_dirac_vector_t& vector) {
		vectorspace.push_back(&vector);
	}

	//Remove all the vectors in the vector space
	void clearVectorSpace() {
		vectorspace.clear();
	}

protected:
	std::vector<reduced_dirac_vector_t const*> vectorspace;

	reduced_dirac_vector_t tmp_input;
	reduced_dirac_vector_t tmp_output;
};

}

#endif

