#ifndef BLOCKBASIS_H
#define BLOCKBASIS_H
#include "Environment.h"

namespace Update {

class BlockBasis {
public:
	BlockBasis(int _basisDimension);
	BlockBasis(const BlockBasis& toCopy);
	~BlockBasis();

	void orthogonalize();

	void orthogonalize(int index);

	reduced_dirac_vector_t& operator[](int index) {
		return vectorspace[index];
	}

	const reduced_dirac_vector_t& operator[](int index) const {
		return vectorspace[index];
	}

	int size() const {
		return basisDimension;
	}

	void setBasisDimension(int _basisDimension);
	
private:
	int basisDimension;
	reduced_dirac_vector_t* vectorspace;
};

}

#endif

