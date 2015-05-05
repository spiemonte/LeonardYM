/*
 * DeflationInverter.h
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#ifndef DEFLATIONINVERTER_H_
#define DEFLATIONINVERTER_H_
#include "./dirac_operators/DiracOperator.h"
#include "ConjugateGradient.h"

namespace Update {

#ifdef ENABLE_MPI
typedef Lattice::Lattice<short int, Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_index_lattice_t;
#endif
#ifndef ENABLE_MPI
typedef Lattice::Lattice<short int, Lattice::LocalLayout > reduced_index_lattice_t;
#endif

class LeftProjector;
class RightProjector;
class InverseLittleOperator;

class DeflationInverter {
public:
	DeflationInverter(int _recursion = 1);
	DeflationInverter(const DeflationInverter& snd);
	~DeflationInverter();

	bool solve(DiracOperator* dirac, const extended_dirac_vector_t& source, extended_dirac_vector_t& solution);

	void setPrecision(double _epsilon);
	double getPrecision() const;

	unsigned int getLastSteps() const;

	void setBasisDimension(int _basisDimension);
	int getBasisDimension() const;

	void setBlockSize(int _blockSize);
	int getBlockVolume() const;
	void setBlockSizeX(int _blockSizeX);
	void setBlockSizeY(int _blockSizeY);
	void setBlockSizeZ(int _blockSizeZ);
	void setBlockSizeT(int _blockSizeT);

	void setRecursion(int _recursion);

	void generateBasis(DiracOperator* dirac, LeftProjector* upperLeftProjector = 0);
private:
	reduced_dirac_vector_t* localBasis;
	ConjugateGradient* conjugateGradient;
	DeflationInverter* subDeflationInverter;

	LeftProjector* leftProjector;
	RightProjector* rightProjector;
	InverseLittleOperator* inverseLittleOperator;
	

	unsigned int lastStep;
	int basisDimension;
	int totalNumberOfVectors;//TODO
	int totalNumberOfBlocks;
	int blockSizeX;
	int blockSizeY;
	int blockSizeZ;
	int blockSizeT;
	int recursion;
	double precision;

	reduced_index_lattice_t blockIndex;
};

} /* namespace Update */
#endif /* DEFLATIONINVERTER_H_ */
