/*
 * DeflationInverter.h
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#ifndef DEFLATIONINVERTER_H_
#define DEFLATIONINVERTER_H_
#include "./dirac_operators/DiracOperator.h"
#include "BiConjugateGradient.h"

namespace Update {

class DeflationInverter {
public:
	DeflationInverter();
	DeflationInverter(const DeflationInverter& snd);
	~DeflationInverter();

	bool solve(DiracOperator* dirac, const extended_dirac_vector_t& source, extended_dirac_vector_t& solution);

	void setPrecision(double _epsilon);
	double getPrecision() const;

	unsigned int getLastSteps() const;

	void setBasisDimension(int _basisDimension);
	int getBasisDimension() const;

	void setBlockDivision(int _blockDivision);
	int getBlockDivision() const;
private:
	reduced_dirac_vector_t* localBasis;
	BiConjugateGradient* biConjugateGradient;

	unsigned int lastStep;
	int basisDimension;
	int blockDivision;
	int totalNumberOfVectors;
	int totalBlocks;
	double precision;

	void generateBasis(DiracOperator* dirac);
};

} /* namespace Update */
#endif /* DEFLATIONINVERTER_H_ */
