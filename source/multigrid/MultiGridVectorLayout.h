#ifndef MULTIGRIDVECTORLAYOUT_H
#define MULTIGRIDVECTORLAYOUT_H
#include "Environment.h"

namespace Update {

class MultiGridVectorLayout {
public:
	static void initialize();
	static void setBasisDimension(int _basisDimension);
	
	static unsigned int xBlockSize;
	static unsigned int yBlockSize;
	static unsigned int zBlockSize;
	static unsigned int tBlockSize;

	static int totalNumberOfBlocks;
	static int basisDimension;
	static int size;

	static int index(int site) {
		return (*blockIndex)[site];
	}

	static unsigned int instances;
private:
	
	static reduced_index_lattice_t* blockIndex;
};

}

#endif

