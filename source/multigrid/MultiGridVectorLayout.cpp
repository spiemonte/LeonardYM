#include "MultiGridVectorLayout.h"

namespace Update {

unsigned int MultiGridVectorLayout::xBlockSize = 4;
unsigned int MultiGridVectorLayout::yBlockSize = 4;
unsigned int MultiGridVectorLayout::zBlockSize = 4;
unsigned int MultiGridVectorLayout::tBlockSize = 4;

unsigned int MultiGridVectorLayout::instances = 0;

int MultiGridVectorLayout::totalNumberOfBlocks;
int MultiGridVectorLayout::basisDimension = 40;
int MultiGridVectorLayout::size;
	
reduced_index_lattice_t* MultiGridVectorLayout::blockIndex;

void MultiGridVectorLayout::initialize() {
	if (instances != 0 && isOutputProcess()) std::cout << "MultiGridVectorLayout::Warning, there are " << instances << " MultiGridVector classes created, they will not work properly!" << std::endl; 
	
	typedef reduced_index_lattice_t::Layout LT;
	//number of blocks in the x,y,z,t direction
	int numberBX = LT::loc_x/xBlockSize + ((LT::loc_x % xBlockSize) != 0 ? 1 : 0);
	int numberBY = LT::loc_y/yBlockSize + ((LT::loc_y % yBlockSize) != 0 ? 1 : 0);
	int numberBZ = LT::loc_z/zBlockSize + ((LT::loc_z % zBlockSize) != 0 ? 1 : 0);
	int numberBT = LT::loc_t/tBlockSize + ((LT::loc_t % tBlockSize) != 0 ? 1 : 0);

	if (LT::loc_x % xBlockSize != 0 || LT::loc_y % yBlockSize != 0 || LT::loc_z % zBlockSize != 0 || LT::loc_t % tBlockSize != 0) {
		if (isOutputProcess()) std::cout << "MultiGridOperator::Warning, block grid does not evenly match the processor grid: (" << LT::loc_x % xBlockSize << "," << LT::loc_y % yBlockSize << "," << LT::loc_z % zBlockSize << "," << LT::loc_t % tBlockSize << ")" << std::endl;
	}

	totalNumberOfBlocks = numberBX*numberBY*numberBZ*numberBT;
	size = totalNumberOfBlocks*basisDimension;

	blockIndex = new reduced_index_lattice_t;
	
	for (int site = 0; site < blockIndex->localsize; ++site) {
		int x = (LT::globalIndexX(site) % LT::loc_x)/xBlockSize;
		int y = (LT::globalIndexY(site) % LT::loc_y)/yBlockSize;
		int z = (LT::globalIndexZ(site) % LT::loc_z)/zBlockSize;
		int t = (LT::globalIndexT(site) % LT::loc_t)/tBlockSize;

		(*blockIndex)[site] = numberBT*(numberBZ*(numberBY*x + y) + z) + t;
	}
	blockIndex->updateHalo();
}

void MultiGridVectorLayout::setBasisDimension(int _basisDimension) {
	basisDimension = _basisDimension;
	size = basisDimension*totalNumberOfBlocks;
}

}


