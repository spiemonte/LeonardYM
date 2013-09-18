#ifndef LATTICECHUNK_H
#define LATTICECHUNK_H

namespace Lattice {

class LatticeChunk {
	public:
		int id;
		int owner;
		int size;
		int offset;
		
		std::vector<int> sharers;
		std::vector<int> tags;
};

}

#endif