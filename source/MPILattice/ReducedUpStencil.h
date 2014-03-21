#ifndef REDUCEDUPSTENCIL_H
#define REDUCEDUPSTENCIL_H
#include "Site.h"
#include <vector>

namespace Lattice {

class ReducedUpStencil {
	public:
		//The neighbour sites that should be used
		static std::vector<Site> neighbourSites;

		static const int id;
		
		static void initializeNeighbourSites() {
			neighbourSites.clear();
			neighbourSites.push_back(Site(1,0,0,0));
			neighbourSites.push_back(Site(0,1,0,0));
			neighbourSites.push_back(Site(0,0,1,0));
			neighbourSites.push_back(Site(0,0,0,1));
		}
};

}

#endif
