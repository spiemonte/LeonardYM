#ifndef EXTENDEDSTENCIL_H
#define EXTENDEDSTENCIL_H
#include "Site.h"
#include <vector>

namespace Lattice {

class ExtendedStencil {
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
			neighbourSites.push_back(Site(-1,0,0,0));
			neighbourSites.push_back(Site(0,-1,0,0));
			neighbourSites.push_back(Site(0,0,-1,0));
			neighbourSites.push_back(Site(0,0,0,-1));

			neighbourSites.push_back(Site(2,0,0,0));
			neighbourSites.push_back(Site(0,2,0,0));
			neighbourSites.push_back(Site(0,0,2,0));
			neighbourSites.push_back(Site(0,0,0,2));
			neighbourSites.push_back(Site(-2,0,0,0));
			neighbourSites.push_back(Site(0,-2,0,0));
			neighbourSites.push_back(Site(0,0,-2,0));
			neighbourSites.push_back(Site(0,0,0,-2));

			neighbourSites.push_back(Site(1,-1,0,0));
			neighbourSites.push_back(Site(1,0,-1,0));
			neighbourSites.push_back(Site(1,0,0,-1));

			neighbourSites.push_back(Site(-1,1,0,0));
			neighbourSites.push_back(Site(0,1,-1,0));
			neighbourSites.push_back(Site(0,1,0,-1));

			neighbourSites.push_back(Site(-1,0,1,0));
			neighbourSites.push_back(Site(0,-1,1,0));
			neighbourSites.push_back(Site(0,0,1,-1));

			neighbourSites.push_back(Site(-1,0,0,1));
			neighbourSites.push_back(Site(0,-1,0,1));
			neighbourSites.push_back(Site(0,0,-1,1));

			neighbourSites.push_back(Site(1,+1,0,0));
			neighbourSites.push_back(Site(1,0,+1,0));
			neighbourSites.push_back(Site(1,0,0,+1));

			neighbourSites.push_back(Site(+1,1,0,0));
			neighbourSites.push_back(Site(0,1,+1,0));
			neighbourSites.push_back(Site(0,1,0,+1));

			neighbourSites.push_back(Site(+1,0,1,0));
			neighbourSites.push_back(Site(0,+1,1,0));
			neighbourSites.push_back(Site(0,0,1,+1));

			neighbourSites.push_back(Site(+1,0,0,1));
			neighbourSites.push_back(Site(0,+1,0,1));
			neighbourSites.push_back(Site(0,0,+1,1));

			neighbourSites.push_back(Site(2,-1,0,0));
			neighbourSites.push_back(Site(2,0,-1,0));
			neighbourSites.push_back(Site(2,0,0,-1));

			neighbourSites.push_back(Site(-1,2,0,0));
			neighbourSites.push_back(Site(0,2,-1,0));
			neighbourSites.push_back(Site(0,2,0,-1));

			neighbourSites.push_back(Site(-1,0,2,0));
			neighbourSites.push_back(Site(0,-1,2,0));
			neighbourSites.push_back(Site(0,0,2,-1));

			neighbourSites.push_back(Site(-1,0,0,2));
			neighbourSites.push_back(Site(0,-1,0,2));
			neighbourSites.push_back(Site(0,0,-1,2));

			neighbourSites.push_back(Site(2,+1,0,0));
			neighbourSites.push_back(Site(2,0,+1,0));
			neighbourSites.push_back(Site(2,0,0,+1));

			neighbourSites.push_back(Site(+1,2,0,0));
			neighbourSites.push_back(Site(0,2,+1,0));
			neighbourSites.push_back(Site(0,2,0,+1));

			neighbourSites.push_back(Site(+1,0,2,0));
			neighbourSites.push_back(Site(0,+1,2,0));
			neighbourSites.push_back(Site(0,0,2,+1));

			neighbourSites.push_back(Site(+1,0,0,2));
			neighbourSites.push_back(Site(0,+1,0,2));
			neighbourSites.push_back(Site(0,0,+1,2));

			neighbourSites.push_back(Site(-2,+1,0,0));
			neighbourSites.push_back(Site(-2,0,+1,0));
			neighbourSites.push_back(Site(-2,0,0,+1));

			neighbourSites.push_back(Site(+1,-2,0,0));
			neighbourSites.push_back(Site(0,-2,+1,0));
			neighbourSites.push_back(Site(0,-2,0,+1));

			neighbourSites.push_back(Site(+1,0,-2,0));
			neighbourSites.push_back(Site(0,+1,-2,0));
			neighbourSites.push_back(Site(0,0,-2,+1));

			neighbourSites.push_back(Site(+1,0,0,-2));
			neighbourSites.push_back(Site(0,+1,0,-2));
			neighbourSites.push_back(Site(0,0,+1,-2));

			neighbourSites.push_back(Site(-2,-1,0,0));
			neighbourSites.push_back(Site(-2,0,-1,0));
			neighbourSites.push_back(Site(-2,0,0,-1));

			neighbourSites.push_back(Site(-1,-2,0,0));
			neighbourSites.push_back(Site(0,-2,-1,0));
			neighbourSites.push_back(Site(0,-2,0,-1));

			neighbourSites.push_back(Site(-1,0,-2,0));
			neighbourSites.push_back(Site(0,-1,-2,0));
			neighbourSites.push_back(Site(0,0,-2,-1));

			neighbourSites.push_back(Site(-1,0,0,-2));
			neighbourSites.push_back(Site(0,-1,0,-2));
			neighbourSites.push_back(Site(0,0,-1,-2));
		}
};

}

#endif