#ifndef STANDARDSTENCIL_H
#define STANDARDSTENCIL_H
#include "Site.h"
#include <vector>

namespace Lattice {

class StandardStencil {
	public:
		//The neighbour sites that should be used
		static std::vector<Site> neighbourSites;

		static const int id;
		
		static void initializeNeighbourSites() {
			neighbourSites.clear();
			Site dx(1,0,0,0);
			Site dy(0,1,0,0);
			Site dz(0,0,1,0);
			Site dt(0,0,0,1);
			neighbourSites.push_back(dx);
			neighbourSites.push_back(dy);
			neighbourSites.push_back(dz);
			neighbourSites.push_back(dt);
			neighbourSites.push_back(-dx);
			neighbourSites.push_back(-dy);
			neighbourSites.push_back(-dz);
			neighbourSites.push_back(-dt);
			neighbourSites.push_back(dx-dy);
			neighbourSites.push_back(dx-dz);
			neighbourSites.push_back(dx-dt);
			neighbourSites.push_back(dy-dx);
			neighbourSites.push_back(dy-dz);
			neighbourSites.push_back(dy-dt);
			neighbourSites.push_back(dz-dx);
			neighbourSites.push_back(dz-dy);
			neighbourSites.push_back(dz-dt);
			neighbourSites.push_back(dt-dx);
			neighbourSites.push_back(dt-dy);
			neighbourSites.push_back(dt-dz);
		}
};

}

#endif