/*
 * LayoutMap.h
 *
 *  Created on: Feb 27, 2014
 *      Author: spiem_01
 */

#ifndef LAYOUTMAP_H_
#define LAYOUTMAP_H_

namespace Lattice {

template<typename TLayoutFrom, typename TLayoutTo> class LayoutMap {
public:
	LayoutMap() {
		map = layoutFrom.getMapIndex(layoutTo);
	}

	int at(int site) const {
		return map[site];
	}
private:
	int* map;
	TLayoutTo layoutTo;
	TLayoutTo layoutFrom;
};

}



#endif /* LAYOUTMAP_H_ */
