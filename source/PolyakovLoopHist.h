/*
 * PolyakovLoop.h
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

#ifndef POLYAKOVLOOP_H_
#define POLYAKOVLOOP_H_

#include "LatticeSweep.h"
#include "histogrammer.h"

namespace Update {

class PolyakovLoop: public Update::LatticeSweep {
public:
	PolyakovLoop();
	~PolyakovLoop();

	virtual void execute(environment_t& environment);
protected:
	Histogrammer2D* hist2d_;
    Histogrammer* evhist_;
    Histogrammer2D* evhist2d_;
	// add here histogrammer
};

} /* namespace Update */
#endif /* POLYAKOVLOOP_H_ */
