#ifndef POLYAKOVLOOPEIGENVALUES_H_
#define POLYAKOVLOOPEIGENVALUES_H_

#include "LatticeSweep.h"

namespace Update {

class PolyakovLoopEigenvalues : public Update::LatticeSweep {
public:
	PolyakovLoopEigenvalues();
	~PolyakovLoopEigenvalues();

	virtual void execute(environment_t& environment);
};

} /* namespace Update */
#endif /* POLYAKOVLOOP_H_ */
