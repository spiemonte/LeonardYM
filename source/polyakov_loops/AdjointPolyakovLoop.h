#ifndef ADJOINTPOLYAKOVLOOP_H_
#define ADJOINTPOLYAKOVLOOP_H_
#include "LatticeSweep.h"

namespace Update {

class AdjointPolyakovLoop : public LatticeSweep {
public:
	AdjointPolyakovLoop();
	~AdjointPolyakovLoop();

	virtual void execute(environment_t& environment);
};

} /* namespace Update */
#endif /* ADJOINTPOLYAKOVLOOP_H_ */
