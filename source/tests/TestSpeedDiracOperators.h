#ifndef TESTSPEEDDIRACOPERATORS_H_
#define TESTSPEEDDIRACOPERATORS_H_

#include "LatticeSweep.h"

namespace Update {

class TestSpeedDiracOperators : public Update::LatticeSweep {
public:
	TestSpeedDiracOperators();
	~TestSpeedDiracOperators();

	virtual void execute(environment_t& environment);

	static void registerParameters(po::options_description& desc);
};

} /* namespace Update */
#endif /* TESTLINEARALGEBRA_H_ */
