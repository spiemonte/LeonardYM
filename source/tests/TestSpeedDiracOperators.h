#ifndef TESTSPEEDDIRACOPERATORS_H_
#define TESTSPEEDDIRACOPERATORS_H_

#include "LatticeSweep.h"

namespace Update {

class TestSpeedDiracOperators : public Update::LatticeSweep {
public:
	TestSpeedDiracOperators();
	~TestSpeedDiracOperators();

	virtual void execute(environment_t& environment);

	static void registerParameters(std::map<std::string, Option>& desc);
};

} /* namespace Update */
#endif /* TESTLINEARALGEBRA_H_ */
