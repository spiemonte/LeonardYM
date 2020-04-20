#ifndef TESTLINEARALGEBRA_H_
#define TESTLINEARALGEBRA_H_

#include "LatticeSweep.h"

namespace Update {

class TestLinearAlgebra: public Update::LatticeSweep {
public:
	TestLinearAlgebra();
	~TestLinearAlgebra();

	virtual void execute(environment_t& environment);
};

} /* namespace Update */
#endif /* TESTLINEARALGEBRA_H_ */
