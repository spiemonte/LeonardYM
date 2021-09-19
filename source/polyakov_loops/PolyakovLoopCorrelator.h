#ifndef POLYAKOVLOOPCORRELATOR_H_
#define POLYAKOVLOOPCORRELATOR_H_

#include "LatticeSweep.h"

namespace Update {

class PolyakovLoopCorrelator : public Update::LatticeSweep {
public:
	PolyakovLoopCorrelator();
	~PolyakovLoopCorrelator();

	virtual void execute(environment_t& environment);

	static void registerParameters(std::map<std::string, Option>&);
};

} /* namespace Update */
#endif /* POLYAKOVLOOPCORRELATOR_H_ */
