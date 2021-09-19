#ifndef GLUEBALL_H_
#define GLUEBALL_H_
#include "LatticeSweep.h"

namespace Update {

class Glueball : public LatticeSweep {
public:
	Glueball();
	~Glueball();

	virtual void execute(environment_t& environment);

	static void registerParameters(std::map<std::string, Option>& desc);
};

} /* namespace Update */
#endif /* GLUEBALL_H_ */
