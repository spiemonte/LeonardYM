#ifndef FORCE_H_
#define FORCE_H_
#include "Environment.h"

namespace Update {

class Force {
public:
	Force();
	virtual ~Force();

	virtual GaugeGroup force(const environment_t& env, int site, int mu) const = 0;

	virtual void updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env);
};

} /* namespace Update */
#endif /* FORCE_H_ */
