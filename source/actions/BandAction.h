/*
 * BandAction.h
 *
 *  Created on: Jan 11, 2013
 *      Author: spiem_01
 */

#ifndef BANDACTION_H_
#define BANDACTION_H_
#include "Force.h"

#include <vector>
#include <utility>

namespace Update {

class BandAction : public Update::Force {
public:
	BandAction(Force* _subForce);
	~BandAction();

	virtual GaugeGroup force(const environment_t& env, int site, int mu) const;
	virtual void updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env);

	void addBand(int t0, int t1);
private:
	Force* subForce;
	std::vector< std::pair<int, int> > bands;
};

} /* namespace Update */
#endif /* BANDACTION_H_ */
