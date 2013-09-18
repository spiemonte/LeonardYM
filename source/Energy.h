/*
 * Energy.h
 *
 *  Created on: Mar 26, 2012
 *      Author: spiem_01
 */

#ifndef ENERGY_H_
#define ENERGY_H_
#include "Environment.h"

namespace Update {

class Energy {
public:
	Energy();
	virtual ~Energy();

	virtual long_real_t energy(const environment_t& env) = 0;
};

} /* namespace Update */
#endif /* ENERGY_H_ */
