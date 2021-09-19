#ifndef MAXIMALABELIANPROJECTION_H
#define MAXIMALABELIANPROJECTION_H

#include "LatticeSweep.h"

namespace Update {

class MaximalAbelianProjection : public Update::LatticeSweep {
public:
	MaximalAbelianProjection();
	MaximalAbelianProjection(const MaximalAbelianProjection&);
	~MaximalAbelianProjection();

	void execute(environment_t& environment);
};

}

#endif

