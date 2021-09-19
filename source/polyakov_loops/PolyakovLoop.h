#ifndef POLYAKOVLOOP_H_
#define POLYAKOVLOOP_H_

#include "LatticeSweep.h"

namespace Update {

class PolyakovLoop: public Update::LatticeSweep {
public:
	PolyakovLoop();
	~PolyakovLoop();

	virtual void execute(environment_t& environment);
protected:
	void write_polyakov_loop_config(const extended_gauge_lattice_t& polyakov, const std::string& output_name, const std::string& output_directory, int offset, const std::complex<real_t>& average_polyakov_loop) const;
};

} /* namespace Update */
#endif /* POLYAKOVLOOP_H_ */
