#ifndef WILSONLOOP_H_
#define WILSONLOOP_H_
#include "LatticeSweep.h"

namespace Update {

class WilsonLoop : public LatticeSweep {
public:
	WilsonLoop();
	WilsonLoop(const WilsonLoop& toCopy);
	~WilsonLoop();

	virtual void execute(environment_t& environment);

	static void registerParameters(std::map<std::string, Option>& desc);
private:
	reduced_matrix_lattice_t** tWilsonLine;
	reduced_matrix_lattice_t** xWilsonLine;
	int RMax;
	int TMax;
};

} /* namespace Update */
#endif /* WILSONLOOP_H_ */
