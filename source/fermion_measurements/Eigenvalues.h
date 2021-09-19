#ifndef EIGENVALUES_H_
#define EIGENVALUES_H_

#include "LatticeSweep.h"
#include "DiracEigenSolver.h"

namespace Update {

class Eigenvalues : public Update::LatticeSweep {
public:
	Eigenvalues();
	Eigenvalues(const Eigenvalues& toCopy);
	~Eigenvalues();

	virtual void execute(environment_t& environment);

	static void registerParameters(std::map<std::string, Option>& desc);
	
private:
	DiracEigenSolver* diracEigenSolver;
};

} /* namespace Update */
#endif /* EIGENVALUES_H_ */
