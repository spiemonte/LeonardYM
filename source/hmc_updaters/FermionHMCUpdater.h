#ifndef FERMIONHCMUPDATER_H_
#define FERMIONHCMUPDATER_H_

#include "HMCUpdater.h"

namespace Update {

class FermionHMCUpdater : public HMCUpdater {
public:
	FermionHMCUpdater();
	~FermionHMCUpdater();

	void generateGaussianDiracVector(extended_dirac_vector_t& vector);

private:
#ifndef MULTITHREADING
	//The generator of random numbers
	random_generator_t randomGenerator;
	//The generator of normal random numbers
	random_normal_generator_t randomNormal;
#endif
#ifdef MULTITHREADING
	//The generators of random numbers
	random_generator_t** randomGenerator;
	//The generators of normal random numbers
	random_normal_generator_t** randomNormal;
#endif
};

} /* namespace Update */
#endif /* FERMIONHCMUPDATER_H_ */
