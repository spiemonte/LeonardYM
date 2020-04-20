#ifndef RANDOMSCALARUPDATER_H_
#define RANDOMSCALARUPDATER_H_
#include "LatticeSweep.h"
#include "Environment.h"
#include "utils/RandomSeed.h"

namespace Update {

class RandomScalarUpdater : public Update::LatticeSweep {
public:
	RandomScalarUpdater();
	RandomScalarUpdater(const RandomScalarUpdater&);
	~RandomScalarUpdater();
	
	virtual void execute(environment_t& environment);
protected:
	void generateGaussianAdjointScalar(extended_adjoint_color_vector_t& vector);
	void generateGaussianFundamentalScalar(extended_color_vector_t& vector);

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
#endif /* SCALARUPDATER_H_ */
