#ifndef FUNDAMENTALMETROPOLISSCALARUPDATER_H_
#define FUNDAMENTALMETROPOLISSCALARUPDATER_H_
#include "LatticeSweep.h"
#include "Environment.h"
#include "utils/RandomSeed.h"

namespace Update {

class FundamentalMetropolisScalarUpdater : public Update::LatticeSweep {
public:
	FundamentalMetropolisScalarUpdater();
	FundamentalMetropolisScalarUpdater(const FundamentalMetropolisScalarUpdater& toCopy);
	~FundamentalMetropolisScalarUpdater();
	
	virtual void execute(environment_t& environment);

	static void registerParameters(po::options_description& desc);

private:
#ifndef MULTITHREADING
	//The generator of random numbers
	random_generator_t randomGenerator;
	//The generator of normal random numbers
	random_normal_generator_t randomNormal;
	//The generator of random numbers, uniform distribution [0,1]
        random_uniform_generator_t randomUniform;
#endif
#ifdef MULTITHREADING
	//The generators of random numbers
	random_generator_t** randomGenerator;
	//The generators of normal random numbers
	random_normal_generator_t** randomNormal;
	//The generator of random numbers, uniform distribution [0,1]
        random_uniform_generator_t** randomUniform;
#endif
};

} /* namespace Update */
#endif /* METROPOLISSCALARUPDATER_H_ */
