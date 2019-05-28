#ifndef GAUGEFIXING_H
#define GAUGEFIXING_H

#include "LatticeSweep.h"
#include "utils/RandomSeed.h"

namespace Update {

class GaugeFixing {
public:
	GaugeFixing();
	GaugeFixing(const GaugeFixing&);
	~GaugeFixing();

	//Create a random gauge transformation
	void generateRandomGaugeTransformation(extended_matrix_lattice_t& gauge_transformation, real_t epsilon);

	//Apply a given gauge transformation gauge_transformation to the lattice
	void transform(extended_gauge_lattice_t& lattice, const extended_matrix_lattice_t& gauge_transformation) const;

	//Compute the A_\mu field from a given gauge fixed configuration
	void getLieAlgebraField(extended_gauge_lattice_t& Afield, const extended_gauge_lattice_t& lattice) const;

	//Metropolis for parellel tempering
	bool metropolis(real_t delta);

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
	//The generator of random numbers
	random_generator_t randomGeneratorUniform;
	//The generator of random numbers, uniform distribution [0,1]
	random_uniform_generator_t randomUniform;
};

}

#endif

