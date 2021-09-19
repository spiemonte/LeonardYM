#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "RandomSeed.h"

int Update::RandomSeed::counter = -1;
std::mt19937 Update::RandomSeed::rng;
std::uniform_int_distribution<> Update::RandomSeed::dist = std::uniform_int_distribution<>(-1000000000,1000000000);


namespace Update {

RandomSeed::RandomSeed() { }

RandomSeed::~RandomSeed() { }

int RandomSeed::randomSeed() {
	variate_generator<std::mt19937&, std::uniform_int_distribution<> > randomInteger(rng, dist);
	++counter;
#ifndef ENABLE_MPI
	return abs(seed[counter % 23]+time(NULL)+randomInteger());
#endif
#ifdef ENABLE_MPI
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	return abs(seed[counter % 128]+mpiseed[world_rank % 512]+((12*world_rank + 13) % (512 - 1))+time(NULL)+randomInteger());
#endif
}

variate_generator<std::mt19937&, std::normal_distribution<> > RandomSeed::getNormalNumberGenerator(std::mt19937& gen, double sd) {
	std::normal_distribution<> normal(0.,sd);
	return variate_generator<std::mt19937&, std::normal_distribution<> >(gen, normal);
}

variate_generator<std::mt19937&, std::uniform_real_distribution<> > RandomSeed::getRandomNumberGenerator(std::mt19937& gen) {
	std::uniform_real_distribution<> dist(0.0,1.0);
	return variate_generator<std::mt19937&, std::uniform_real_distribution<> >(gen, dist);
}

variate_generator<std::mt19937&, std::uniform_int_distribution<> > RandomSeed::getRandomIntegerGenerator(std::mt19937& gen) {
	std::uniform_int_distribution<> dist(0,1);
	return variate_generator<std::mt19937&, std::uniform_int_distribution<> >(gen, dist);
}

} /* namespace Update */
