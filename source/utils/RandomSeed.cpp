/*
 * RandomSeed.cpp
 *
 *  Created on: Feb 28, 2012
 *      Author: spiem_01
 */

#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "RandomSeed.h"


namespace Update {

RandomSeed::RandomSeed() { }

RandomSeed::~RandomSeed() { }

int RandomSeed::randomSeed() {
	boost::variate_generator<boost::mt19937&, boost::uniform_int<> > randomInteger(rng, dist);
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

boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > RandomSeed::getNormalNumberGenerator(boost::mt19937& gen, double sd) {
	boost::normal_distribution<> normal(0.,sd);
	return boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >(gen, normal);
}

boost::variate_generator<boost::mt19937&, boost::uniform_01<> > RandomSeed::getRandomNumberGenerator(boost::mt19937& gen) {
	boost::uniform_01<> dist;
	return boost::variate_generator<boost::mt19937&, boost::uniform_01<> >(gen, dist);
}

boost::variate_generator<boost::mt19937&, boost::uniform_int<> > RandomSeed::getRandomIntegerGenerator(boost::mt19937& gen) {
	boost::uniform_int<> dist(0,1);
	return boost::variate_generator<boost::mt19937&, boost::uniform_int<> >(gen, dist);
}

} /* namespace Update */
