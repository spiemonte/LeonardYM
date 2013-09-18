/*
 * RandomSeed.h
 *
 *  Created on: Feb 28, 2012
 *      Author: spiem_01
 */

#ifndef RANDOMSEED_H_
#define RANDOMSEED_H_
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>


namespace Update {

const int seed[] = {122957928, 317424862, 248759476, 149829867, 578467498, 537848001, 37523957, 24972985, 18927985, 32955921, 234290626, 108757876, 56729867, 468469498, 487529230, 132616280, 67839920, 84193034, 48684633, 7345549, 3648226, 19768548, 89247769};
const int mpiseed[] = {89376783, 456739568, 96748667, 138626485, 835883486, 113649871, 355747, 1746174, 16361136, 9856722, 13437, 83568373, 27682347, 175847892, 123546748, 12553764, 139579578, 13882249, 9696739, 67153573, 7719687, 3496733, 76191687};

class RandomSeed {
public:
	RandomSeed();
	~RandomSeed();

	static int randomSeed();
	static boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > getNormalNumberGenerator(boost::mt19937& gen, double sd = 1.);
	static boost::variate_generator<boost::mt19937&, boost::uniform_01<> > getRandomNumberGenerator(boost::mt19937& gen);
	static boost::variate_generator<boost::mt19937&, boost::uniform_int<> > getRandomIntegerGenerator(boost::mt19937& gen);
private:
	static int counter;
	static boost::mt19937 rng;
	static boost::uniform_int<> dist;
};

typedef boost::mt19937 random_generator_t;
typedef boost::variate_generator<boost::mt19937&, boost::uniform_01<> > random_uniform_generator_t;
typedef boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > random_normal_generator_t;
typedef boost::variate_generator<boost::mt19937&, boost::uniform_int<> > random_integer_generator_t;

} /* namespace Update */
#endif /* RANDOMSEED_H_ */
