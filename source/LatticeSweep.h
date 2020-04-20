#ifndef LATTICESWEEP_H_
#define LATTICESWEEP_H_
#include "Environment.h"

namespace po = boost::program_options;

namespace Update {

class LatticeSweep {
public:
	LatticeSweep();
	LatticeSweep(unsigned int _numberTimes, unsigned int _sweepToJump);
	virtual ~LatticeSweep();

	//Factory method to create a given Sweep
	static LatticeSweep* getInstance(const std::string& name);

	//Read the sweep list from the configuration file
	static LatticeSweep* read(const std::string& toRead);

	//Call the sweep as many times as specified by the configuration script
	void call(environment_t& environment);

	//virtual method to execute the sweep on the environment
	virtual void execute(environment_t& environment) = 0;

	//Set the number of calls of the sweep
	void setNumberTimes(unsigned int _numberTimes);
	void setSweepToJump(unsigned int _sweepToJump);

	static void printSweepsName();

	//Register all parameters of all sweeps
	static void addParameters(po::options_description& desc);

	//Static function to be overridden by the sweep implementation
	static void registerParameters(po::options_description& desc);
private:
	unsigned int numberTimes;
	unsigned int sweepToJump;
};

} /* namespace Update */
#endif /* LATTICESWEEP_H_ */
