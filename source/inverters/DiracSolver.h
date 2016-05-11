#ifndef DIRACSOLVER_H
#define DIRACSOLVER_H

namespace Update {

class DiracSolver {
public:
	DiracSolver();
	virtual ~DiracSolver();

	static DiracSolver* getInstance(const std::string& solverName, const std::string& basename, const StorageParameters& sp);
	void registerParameters(po::options_description& desc, const std::string& basename);

	virtual bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, reduced_dirac_vector_t const* initial_guess = 0) = 0;

#ifdef ENABLE_MPI
	virtual bool solve(DiracOperator* dirac, const extended_dirac_vector_t& source, extended_dirac_vector_t& solution, extended_dirac_vector_t const* initial_guess = 0) = 0;
#endif

	void setPrecision(double _epsilon);
	double getPrecision() const;

	double getLastError() const;
	unsigned int getLastSteps() const;

	void setMaximumSteps(unsigned int _maxSteps);
	unsigned int getMaximumSteps() const;

protected:
	unsigned int maxSteps;
	double epsilon;
	double lastError;
	unsigned int lastSteps;
	
};

}

#endif
