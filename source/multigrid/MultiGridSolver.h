#ifndef MULTIGRIDSOLVER_H
#define MULTIGRIDSOLVER_H
#include "BlockBasis.h"
#include "MultiGridBiConjugateGradient.h"
#include "dirac_operators/SAPPreconditioner.h"
#include "inverters/Solver.h"
#include <vector>

namespace Update {

class MultiGridSolver : public Solver {
	public:
		using Solver::solve;

		MultiGridSolver(int basisDimension, const std::vector<unsigned int>& _blockSize, BlockDiracOperator* _blackBlockDiracOperator, BlockDiracOperator* _redBlockDiracOperator);
		
		bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, reduced_dirac_vector_t const* initial_guess = 0);
		
		void initializeBasis(DiracOperator* dirac);

		void setSAPIterations(int _SAPIterantions);
		int getSAPIterations() const;

		void setSAPMaxSteps(int _SAPIterantions);
		int getSAPMaxSteps() const;

		void setSAPPrecision(real_t _SAPPrecision);
		real_t getSAPPrecision() const;

		void setGMRESIterations(int _GMRESIterations);
		int getGMRESIterations() const;

		void setGMRESPrecision(real_t _GMRESPrecision);
		real_t getGMRESPrecision() const;
		
		void setBasisDimension(unsigned int dim);
		unsigned int getBasisDimension() const;
	
		BlockBasis* getBasis();

	private:
		BlockBasis blockBasis;
		std::vector<unsigned int> blockSize;

		BlockDiracOperator* blackBlockDiracOperator;
		BlockDiracOperator* redBlockDiracOperator;
		ComplementBlockDiracOperator* K;

		MultiGridBiConjugateGradientSolver* biMgSolver;

		int SAPIterantions;
		int SAPMaxSteps;
		real_t SAPPrecision;

		int GMRESIterations;
		real_t GMRESPrecision;

		int BiMGIterations;
		real_t BiMGPrecision;

		reduced_dirac_vector_t r;
		reduced_dirac_vector_t c[15];
		reduced_dirac_vector_t u[15];
		reduced_dirac_vector_t mg_inverse;
		reduced_dirac_vector_t source_sap;
};

}

#endif

