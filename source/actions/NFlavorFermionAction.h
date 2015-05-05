/*
 * NFlavorFermionAction.h
 *
 *  Created on: Aug 8, 2012
 *      Author: spiem_01
 */

#ifndef NFLAVORFERMIONACTION_H_
#define NFLAVORFERMIONACTION_H_

#include "dirac_functions/RationalApproximation.h"
#include "dirac_operators/DiracOperator.h"
#include "FermionicAction.h"

#include <vector>

namespace Update {

class NFlavorFermionAction : public FermionicAction {
public:
	NFlavorFermionAction(DiracOperator* _squareDiracOperator, DiracOperator* _diracOperator, const std::vector<RationalApproximation>& _rationalApproximations);
	virtual ~NFlavorFermionAction();

	virtual GaugeGroup force(const environment_t& env, int site, int mu) const;

	virtual void updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env);

	virtual long_real_t energy(const environment_t& env);

	void addPseudoFermion(extended_dirac_vector_t* _pseudofermion);
	void cleanPseudoFermions();
	std::vector<extended_dirac_vector_t*> getPseudoFermion() const;

	double getForcePrecision() const;
	void setForcePrecision(double precision);

	int getForceMaxIterations() const;
	void setForceMaxIterations(int iterations);
private:
	NFlavorFermionAction(const NFlavorFermionAction& ) : FermionicAction(NULL) { }

	DiracOperator* squareDiracOperator;
	//The fermion force
	FermionForce* fermionForce;
	//The precision used in the force inversion
	double forcePrecision;
	//The maximum number of iterations used by the cg solver
	int maxIterations;
	//The pseudofermions field
	std::vector<extended_dirac_vector_t*> pseudofermions;
	//The rational function approximation
	std::vector<RationalApproximation> rationalApproximations;
	//Static vector of the dirac_vector needed for the calculation of the force
	std::vector< std::vector<extended_dirac_vector_t> > Xs;
	std::vector< std::vector<extended_dirac_vector_t> > Ys;
	extended_dirac_vector_t* tmp_pseudofermion;

	MultishiftSolver* multishiftSolver;
};

} /* namespace Update */
#endif /* NFLAVORFERMIONACTION_H_ */
