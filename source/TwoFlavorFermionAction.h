/*
 * TwoFlavorFermionAction.h
 *
 *  Created on: Aug 8, 2012
 *      Author: spiem_01
 */

#ifndef TWOFLAVORFERMIONACTION_H_
#define TWOFLAVORFERMIONACTION_H_

#include "FermionForce.h"
#include "GaugeAction.h"
#include "dirac_operators/DiracOperator.h"
#include "Energy.h"
#include "FermionicAction.h"

namespace Update {

class TwoFlavorFermionAction : public FermionicAction {
public:
	TwoFlavorFermionAction(DiracOperator* _diracOperator);
	~TwoFlavorFermionAction();

	virtual GaugeGroup force(const environment_t& env, int site, int mu) const;

	virtual void updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env);

	virtual long_real_t energy(const environment_t& env);

	void setPseudoFermion(extended_dirac_vector_t* _pseudofermion);
	extended_dirac_vector_t* getPseudoFermion() const;

	double getForcePrecision() const;
	void setForcePrecision(double precision);
private:
	//The fermion force
	FermionForce* fermionForce;
	//The precision for the force
	double forcePrecision;
	//The pseudofermion field
	extended_dirac_vector_t* pseudofermion;
	//The vector needed for the calculation of the force
	extended_dirac_vector_t X;
	extended_dirac_vector_t Y;
};

} /* namespace Update */
#endif /* TWOFLAVORFERMIONACTION_H_ */
