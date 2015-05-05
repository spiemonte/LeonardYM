/*
 * FermionicAction.h
 *
 *  Created on: Mar 28, 2012
 *      Author: spiem_01
 */

#ifndef FERMIONICACTION_H_
#define FERMIONICACTION_H_

#include "Energy.h"
#include "Force.h"
#include "dirac_operators/DiracOperator.h"

namespace Update {

class FermionicAction : public Update::Energy, public Update::Force {
public:
	FermionicAction(DiracOperator* _diracOperator);
	~FermionicAction();

	DiracOperator* getDiracOperator() const;
	void setDiracOperator(DiracOperator* _diracOperator);
protected:
	DiracOperator* diracOperator;
};

} /* namespace Update */
#endif /* FERMIONICACTION_H_ */
