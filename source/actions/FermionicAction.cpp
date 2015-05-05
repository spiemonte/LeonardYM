/*
 * FermionicAction.cpp
 *
 *  Created on: Mar 28, 2012
 *      Author: spiem_01
 */

#include "FermionicAction.h"

namespace Update {

FermionicAction::FermionicAction(DiracOperator* _diracOperator) : diracOperator(_diracOperator) { }

FermionicAction::~FermionicAction() {
	delete diracOperator;
}

DiracOperator* FermionicAction::getDiracOperator() const {
	return diracOperator;
}

void FermionicAction::setDiracOperator(DiracOperator* _diracOperator) {
	diracOperator = _diracOperator;
}

} /* namespace Update */
