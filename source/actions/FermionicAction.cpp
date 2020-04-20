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
