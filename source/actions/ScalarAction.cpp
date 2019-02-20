#include "ScalarAction.h"

namespace Update {

ScalarAction::ScalarAction(real_t _mu, real_t _lambda, real_t _lambda_8) : m(_mu), lambda(_lambda), lambda_8(_lambda_8) { }

ScalarAction::~ScalarAction() { }

void ScalarAction::setLambda(const real_t& _lambda) {
	lambda = _lambda;
}

real_t ScalarAction::getLambda() const {
	return lambda;
}

void ScalarAction::setAdjointLambda(const real_t& _lambda_8) {
	lambda_8 = _lambda_8;
}

real_t ScalarAction::getAdjointLambda() const {
	return lambda_8;
}

void ScalarAction::setMu(const real_t& _mu) {
	m = _mu;
}

real_t ScalarAction::getMu() const {
	return m;
}

}

