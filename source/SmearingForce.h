/*
 * SmearingForce.h
 *
 *  Created on: Dec 7, 2012
 *      Author: spiem_01
 */

#ifndef SMEARINGFORCE_H_
#define SMEARINGFORCE_H_

namespace Update {

class SmearingForce {
public:
	SmearingForce();
	~SmearingForce();

	FermionicForceMatrix derivative(environment_t& env, const FermionicForceMatrix& unsmearedDerivative, int site, unsigned int mu);
private:

};

} /* namespace Update */
#endif /* SMEARINGFORCE_H_ */
