#ifndef SCALARACTION_H
#define SCALARACTION_H
#include "Environment.h"
#include "hmc_forces/Force.h"
#include "Energy.h"

namespace Update {

class ScalarAction :  public Energy, public Force {
public:
        ScalarAction(real_t _mu, real_t _lambda, real_t _lambda_8);
        ~ScalarAction();

#ifndef __IBMCPP__
        using Force::force;
#endif
        void setLambda(const real_t&);
        real_t getLambda() const;

	void setAdjointLambda(const real_t&);
        real_t getAdjointLambda() const;

        void setMu(const real_t&);
        real_t getMu() const;
protected:
        real_t m;
        real_t lambda;
	real_t lambda_8;
};

} /* namespace Update */

#endif /* ADJOINTSCALARACTION_H_ */

