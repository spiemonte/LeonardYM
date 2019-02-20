#ifndef FUNDAMENTALSCALARACTION_H
#define FUNDAMENTALSCALARACTION_H
#include "Environment.h"
#include "ScalarAction.h"

namespace Update {

class FundamentalScalarAction :  public ScalarAction {
public:
        FundamentalScalarAction(real_t _mu, real_t _lambda, real_t _lambda_8);
        ~FundamentalScalarAction();

#ifndef __IBMCPP__
        using Force::force;
#endif
	FundamentalVector getKineticCoupling(const extended_gauge_lattice_t& fundamentalLinkConfiguration, const extended_color_vector_t& scalar_field, int site) const;
	real_t deltaEnergy(const FundamentalVector& kineticCoupling, const FundamentalVector& old_field, const FundamentalVector& new_field) const;

	virtual GaugeGroup force(const environment_t& env, int site, int mu) const;
        virtual void updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env);
        virtual long_real_t energy(const environment_t& env);
private:
	LieGenerator<GaugeGroup> gaugeLieGenerator;

};

} /* namespace Update */

#endif /* ADJOINTSCALARACTION_H_ */

