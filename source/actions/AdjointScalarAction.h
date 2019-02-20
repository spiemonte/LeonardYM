#ifndef ADJOINTSCALARACTION_H
#define ADJOINTSCALARACTION_H
#include "Environment.h"
#include "ScalarAction.h"

namespace Update {

class AdjointScalarAction :  public ScalarAction {
public:
        AdjointScalarAction(real_t _mu, real_t _lambda, real_t _lambda_8);
        ~AdjointScalarAction();

#ifndef __IBMCPP__
        using Force::force;
#endif
	AdjointVector getKineticCoupling(const extended_adjoint_lattice_t& adjointLinkConfiguration, const extended_adjoint_color_vector_t& scalar_field, int site) const;
	real_t deltaEnergy(const AdjointVector& kineticCoupling, const AdjointVector& old_field, const AdjointVector& new_field) const;

	virtual GaugeGroup force(const environment_t& env, int site, int mu) const;
        virtual void updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env);
        virtual long_real_t energy(const environment_t& env);

private:
	LieGenerator<AdjointGroup> adjointLieGenerator;

        //virtual void updateForce(extended_adjoint_color_vector_t& forceLattice, const environment_t& env) const;
};

} /* namespace Update */

#endif /* ADJOINTSCALARACTION_H_ */

