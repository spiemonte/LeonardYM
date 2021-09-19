#ifndef SCALARACTION_H
#define SCALARACTION_H
#include "Environment.h"
#include "hmc_forces/Force.h"
#include "Energy.h"

namespace Update {

class ScalarAction : public Energy, public Force {
public:
        ScalarAction();
        ScalarAction(const real_t& _mass_adjoint, const real_t& _mass_fundamental, const real_t& _lambda_fundamental, const real_t& _lambda_adjoint, const real_t& _lambda_mixed, const real_t& _lambda_8);
        ~ScalarAction();

#ifndef __IBMCPP__
        using Force::force;
#endif
	AdjointRealVector getKineticCoupling(const extended_adjoint_lattice_t& adjointLinkConfiguration, const extended_adjoint_real_color_vector_t& scalar_field, int site) const;
        FundamentalVector getKineticCoupling(const extended_gauge_lattice_t& linkConfiguration, const extended_color_vector_t& scalar_field, int site) const;
	real_t deltaEnergy(const AdjointRealVector& kineticCoupling, const AdjointRealVector& old_field, const AdjointRealVector& new_field) const;

        real_t energy(const std::vector<AdjointRealVector const*>& adjoint_fields, const std::vector<FundamentalVector const*>& fundamental_fields, const AdjointRealVector& kineticCoupling, unsigned int index_field_kinetic_coupling) const;
        real_t energy(const std::vector<AdjointRealVector const*>& adjoint_fields, const std::vector<FundamentalVector const*>& fundamental_fields, const FundamentalVector& kineticCoupling, unsigned int index_field_kinetic_coupling) const;
        real_t potentialEnergy(const std::vector<AdjointRealVector const*>& adjoint_fields, const std::vector<FundamentalVector const*>& fundamental_fields) const;

        real_t energy(const std::vector<AdjointRealVector*>& adjoint_fields, const std::vector<FundamentalVector*>& fundamental_fields, const AdjointRealVector& kineticCoupling, unsigned int index_field_kinetic_coupling) const;
        real_t energy(const std::vector<AdjointRealVector*>& adjoint_fields, const std::vector<FundamentalVector*>& fundamental_fields, const FundamentalVector& kineticCoupling, unsigned int index_field_kinetic_coupling) const;
        real_t potentialEnergy(const std::vector<AdjointRealVector*>& adjoint_fields, const std::vector<FundamentalVector*>& fundamental_fields) const;

        std::vector<AdjointRealVector*> scalarFieldsAtSite(std::vector<extended_adjoint_real_color_vector_t>& adjoint_scalar_fields, unsigned int site) const;
        std::vector<FundamentalVector*> scalarFieldsAtSite(std::vector<extended_color_vector_t>& fundamental_scalar_fields, unsigned int site) const;

        std::vector<AdjointRealVector const*> scalarFieldsAtSite(const std::vector<extended_adjoint_real_color_vector_t>& adjoint_scalar_fields, unsigned int site) const;
        std::vector<FundamentalVector const*> scalarFieldsAtSite(const std::vector<extended_color_vector_t>& fundamental_scalar_fields, unsigned int site) const;

	virtual GaugeGroup force(const environment_t& env, int site, int mu) const;
        virtual void updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env);
        virtual long_real_t energy(const environment_t& env);

private:
	LieGenerator<AdjointGroup> adjointLieGenerator;
        LieGenerator<GaugeGroup> gaugeLieGenerator;

        real_t mass_adjoint;
        real_t mass_fundamental;
        real_t lambda_fundamental;
        real_t lambda_adjoint;
        real_t lambda_mixed;
        real_t lambda_8;

        //virtual void updateForce(extended_adjoint_color_vector_t& forceLattice, const environment_t& env) const;
};

} /* namespace Update */

#endif /* ADJOINTSCALARACTION_H_ */

