#include "FundamentalScalarAction.h"

namespace Update {

FundamentalScalarAction::FundamentalScalarAction(real_t _mu, real_t _lambda, real_t _lambda_8) : ScalarAction(_mu, _lambda, _lambda_8), gaugeLieGenerator() { }

FundamentalScalarAction::~FundamentalScalarAction() { }

FundamentalVector FundamentalScalarAction::getKineticCoupling(const extended_gauge_lattice_t& fundamentalLinkConfiguration, const extended_color_vector_t& scalar_field, int site) const {
	FundamentalVector kinetic_coupling;
	set_to_zero(kinetic_coupling);
	typedef extended_color_vector_t LT;

	for (unsigned int mu = 0; mu < 4; ++mu) {
		//kinetic_coupling += - adjointLinkConfiguration[site][mu]*adjointLinkConfiguration[LT::sup(site,mu)][mu]*scalar_field[LT::sup(LT::sup(site,mu),mu)] - htrans(adjointLinkConfiguration[LT::sdn(site,mu)][mu])*htrans(adjointLinkConfiguration[LT::sdn(LT::sdn(site,mu),mu)][mu])*scalar_field[LT::sdn(LT::sdn(site,mu),mu)];
		kinetic_coupling += - fundamentalLinkConfiguration[site][mu]*scalar_field[LT::sup(site,mu)] - htrans(fundamentalLinkConfiguration[LT::sdn(site,mu)][mu])*scalar_field[LT::sdn(site,mu)];
	}

	return kinetic_coupling;
}

real_t FundamentalScalarAction::deltaEnergy(const FundamentalVector& kineticCoupling, const FundamentalVector& old_field, const FundamentalVector& new_field) const {
	real_t old_square = real(vector_dot(old_field,old_field));
	real_t old_potential = (-m-8.)*old_square/2. - (lambda/24.)*old_square*old_square;
	real_t old_adjoint_potential = 0.;
	//For every generator
	for (unsigned int i = 0; i < gaugeLieGenerator.numberGenerators(); ++i) {
		old_square = real(vector_dot(old_field,gaugeLieGenerator.get(i)*old_field));
		old_adjoint_potential -= (lambda_8/24.)*old_square*old_square;
	}
	real_t old_kinetic = -real(vector_dot(kineticCoupling,old_field));

	real_t new_square = real(vector_dot(new_field,new_field));
        real_t new_potential = (-m-8.)*new_square/2. - (lambda/24.)*new_square*new_square;
	real_t new_adjoint_potential = 0.;
	//For every generator
	for (unsigned int i = 0; i < gaugeLieGenerator.numberGenerators(); ++i) {
		new_square = real(vector_dot(new_field,gaugeLieGenerator.get(i)*new_field));
		new_adjoint_potential -= (lambda_8/24.)*new_square*new_square;
	}
	real_t new_kinetic = -real(vector_dot(kineticCoupling,new_field));

	return -(new_kinetic + new_potential + new_adjoint_potential - old_kinetic - old_potential - old_adjoint_potential);
}

long_real_t FundamentalScalarAction::energy(const environment_t& env) {
	long_real_t energy = 0.;

	typedef extended_color_vector_t LT;
	std::vector<extended_color_vector_t>::const_iterator scalar_field;

	for (scalar_field = env.fundamental_scalar_fields.begin(); scalar_field < env.fundamental_scalar_fields.end(); ++scalar_field) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
#pragma omp parallel for reduction(+:energy)
			for (int site = 0; site < scalar_field->localsize; ++site) {
				//FundamentalVector dmu = (env.getFundamentalLattice()[site][mu]*(*scalar_field)[LT::sup(site,mu)] - htrans(env.getFundamentalLattice()[LT::sdn(site,mu)][mu])*(*scalar_field)[LT::sdn(site,mu)]);
				FundamentalVector dmu = (env.getFundamentalLattice()[site][mu]*(*scalar_field)[LT::sup(site,mu)] - (*scalar_field)[site]);

				for (int c = 0; c < numberColors; ++c) {
					energy += 0.5*real(conj(dmu[c])*(dmu[c]));
				}
			}

		}

#pragma omp parallel for reduction(+:energy)
		for (int site = 0; site < scalar_field->localsize; ++site) {
			real_t square = 0;
			for (int c = 0; c < numberColors; ++c) {
				real_t sq = real(conj((*scalar_field)[site][c])*((*scalar_field)[site][c]));
				square += sq;
				energy += m*sq/2.;
			}
			energy += (lambda/24.)*square*square;

			//For every generator
			for (unsigned int i = 0; i < gaugeLieGenerator.numberGenerators(); ++i) {
				real_t sq = real(vector_dot((*scalar_field)[site],gaugeLieGenerator.get(i)*(*scalar_field)[site]));
				energy += (lambda_8/24.)*sq*sq;
			}
		}
	}

	reduceAllSum(energy);

	return energy;
}

void FundamentalScalarAction::updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env) {
	typedef extended_color_vector_t LT;

	std::vector<extended_color_vector_t>::const_iterator scalar_field;

#pragma omp parallel for
	for (int site = 0; site < forceLattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			set_to_zero(forceLattice[site][mu]);
		}
	}

	for (scalar_field = env.fundamental_scalar_fields.begin(); scalar_field < env.fundamental_scalar_fields.end(); ++scalar_field) {
		for (unsigned int mu = 0; mu < 4; ++mu) {

#pragma omp parallel for
			for (int site = 0; site < scalar_field->localsize; ++site) {
					//For every generator
					for (unsigned int i = 0; i < gaugeLieGenerator.numberGenerators(); ++i) {
		forceLattice[site][mu] -= 0.5*vector_dot((*scalar_field)[LT::sup(site,mu)], htrans(env.getFundamentalLattice()[site][mu])*gaugeLieGenerator.get(i)*(*scalar_field)[site])*gaugeLieGenerator.get(i);
		forceLattice[site][mu] += 0.5*vector_dot((*scalar_field)[site], gaugeLieGenerator.get(i)*env.getFundamentalLattice()[site][mu]*(*scalar_field)[LT::sup(site,mu)])*gaugeLieGenerator.get(i);
		//forceLattice[site][mu] -= 0.5*vector_dot((*scalar_field)[LT::sup(site,mu)], htrans(env.getFundamentalLattice()[site][mu])*fermionLieGenerator.get(i)*htrans(env.getFundamentalLattice()[LT::sdn(site,mu)][mu])*(*scalar_field)[LT::sdn(site,mu)])*gaugeLieGenerator.get(i);
		//forceLattice[site][mu] -= 0.5*vector_dot((*scalar_field)[LT::sup(LT::sup(site,mu),mu)], htrans(env.getFundamentalLattice()[LT::sup(site,mu)][mu])*htrans(env.getFundamentalLattice()[site][mu])*fermionLieGenerator.get(i)*(*scalar_field)[site])*gaugeLieGenerator.get(i);
		//forceLattice[site][mu] += 0.5*vector_dot((*scalar_field)[LT::sdn(site,mu)],env.getFundamentalLattice()[LT::sdn(site,mu)][mu]*fermionLieGenerator.get(i)*env.getFundamentalLattice()[site][mu]*(*scalar_field)[LT::sup(site,mu)])*gaugeLieGenerator.get(i);
		//forceLattice[site][mu] += 0.5*vector_dot((*scalar_field)[site], fermionLieGenerator.get(i)*env.getFundamentalLattice()[site][mu]*env.getFundamentalLattice()[LT::sup(site,mu)][mu]*(*scalar_field)[LT::sup(LT::sup(site,mu),mu)])*gaugeLieGenerator.get(i);	
					}
			}
		}
	}

	forceLattice.updateHalo();
}

GaugeGroup FundamentalScalarAction::force(const environment_t& env, int site, int mu) const {
	typedef extended_color_vector_t LT;

	GaugeGroup result;
	set_to_zero(result);

	std::vector<extended_color_vector_t>::const_iterator scalar_field;
	for (scalar_field = env.fundamental_scalar_fields.begin(); scalar_field < env.fundamental_scalar_fields.end(); ++scalar_field) {
		//For every generator
		for (unsigned int i = 0; i < gaugeLieGenerator.numberGenerators(); ++i) {
			result -= 0.5*vector_dot((*scalar_field)[LT::sup(site,mu)], htrans(env.getFundamentalLattice()[site][mu])*gaugeLieGenerator.get(i)*(*scalar_field)[site])*gaugeLieGenerator.get(i);
                	result += 0.5*vector_dot((*scalar_field)[site], gaugeLieGenerator.get(i)*env.getFundamentalLattice()[site][mu]*(*scalar_field)[LT::sup(site,mu)])*gaugeLieGenerator.get(i);
			/*result -= 0.5*vector_dot((*scalar_field)[LT::sup(site,mu)], htrans(env.getFundamentalLattice()[site][mu])*fermionLieGenerator.get(i)*htrans(env.getFundamentalLattice()[LT::sdn(site,mu)][mu])*(*scalar_field)[LT::sdn(site,mu)])*gaugeLieGenerator.get(i);
                	result -= 0.5*vector_dot((*scalar_field)[LT::sup(LT::sup(site,mu),mu)], htrans(env.getFundamentalLattice()[LT::sup(site,mu)][mu])*htrans(env.getFundamentalLattice()[site][mu])*fermionLieGenerator.get(i)*(*scalar_field)[site])*gaugeLieGenerator.get(i);
                	result += 0.5*vector_dot((*scalar_field)[LT::sdn(site,mu)],env.getFundamentalLattice()[LT::sdn(site,mu)][mu]*fermionLieGenerator.get(i)*env.getFundamentalLattice()[site][mu]*(*scalar_field)[LT::sup(site,mu)])*gaugeLieGenerator.get(i);
                	result += 0.5*vector_dot((*scalar_field)[site], fermionLieGenerator.get(i)*env.getFundamentalLattice()[site][mu]*env.getFundamentalLattice()[LT::sup(site,mu)][mu]*(*scalar_field)[LT::sup(LT::sup(site,mu),mu)])*gaugeLieGenerator.get(i);*/
                }
	}

	return result;
}

}

