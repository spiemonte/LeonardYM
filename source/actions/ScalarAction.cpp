#include "ScalarAction.h"

namespace Update {

ScalarAction::ScalarAction() : mass_adjoint(0), mass_fundamental(0), lambda_fundamental(0), lambda_adjoint(0), lambda_mixed(0), lambda_8(0) { }

ScalarAction::ScalarAction(const real_t& _mass_adjoint, const real_t& _mass_fundamental, const real_t& _lambda_fundamental, const real_t& _lambda_adjoint, const real_t& _lambda_mixed, const real_t& _lambda_8) : mass_adjoint(_mass_adjoint), mass_fundamental(_mass_fundamental), lambda_fundamental(_lambda_fundamental), lambda_adjoint(_lambda_adjoint), lambda_mixed(_lambda_mixed), lambda_8(_lambda_8) { }

ScalarAction::~ScalarAction() { }

AdjointRealVector ScalarAction::getKineticCoupling(const extended_adjoint_lattice_t& adjointLinkConfiguration, const extended_adjoint_real_color_vector_t& scalar_field, int site) const {
	AdjointRealVector kinetic_coupling;
	set_to_zero(kinetic_coupling);
	typedef extended_adjoint_real_color_vector_t LT;

	for (unsigned int mu = 0; mu < 4; ++mu) {
		//kinetic_coupling += - adjointLinkConfiguration[site][mu]*adjointLinkConfiguration[LT::sup(site,mu)][mu]*scalar_field[LT::sup(LT::sup(site,mu),mu)] - htrans(adjointLinkConfiguration[LT::sdn(site,mu)][mu])*htrans(adjointLinkConfiguration[LT::sdn(LT::sdn(site,mu),mu)][mu])*scalar_field[LT::sdn(LT::sdn(site,mu),mu)];
		kinetic_coupling += - adjointLinkConfiguration[site][mu]*scalar_field[LT::sup(site,mu)] - htrans(adjointLinkConfiguration[LT::sdn(site,mu)][mu])*scalar_field[LT::sdn(site,mu)];
	}

	return kinetic_coupling;
}

FundamentalVector ScalarAction::getKineticCoupling(const extended_gauge_lattice_t& linkConfiguration, const extended_color_vector_t& scalar_field, int site) const {
	FundamentalVector kinetic_coupling;
	set_to_zero(kinetic_coupling);
	typedef extended_color_vector_t LT;

	for (unsigned int mu = 0; mu < 4; ++mu) {
		//kinetic_coupling += - adjointLinkConfiguration[site][mu]*adjointLinkConfiguration[LT::sup(site,mu)][mu]*scalar_field[LT::sup(LT::sup(site,mu),mu)] - htrans(adjointLinkConfiguration[LT::sdn(site,mu)][mu])*htrans(adjointLinkConfiguration[LT::sdn(LT::sdn(site,mu),mu)][mu])*scalar_field[LT::sdn(LT::sdn(site,mu),mu)];
		kinetic_coupling += - linkConfiguration[site][mu]*scalar_field[LT::sup(site,mu)] - htrans(linkConfiguration[LT::sdn(site,mu)][mu])*scalar_field[LT::sdn(site,mu)];
	}

	return kinetic_coupling;
}

real_t ScalarAction::energy(const std::vector<AdjointRealVector const*>& adjoint_fields, const std::vector<FundamentalVector const*>& fundamental_fields, const AdjointRealVector& kineticCoupling, unsigned int index_field_kinetic_coupling) const {
	real_t kinetic = -vector_dot(kineticCoupling,*adjoint_fields[index_field_kinetic_coupling]);

	real_t potential = this->potentialEnergy(adjoint_fields, fundamental_fields);

	return kinetic + potential;
}

real_t ScalarAction::energy(const std::vector<AdjointRealVector const*>& adjoint_fields, const std::vector<FundamentalVector const*>& fundamental_fields, const FundamentalVector& kineticCoupling, unsigned int index_field_kinetic_coupling) const {
	real_t kinetic = -real(vector_dot(kineticCoupling,*fundamental_fields[index_field_kinetic_coupling]));

	real_t potential = this->potentialEnergy(adjoint_fields, fundamental_fields);

	return kinetic + potential;
}

real_t ScalarAction::potentialEnergy(const std::vector<AdjointRealVector const*>& adjoint_fields, const std::vector<FundamentalVector const*>& fundamental_fields) const {
	real_t square_fundamental = 0.;
	for (unsigned int i = 0; i < fundamental_fields.size(); ++i) {
		square_fundamental += real(vector_dot(*fundamental_fields[i],*fundamental_fields[i]));
	}

	real_t square_adjoint = 0.;
	for (unsigned int i = 0; i < adjoint_fields.size(); ++i) {
		square_adjoint += vector_dot(*adjoint_fields[i],*adjoint_fields[i]);
	}

	real_t quartic_traced_adjoint = 0.;
	//For all Lie algebra generators
	for (unsigned int j = 0; j < adjointLieGenerator.numberGenerators(); ++j) {
		real_t square_traced_adjoint = 0.;
		for (unsigned int i = 0; i < adjoint_fields.size(); ++i) {
			square_traced_adjoint += real(vector_dot(*adjoint_fields[i],adjointLieGenerator.get(j)*(*adjoint_fields[i])));
		}
		quartic_traced_adjoint += (lambda_8/24.)*square_traced_adjoint*square_traced_adjoint;
	}

	real_t mass_term = (mass_adjoint-8.)*square_adjoint/2. + (mass_fundamental-8.)*square_fundamental/2;

	return mass_term + quartic_traced_adjoint + (lambda_fundamental/24.)*square_fundamental*square_fundamental  + (lambda_adjoint/24.)*square_adjoint*square_adjoint + (lambda_mixed/24.)*square_adjoint*square_fundamental;
}

real_t ScalarAction::energy(const std::vector<AdjointRealVector*>& adjoint_fields, const std::vector<FundamentalVector*>& fundamental_fields, const AdjointRealVector& kineticCoupling, unsigned int index_field_kinetic_coupling) const {
	real_t kinetic = -vector_dot(kineticCoupling,*adjoint_fields[index_field_kinetic_coupling]);

	real_t potential = this->potentialEnergy(adjoint_fields, fundamental_fields);

	return kinetic + potential;
}

real_t ScalarAction::energy(const std::vector<AdjointRealVector*>& adjoint_fields, const std::vector<FundamentalVector*>& fundamental_fields, const FundamentalVector& kineticCoupling, unsigned int index_field_kinetic_coupling) const {
	real_t kinetic = -real(vector_dot(kineticCoupling,*fundamental_fields[index_field_kinetic_coupling]));

	real_t potential = this->potentialEnergy(adjoint_fields, fundamental_fields);

	return kinetic + potential;
}

real_t ScalarAction::potentialEnergy(const std::vector<AdjointRealVector*>& adjoint_fields, const std::vector<FundamentalVector*>& fundamental_fields) const {
	real_t square_fundamental = 0.;
	for (unsigned int i = 0; i < fundamental_fields.size(); ++i) {
		square_fundamental += real(vector_dot(*fundamental_fields[i],*fundamental_fields[i]));
	}

	real_t square_adjoint = 0.;
	for (unsigned int i = 0; i < adjoint_fields.size(); ++i) {
		square_adjoint += vector_dot(*adjoint_fields[i],*adjoint_fields[i]);
	}

	real_t quartic_traced_adjoint = 0.;
	//For all Lie algebra generators
	for (unsigned int j = 0; j < adjointLieGenerator.numberGenerators(); ++j) {
		real_t square_traced_adjoint = 0.;
		for (unsigned int i = 0; i < adjoint_fields.size(); ++i) {
			square_traced_adjoint += real(vector_dot(*adjoint_fields[i],adjointLieGenerator.get(j)*(*adjoint_fields[i])));
		}
		quartic_traced_adjoint += (lambda_8/24.)*square_traced_adjoint*square_traced_adjoint;
	}

	real_t mass_term = (mass_adjoint-8.)*square_adjoint/2. + (mass_fundamental-8.)*square_fundamental/2;

	return mass_term + quartic_traced_adjoint + (lambda_fundamental/24.)*square_fundamental*square_fundamental  + (lambda_adjoint/24.)*square_adjoint*square_adjoint + (lambda_mixed/24.)*square_adjoint*square_fundamental;
}

long_real_t ScalarAction::energy(const environment_t& env) {
	long_real_t energy = 0.;

	typedef extended_adjoint_real_color_vector_t LT;
	std::vector<extended_adjoint_real_color_vector_t>::const_iterator adjoint_scalar_field;

	for (adjoint_scalar_field = env.adjoint_scalar_fields.begin(); adjoint_scalar_field < env.adjoint_scalar_fields.end(); ++adjoint_scalar_field) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
#pragma omp parallel for reduction(+:energy)
			for (int site = 0; site < adjoint_scalar_field->localsize; ++site) {
				AdjointRealVector dmu = env.getAdjointLattice()[site][mu]*(*adjoint_scalar_field)[LT::sup(site,mu)];

				for (int c = 0; c < numberColors*numberColors - 1; ++c) {
					//Only nearest neighbors couplings, mass and quadratic terms set in the potential
					energy += ((*adjoint_scalar_field)[site][c]*dmu[c]);
				}
			}

		}
	}

	std::vector<extended_color_vector_t>::const_iterator fundamental_scalar_field;

	for (fundamental_scalar_field = env.fundamental_scalar_fields.begin(); fundamental_scalar_field < env.fundamental_scalar_fields.end(); ++fundamental_scalar_field) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
#pragma omp parallel for reduction(+:energy)
			for (int site = 0; site < fundamental_scalar_field->localsize; ++site) {
				FundamentalVector dmu = env.getFundamentalLattice()[site][mu]*(*fundamental_scalar_field)[LT::sup(site,mu)];

				for (int c = 0; c < numberColors; ++c) {
					//Mass and quadratic terms set in the potential
					energy += real(conj((*fundamental_scalar_field)[site][c])*dmu[c]);
				}
			}
		}
	}

#pragma omp parallel for reduction(+:energy)
	for (int site = 0; site < env.getAdjointLattice().localsize; ++site) {
		energy += this->potentialEnergy(scalarFieldsAtSite(env.adjoint_scalar_fields, site), scalarFieldsAtSite(env.fundamental_scalar_fields, site));
	}
	
	reduceAllSum(energy);

	return energy;
}

std::vector<AdjointRealVector*> ScalarAction::scalarFieldsAtSite(std::vector<extended_adjoint_real_color_vector_t>& adjoint_scalar_fields, unsigned int site) const {
	std::vector<AdjointRealVector*> adjoint_fields;

	for (unsigned int n = 0; n < adjoint_scalar_fields.size(); ++n) {
		adjoint_fields.push_back(&(adjoint_scalar_fields[n][site]));
	}

	return adjoint_fields;
}

std::vector<FundamentalVector*> ScalarAction::scalarFieldsAtSite(std::vector<extended_color_vector_t>& fundamental_scalar_fields, unsigned int site) const {
	std::vector<FundamentalVector*> fundamental_fields;

	for (unsigned int n = 0; n < fundamental_scalar_fields.size(); ++n) {
		fundamental_fields.push_back(&(fundamental_scalar_fields[n][site]));
	}

	return fundamental_fields;
}

std::vector<AdjointRealVector const*> ScalarAction::scalarFieldsAtSite(const std::vector<extended_adjoint_real_color_vector_t>& adjoint_scalar_fields, unsigned int site) const {
	std::vector<AdjointRealVector const*> adjoint_fields;

	for (unsigned int n = 0; n < adjoint_scalar_fields.size(); ++n) {
		adjoint_fields.push_back(&(adjoint_scalar_fields[n][site]));
	}

	return adjoint_fields;
}

std::vector<FundamentalVector const*> ScalarAction::scalarFieldsAtSite(const std::vector<extended_color_vector_t>& fundamental_scalar_fields, unsigned int site) const {
	std::vector<FundamentalVector const*> fundamental_fields;

	for (unsigned int n = 0; n < fundamental_scalar_fields.size(); ++n) {
		fundamental_fields.push_back(&(fundamental_scalar_fields[n][site]));
	}

	return fundamental_fields;
}

void ScalarAction::updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env) {
	typedef extended_adjoint_real_color_vector_t LT;

	std::vector<extended_adjoint_real_color_vector_t>::const_iterator adjoint_scalar_field;

#pragma omp parallel for
	for (int site = 0; site < forceLattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			set_to_zero(forceLattice[site][mu]);
		}
	}

	for (adjoint_scalar_field = env.adjoint_scalar_fields.begin(); adjoint_scalar_field < env.adjoint_scalar_fields.end(); ++adjoint_scalar_field) {
#pragma omp parallel for
		for (int site = 0; site < adjoint_scalar_field->localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
					//For every generator
					for (unsigned int i = 0; i < adjointLieGenerator.numberGenerators(); ++i) {
		forceLattice[site][mu] += std::complex<real_t>(0., imag(vector_dot((*adjoint_scalar_field)[LT::sup(site,mu)], htrans(env.getAdjointLattice()[site][mu])*adjointLieGenerator.get(i)*(*adjoint_scalar_field)[site])))*gaugeLieGenerator.get(i);
		//forceLattice[site][mu] -= (0.5*vector_dot((*scalar_field)[LT::sup(site,mu)], htrans(env.getAdjointLattice()[site][mu])*adjointLieGenerator.get(i)*(*scalar_field)[site]))*gaugeLieGenerator.get(i);
		//forceLattice[site][mu] += (0.5*vector_dot((*scalar_field)[site], adjointLieGenerator.get(i)*env.getAdjointLattice()[site][mu]*(*scalar_field)[LT::sup(site,mu)]))*gaugeLieGenerator.get(i);
		//forceLattice[site][mu] -= 0.5*vector_dot((*scalar_field)[LT::sup(site,mu)], htrans(env.getAdjointLattice()[site][mu])*adjointLieGenerator.get(i)*htrans(env.getAdjointLattice()[LT::sdn(site,mu)][mu])*(*scalar_field)[LT::sdn(site,mu)])*gaugeLieGenerator.get(i);
		//forceLattice[site][mu] -= 0.5*vector_dot((*scalar_field)[LT::sup(LT::sup(site,mu),mu)], htrans(env.getAdjointLattice()[LT::sup(site,mu)][mu])*htrans(env.getAdjointLattice()[site][mu])*adjointLieGenerator.get(i)*(*scalar_field)[site])*gaugeLieGenerator.get(i);
		//forceLattice[site][mu] += 0.5*vector_dot((*scalar_field)[LT::sdn(site,mu)],env.getAdjointLattice()[LT::sdn(site,mu)][mu]*adjointLieGenerator.get(i)*env.getAdjointLattice()[site][mu]*(*scalar_field)[LT::sup(site,mu)])*gaugeLieGenerator.get(i);
		//forceLattice[site][mu] += 0.5*vector_dot((*scalar_field)[site], adjointLieGenerator.get(i)*env.getAdjointLattice()[site][mu]*env.getAdjointLattice()[LT::sup(site,mu)][mu]*(*scalar_field)[LT::sup(LT::sup(site,mu),mu)])*gaugeLieGenerator.get(i);	
					}
			}
		}
	}

	std::vector<extended_color_vector_t>::const_iterator fundamental_scalar_field;

	for (fundamental_scalar_field = env.fundamental_scalar_fields.begin(); fundamental_scalar_field < env.fundamental_scalar_fields.end(); ++fundamental_scalar_field) {
		for (unsigned int mu = 0; mu < 4; ++mu) {

#pragma omp parallel for
			for (int site = 0; site < fundamental_scalar_field->localsize; ++site) {
					//For every generator
					for (unsigned int i = 0; i < gaugeLieGenerator.numberGenerators(); ++i) {
		forceLattice[site][mu] += 0.5*vector_dot((*fundamental_scalar_field)[LT::sup(site,mu)], htrans(env.getFundamentalLattice()[site][mu])*gaugeLieGenerator.get(i)*(*fundamental_scalar_field)[site])*gaugeLieGenerator.get(i);
		forceLattice[site][mu] -= 0.5*vector_dot((*fundamental_scalar_field)[site], gaugeLieGenerator.get(i)*env.getFundamentalLattice()[site][mu]*(*fundamental_scalar_field)[LT::sup(site,mu)])*gaugeLieGenerator.get(i);
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

GaugeGroup ScalarAction::force(const environment_t& env, int site, int mu) const {
	typedef extended_adjoint_real_color_vector_t LT;		

	GaugeGroup result;
	set_to_zero(result);

	std::vector<extended_adjoint_real_color_vector_t>::const_iterator adjoint_scalar_field;
	for (adjoint_scalar_field = env.adjoint_scalar_fields.begin(); adjoint_scalar_field < env.adjoint_scalar_fields.end(); ++adjoint_scalar_field) {
		//For every generator
		for (unsigned int i = 0; i < adjointLieGenerator.numberGenerators(); ++i) {
			result += 0.5*vector_dot((*adjoint_scalar_field)[LT::sup(site,mu)], htrans(env.getAdjointLattice()[site][mu])*adjointLieGenerator.get(i)*(*adjoint_scalar_field)[site])*gaugeLieGenerator.get(i);
            result -= 0.5*vector_dot((*adjoint_scalar_field)[site], adjointLieGenerator.get(i)*env.getAdjointLattice()[site][mu]*(*adjoint_scalar_field)[LT::sup(site,mu)])*gaugeLieGenerator.get(i);
			/*result -= 0.5*vector_dot((*scalar_field)[LT::sup(site,mu)], htrans(env.getAdjointLattice()[site][mu])*adjointLieGenerator.get(i)*htrans(env.getAdjointLattice()[LT::sdn(site,mu)][mu])*(*scalar_field)[LT::sdn(site,mu)])*gaugeLieGenerator.get(i);
                	result -= 0.5*vector_dot((*scalar_field)[LT::sup(LT::sup(site,mu),mu)], htrans(env.getAdjointLattice()[LT::sup(site,mu)][mu])*htrans(env.getAdjointLattice()[site][mu])*adjointLieGenerator.get(i)*(*scalar_field)[site])*gaugeLieGenerator.get(i);
                	result += 0.5*vector_dot((*scalar_field)[LT::sdn(site,mu)],env.getAdjointLattice()[LT::sdn(site,mu)][mu]*adjointLieGenerator.get(i)*env.getAdjointLattice()[site][mu]*(*scalar_field)[LT::sup(site,mu)])*gaugeLieGenerator.get(i);
                	result += 0.5*vector_dot((*scalar_field)[site], adjointLieGenerator.get(i)*env.getAdjointLattice()[site][mu]*env.getAdjointLattice()[LT::sup(site,mu)][mu]*(*scalar_field)[LT::sup(LT::sup(site,mu),mu)])*gaugeLieGenerator.get(i);*/
        }
	}

	std::vector<extended_color_vector_t>::const_iterator fundamental_scalar_field;
	for (fundamental_scalar_field = env.fundamental_scalar_fields.begin(); fundamental_scalar_field < env.fundamental_scalar_fields.end(); ++fundamental_scalar_field) {
		//For every generator
		for (unsigned int i = 0; i < gaugeLieGenerator.numberGenerators(); ++i) {
			result += 0.5*vector_dot((*fundamental_scalar_field)[LT::sup(site,mu)], htrans(env.getFundamentalLattice()[site][mu])*gaugeLieGenerator.get(i)*(*fundamental_scalar_field)[site])*gaugeLieGenerator.get(i);
            result -= 0.5*vector_dot((*fundamental_scalar_field)[site], gaugeLieGenerator.get(i)*env.getFundamentalLattice()[site][mu]*(*fundamental_scalar_field)[LT::sup(site,mu)])*gaugeLieGenerator.get(i);
			/*result -= 0.5*vector_dot((*scalar_field)[LT::sup(site,mu)], htrans(env.getFundamentalLattice()[site][mu])*fermionLieGenerator.get(i)*htrans(env.getFundamentalLattice()[LT::sdn(site,mu)][mu])*(*scalar_field)[LT::sdn(site,mu)])*gaugeLieGenerator.get(i);
                	result -= 0.5*vector_dot((*scalar_field)[LT::sup(LT::sup(site,mu),mu)], htrans(env.getFundamentalLattice()[LT::sup(site,mu)][mu])*htrans(env.getFundamentalLattice()[site][mu])*fermionLieGenerator.get(i)*(*scalar_field)[site])*gaugeLieGenerator.get(i);
                	result += 0.5*vector_dot((*scalar_field)[LT::sdn(site,mu)],env.getFundamentalLattice()[LT::sdn(site,mu)][mu]*fermionLieGenerator.get(i)*env.getFundamentalLattice()[site][mu]*(*scalar_field)[LT::sup(site,mu)])*gaugeLieGenerator.get(i);
                	result += 0.5*vector_dot((*scalar_field)[site], fermionLieGenerator.get(i)*env.getFundamentalLattice()[site][mu]*env.getFundamentalLattice()[LT::sup(site,mu)][mu]*(*scalar_field)[LT::sup(LT::sup(site,mu),mu)])*gaugeLieGenerator.get(i);*/
        }
	}

	return result;
}

/*void AdjointScalarAction::updateForce(extended_adjoint_real_color_vector_t& force, const environment_t& env) const {
	typedef extended_adjoint_color_vector_t LT;

	std::vector<extended_adjoint_real_color_vector_t>::const_iterator scalar_field;
	for (scalar_field = env.adjoint_scalar_fields.begin(); scalar_field < env.adjoint_scalar_fields.end(); ++scalar_field) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			extended_adjoint_real_color_vector_t dmu;

#pragma omp parallel for
			for (int site = 0; site < dmu.localsize; ++site) {
				dmu[site] = (env.getAdjointLattice()[site][mu]*(*scalar_field)[LT::sup(site,mu)] - htrans(env.getAdjointLattice()[LT::sdn(site,mu)][mu])*(*scalar_field)[LT::sdn(site,mu)]);
			}

			dmu.updateHalo();


#pragma omp parallel for
			for (int site = 0; site < dmu.localsize; ++site) {
				force[site] += env.getAdjointLattice()[site][mu]*dmu[LT::sup(site,mu)];
				force[site] -= htrans(env.getAdjointLattice()[LT::sdn(site,mu)][mu])*dmu[LT::sdn(site,mu)];
			}
		}
	 
#pragma omp parallel for
		for (int site = 0; site < force.localsize; ++site) {
			real_t square = real(vector_dot((*scalar_field)[site],(*scalar_field)[site]));
			for (int c = 0; c < numberColors*numberColors - 1; ++c) {
					force[site][c] += m*(*scalar_field)[site][c]+(lambda/6)*(*scalar_field)[site][c]*square;
					force[site][c] += m*conj((*scalar_field)[site][c])+(lambda/6)*conj((*scalar_field)[site][c])*square;
			}
		}
	}

	force.updateHalo();
}*/

}

