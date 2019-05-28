#include "GluinoGlue.h"
#include "utils/LieGenerators.h"
#include "utils/Gamma.h"
#include "utils/StoutSmearing.h"
#include "io/GlobalOutput.h"
#include "algebra_utils/AlgebraUtils.h"
#include "utils/MultiThreadSummator.h"

namespace Update {

GluinoGlue::GluinoGlue() : LatticeSweep(), StochasticEstimator(), diracOperator(0), biConjugateGradient(0) { }

GluinoGlue::GluinoGlue(const GluinoGlue& toCopy) : LatticeSweep(toCopy), StochasticEstimator(toCopy), diracOperator(0), biConjugateGradient(0) { }

GluinoGlue::~GluinoGlue() {
	if (diracOperator) delete diracOperator;
	if (biConjugateGradient) delete biConjugateGradient;
}

void GluinoGlue::execute(environment_t& environment) {
	typedef extended_dirac_vector_t::Layout Layout;
	LieGenerator<GaugeGroup> tau;

	if (diracOperator == 0) {
		diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	}

	extended_gauge_lattice_t lattice;

	try {
		unsigned int numberLevelSmearing = environment.configurations.get<unsigned int>("GluinoGlue::level_stout_smearing");
		double smearingRho = environment.configurations.get<double>("GluinoGlue::rho_stout_smearing");
		StoutSmearing stoutSmearing;
		stoutSmearing.spatialSmearing(environment.gaugeLinkConfiguration, lattice, numberLevelSmearing, smearingRho);
	} catch (NotFoundOption& ex) {
		if (isOutputProcess()) std::cout << "GluinoGlue::No smearing options found, proceeding without!" << std::endl;
	}

	extended_fermion_lattice_t result = environment.getFermionLattice();

	try {
		unsigned int t = environment.configurations.get<unsigned int>("GluinoGlue::t_source_origin");
		extended_fermion_lattice_t swaplinkconfig;
		typedef extended_fermion_lattice_t LT;
		if (t != 0) {
			for (unsigned int n = 0; n < t; ++n) {
				//We do a swap
				for(int site = 0; site < (result.localsize); ++site){
					for (unsigned int mu = 0; mu < 4; ++mu) swaplinkconfig[site][mu] = result[site][mu];
				}
				swaplinkconfig.updateHalo();
				//We wrap
				for(int site = 0; site < (result.localsize); ++site){
					for (unsigned int mu = 0; mu < 4; ++mu) result[site][mu] = swaplinkconfig[LT::sup(site,3)][mu];
				}
				result.updateHalo();
			}
		}
	} catch (NotFoundOption& ex) {
	}

	diracOperator->setLattice(result);
	diracOperator->setGamma5(false);

	if (biConjugateGradient == 0) biConjugateGradient = new BiConjugateGradient();
	biConjugateGradient->setPrecision(environment.configurations.get<double>("GluinoGlue::inverter_precision"));
	biConjugateGradient->setMaximumSteps(environment.configurations.get<unsigned int>("GluinoGlue::inverter_max_steps"));
	
	MultiThreadSummator<long_real_t>* gluinoGlueCorrelator = new MultiThreadSummator<long_real_t>[Layout::glob_t];
	for (int t = 0; t < Layout::glob_t; ++t) {
		gluinoGlueCorrelator[t].reset();
	}

	
	for (unsigned int alpha = 0; alpha < 4; ++alpha) {
		//First we generate the source
#pragma omp parallel for
		for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (int c = 0; c < diracVectorLength; ++c) {
					source[site][mu][c] = 0.;
					if (Layout::globalIndexT(site) == 0 && Layout::globalIndexX(site) == 0 && Layout::globalIndexY(site) == 0 && Layout::globalIndexZ(site) == 0) {
						for (unsigned int i = 0; i < diracVectorLength; ++i) {
							for (unsigned int j = 0; j < diracVectorLength; ++j) {
								if (Sigma::sigma(i,j,alpha,mu) != static_cast<real_t>(0.)) source[site][mu][c] += Sigma::sigma(i,j,alpha,mu)*trace(cloverPlaquette(lattice,site,i,j)*tau.get(c));
							}
						}
					}
				}
			}
		}

		biConjugateGradient->solve(diracOperator, source, psi);

		for (int t = 0; t < Layout::glob_t; ++t) {
#pragma omp parallel for
			for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int c = 0; c < diracVectorLength; ++c) {
						eta[site][mu][c] = 0.;
						if (Layout::globalIndexT(site) == t) {
							for (unsigned int i = 0; i < diracVectorLength; ++i) {
								for (unsigned int j = 0; j < diracVectorLength; ++j) {
									if (Sigma::sigma(i,j,alpha,mu) != static_cast<real_t>(0.)) eta[site][mu][c] += Sigma::sigma(i,j,alpha,mu)*trace(cloverPlaquette(lattice,site,i,j)*tau.get(c));
								}
							}
						}
					}
				}
			}

#pragma omp parallel for
			for (int site = 0; site < Layout::localsize; ++site) {
				if (Layout::globalIndexT(site) == 0) {
					for (unsigned int beta = 0; beta < 4; ++beta) {
						for (int a = 0; a < diracVectorLength; ++a) {
							gluinoGlueCorrelator[t].add(real(eta[site][beta][a]*psi[site][beta][a]));
						}
					}
				}
			}
		}
	}

	for (int t = 0; t < Layout::glob_t; ++t) {
		gluinoGlueCorrelator[t].computeResult();
	}

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();

		output->push("gluinoglue");
		for (int t = 0; t < Layout::glob_t; ++t) {
			std::cout << "GluinoGlue::Correlator at t " << t << " is " << -gluinoGlueCorrelator[t].getResult() << std::endl;

			output->write("gluinoglue", -gluinoGlueCorrelator[t].getResult());
		}
		output->pop("gluinoglue");
	}

	delete[] gluinoGlueCorrelator;
	
}

GaugeGroup GluinoGlue::cloverPlaquette(const extended_gauge_lattice_t& lattice, int site, int mu, int nu) {
	typedef extended_dirac_vector_t LT;
	
	GaugeGroup result;
	set_to_zero(result);
	if (mu != nu) {
		result += lattice[site][mu]*lattice[LT::sup(site,mu)][nu]*htrans(lattice[LT::sup(site,nu)][mu])*htrans(lattice[site][nu]);
		result += htrans(lattice[LT::sdn(site,nu)][nu])*lattice[LT::sdn(site,nu)][mu]*lattice[LT::sdn(LT::sup(site,mu),nu)][nu]*htrans(lattice[site][mu]);
		result += lattice[site][nu]*htrans(lattice[LT::sdn(LT::sup(site,nu),mu)][mu])*htrans(lattice[LT::sdn(site,mu)][nu])*lattice[LT::sdn(site,mu)][mu];
		result += htrans(lattice[LT::sdn(site,mu)][mu])*htrans(lattice[LT::sdn(LT::sdn(site,mu),nu)][nu])*lattice[LT::sdn(LT::sdn(site,mu),nu)][mu]*lattice[LT::sdn(site,nu)][nu];
	}
	GaugeGroup hresult = htrans(result);
	GaugeGroup cplaq = std::complex<real_t>(0,+1./8.)*(result-hresult);
	return cplaq;
}

void GluinoGlue::registerParameters(po::options_description& desc) {
	static bool single = true;
	if (single) desc.add_options()
		("GluinoGlue::inverter_precision", po::value<real_t>()->default_value(0.0000000001), "set the inverter precision")
		("GluinoGlue::inverter_max_steps", po::value<unsigned int>()->default_value(10000), "maximum number of inverter steps")
		("GluinoGlue::t_source_origin", po::value<unsigned int>()->default_value(0), "T origin for the wall source")
		("GluinoGlue::rho_stout_smearing", po::value<real_t>(), "set the stout smearing parameter")
		("GluinoGlue::levels_stout_smearing", po::value<unsigned int>(), "levels of stout smearing")
		;
	single = false;
}

}

