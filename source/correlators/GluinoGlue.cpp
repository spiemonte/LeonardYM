#include "GluinoGlue.h"
#include "utils/LieGenerators.h"
#include "utils/Gamma.h"
#include "utils/StoutSmearing.h"
#include "io/GlobalOutput.h"
#include "algebra_utils/AlgebraUtils.h"
#include "utils/MultiThreadSummator.h"
#include "dirac_operators/Propagator.h"
#include "utils/RandomGaugeTransformation.h"
#include "wilson_loops/Plaquette.h"

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
	extended_fermion_lattice_t smearedFermionLattice = environment.getFermionLattice();

	try {
		unsigned int numberLevelSmearing = environment.configurations.get<unsigned int>("GluinoGlue::stout_smearing_levels");
		double smearingRho = environment.configurations.get<double>("GluinoGlue::stout_smearing_rho");
		StoutSmearing stoutSmearing;
		stoutSmearing.spatialSmearing(environment.gaugeLinkConfiguration, lattice, numberLevelSmearing, smearingRho);
	} catch (NotFoundOption& ex) {
		if (isOutputProcess()) std::cout << "GluinoGlue::No smearing options found, proceeding without!" << std::endl;
		lattice = environment.gaugeLinkConfiguration;
	}

	extended_fermion_lattice_t fermionLattice = environment.getFermionLattice();

	try {
		unsigned int t = environment.configurations.get<unsigned int>("GluinoGlue::t_source_origin");
		
		if (t != 0) {
			extended_fermion_lattice_t swaplinkconfig;
			typedef extended_fermion_lattice_t LT;
			for (unsigned int n = 0; n < t; ++n) {
				//We do a swap
				for(int site = 0; site < (fermionLattice.localsize); ++site){
					for (unsigned int mu = 0; mu < 4; ++mu) swaplinkconfig[site][mu] = fermionLattice[site][mu];
				}
				swaplinkconfig.updateHalo();
				//We wrap
				for(int site = 0; site < (fermionLattice.localsize); ++site){
					for (unsigned int mu = 0; mu < 4; ++mu) fermionLattice[site][mu] = swaplinkconfig[LT::sup(site,3)][mu];
				}
				fermionLattice.updateHalo();
			}
		}
		
		if (t != 0) {
			extended_gauge_lattice_t swaplinkconfig;
			typedef extended_gauge_lattice_t LT;
			for (unsigned int n = 0; n < t; ++n) {
				//We do a swap
				for(int site = 0; site < (lattice.localsize); ++site){
					for (unsigned int mu = 0; mu < 4; ++mu) swaplinkconfig[site][mu] = lattice[site][mu];
				}
				swaplinkconfig.updateHalo();
				//We wrap
				for(int site = 0; site < (lattice.localsize); ++site){
					for (unsigned int mu = 0; mu < 4; ++mu) lattice[site][mu] = swaplinkconfig[LT::sup(site,3)][mu];
				}
				lattice.updateHalo();
			}
		}
	} catch (NotFoundOption& ex) {
	}

	diracOperator->setLattice(fermionLattice);
	diracOperator->setGamma5(false);

	if (biConjugateGradient == 0) biConjugateGradient = new BiConjugateGradient();
	biConjugateGradient->setPrecision(environment.configurations.get<double>("GluinoGlue::inverter_precision"));
	biConjugateGradient->setMaximumSteps(environment.configurations.get<unsigned int>("GluinoGlue::inverter_max_steps"));

	real_t source_smearing_rho = environment.configurations.get<double>("GluinoGlue::source_smearing_rho");
	unsigned int source_smearing_levels = environment.configurations.get<unsigned int>("GluinoGlue::source_smearing_levels");
	
	MultiThreadSummator<long_real_t>* gluinoGlueCorrelator = new MultiThreadSummator<long_real_t>[Layout::glob_t];
	for (int t = 0; t < Layout::glob_t; ++t) {
		gluinoGlueCorrelator[t].reset();
	}

	extended_dirac_vector_t tmp;

	int inversionSteps = 0;

	//First we generate the propagator
	for (unsigned int alpha = 0; alpha < 4; ++alpha) {			
		//First we generate the source
#pragma omp parallel for
		for (int site = 0; site < source.localsize; ++site) {
			for (unsigned int rho = 0; rho < 4; ++rho) {
				set_to_zero(source[site][rho]);
				if (Layout::globalIndexT(site) == 0) {
					for (unsigned int k = 0; k < 3; ++k) {
						for (unsigned int l = 0; l < 3; ++l) {
							if (Sigma::sigma(k,l,alpha,rho) != static_cast<real_t>(0.)) {
								for (int b = 0; b < diracVectorLength; ++b) {
									source[site][rho][b] += Sigma::sigma(k,l,alpha,rho)*trace(cloverPlaquette(lattice,site,k,l)*tau.get(b));
								}
							}
						}
					}
				}
			}
		}

		this->smearSource(source, smearedFermionLattice, source_smearing_levels, source_smearing_rho);

	
		biConjugateGradient->solve(diracOperator, source, tmp);
		Propagator::constructPropagator(diracOperator, tmp, psi[alpha]);

		this->smearSource(psi[alpha], smearedFermionLattice, source_smearing_levels, source_smearing_rho);

		inversionSteps += biConjugateGradient->getLastSteps();
	}

	if (isOutputProcess()) std::cout << "GluinoGlue::Correlators computed with " << inversionSteps << " inversion steps" << std::endl;

	//Now we compute the contraction at the sink
#pragma omp parallel for
	for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
		for (unsigned int alpha = 0; alpha < 4; ++alpha) {
			for (unsigned int beta = 0; beta < 4; ++beta) {
				for (unsigned int i = 0; i < 3; ++i) {
					for (unsigned int j = 0; j < 3; ++j) {
						if (Sigma::sigma(i,j,alpha,beta) != static_cast<real_t>(0.)) {
							for (int a = 0; a < diracVectorLength; ++a) {
								gluinoGlueCorrelator[Layout::globalIndexT(site)].add(real( Sigma::sigma(i,j,alpha,beta)*trace(cloverPlaquette(lattice,site,i,j)*tau.get(a))*psi[alpha][site][beta][a] ));
							}
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
			std::cout << "GluinoGlue::Correlator at t " << t << " is " << -gluinoGlueCorrelator[t].getResult()/4. << std::endl;

			output->write("gluinoglue", -gluinoGlueCorrelator[t].getResult()/4.);
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
	GaugeGroup cplaq = std::complex<real_t>(0.,-1./8.)*(result-hresult);
	return cplaq;
}

void GluinoGlue::registerParameters(po::options_description& desc) {
	static bool single = true;
	if (single) desc.add_options()
		("GluinoGlue::inverter_precision", po::value<real_t>()->default_value(0.000000000001), "set the inverter precision")
		("GluinoGlue::inverter_max_steps", po::value<unsigned int>()->default_value(10000), "maximum number of inverter steps")
		("GluinoGlue::t_source_origin", po::value<unsigned int>()->default_value(0), "T origin for the wall source")
		("GluinoGlue::stout_smearing_rho", po::value<real_t>()->default_value(0.15), "set the stout smearing parameter")
		("GluinoGlue::stout_smearing_levels", po::value<unsigned int>()->default_value(15), "levels of stout smearing")
		("GluinoGlue::source_smearing_rho", po::value<real_t>()->default_value(0.15), "set the Jacoby smearing rho for the source/sink")
		("GluinoGlue::source_smearing_levels", po::value<unsigned int>()->default_value(5), "set the Jacoby smearing levels for the source/sink")
		;
	single = false;
}

}

