#include "Glueball.h"
#include "io/GlobalOutput.h"
#include "utils/ToString.h"
#include "utils/StoutSmearing.h"
#include "utils/MultiThreadSummator.h"

namespace Update {

Glueball::Glueball() { }

Glueball::~Glueball() { }

void Glueball::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t::Layout Layout;
	typedef extended_gauge_lattice_t LT;

	extended_gauge_lattice_t lattice;
	try {
		unsigned int numberLevelSmearing = environment.configurations.get<unsigned int>("Glueball::stout_smearing_levels");
		double smearingRho = environment.configurations.get<double>("Glueball::stout_smearing_rho");
		extended_gauge_lattice_t smearedConfiguration;
		StoutSmearing stoutSmearing;
		stoutSmearing.spatialSmearing(environment.gaugeLinkConfiguration, lattice, numberLevelSmearing, smearingRho);
	} catch (NotFoundOption& ex) {
		lattice =  environment.gaugeLinkConfiguration;
		if (isOutputProcess()) std::cout << "Glueball::No smearing options found, proceeding without!" << std::endl;
	}

	//We measure the zero momentum projection
	MultiThreadSummator<long_real_t>* zero_glueball = new MultiThreadSummator<long_real_t>[Layout::glob_t];
	MultiThreadSummator<long_real_t>* two_glueball = new MultiThreadSummator<long_real_t>[Layout::glob_t];

	for (int t = 0; t < Layout::glob_t; ++t) {
		zero_glueball[t].reset();
		two_glueball[t].reset();
	}

#pragma omp parallel for
	for (int site = 0; site < Layout::localsize; ++site) {
		long_real_t result_phi_x = real(trace(lattice[site][1]*lattice[LT::sup(site,1)][0]*htrans(lattice[LT::sup(site,0)][1])*htrans(lattice[site][0])));
		long_real_t result_phi_y = real(trace(lattice[site][2]*lattice[LT::sup(site,2)][0]*htrans(lattice[LT::sup(site,0)][2])*htrans(lattice[site][0])));
		long_real_t result_phi_z = real(trace(lattice[site][2]*lattice[LT::sup(site,2)][1]*htrans(lattice[LT::sup(site,1)][2])*htrans(lattice[site][1])));
		
		two_glueball[Layout::globalIndexT(site)].add(result_phi_x - result_phi_y);
		zero_glueball[Layout::globalIndexT(site)].add(result_phi_x + result_phi_y + result_phi_z);
	}

	//We collect the results
	for (int t = 0; t < Layout::glob_t; ++t) {
		zero_glueball[t].computeResult();
		two_glueball[t].computeResult();
	}

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();

		output->push("glueball_zero");
		for (int t = 0; t < Layout::glob_t; ++t) {
			std::cout << "Glueball::0++ Operator at t " << t << " is " << zero_glueball[t].getResult()/Layout::glob_spatial_volume << std::endl;

			output->write("glueball_zero", zero_glueball[t].getResult()/Layout::glob_spatial_volume);
		}
		output->pop("glueball_zero");

		output->push("glueball_two");
		for (int t = 0; t < Layout::glob_t; ++t) {
			std::cout << "Glueball::2++ Operator at t " << t << " is " << two_glueball[t].getResult()/Layout::glob_spatial_volume << std::endl;

			output->write("glueball_two", two_glueball[t].getResult()/Layout::glob_spatial_volume);
		}
		output->pop("glueball_two");
	}

	delete[] zero_glueball;
	delete[] two_glueball;
}

void Glueball::registerParameters(po::options_description& desc) {
	desc.add_options()
		("Glueball::stout_smearing_rho", po::value<real_t>()->default_value(0.15), "set the stout smearing parameter")
		("Glueball::stout_smearing_levels", po::value<unsigned int>()->default_value(10), "levels of stout smearing")
	;
}

} /* namespace Update */
