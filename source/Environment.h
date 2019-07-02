/*
 * Enviroment.h
 *
 *  Created on: Feb 21, 2012
 *      Author: spiem_01
 */

#ifndef ENVIRONMENT_H_
#define ENVIRONMENT_H_

#include "LatticeDefinitions.h"
#include "LatticeBoundaryConditions.h"
#include "MPILattice/MPIUtils.h"

#include "io/StorageParameters.h"
#include "utils/ConvertLattice.h"

namespace Update {

//This class stores the gauge, fermion and scalar field configurations being simulated by the Monte-Carlo updaters
//We also store the configurations and the couplings of the theory
class Environment {
public:
	Environment() : gaugeLinkConfiguration(), fermionicLinkConfiguration(), adjointLinkConfiguration(), sweep(0), iteration(0), measurement(false) { }
	Environment(const boost::program_options::variables_map& vm) : gaugeLinkConfiguration(), fermionicLinkConfiguration(), adjointLinkConfiguration(), configurations(vm), sweep(0), iteration(0), measurement(false) { }
	Environment(const Environment& toCopy);
	Environment& operator=(const Environment& toCopy);
	~Environment() { }

	//This function ensures that all gauge configurations in all representations are correctly matching
	void synchronize();

	//Set the fermion boundary conditions
	template<typename T> void setFermionBc(T& lattice) const {
		try {
			std::string bc = configurations.get<std::string>("boundary_conditions");
			if (bc == "fermion_antiperiodic") {
				switchAntiperiodicBc(lattice);
			}
			else if (bc == "fermion_spatial_antiperiodic") {
				switchSpatialAntiperiodicBc(lattice);
			}
			else if (bc == "fermion_full_antiperiodic") {
				switchFullAntiperiodicBc(lattice);
			}
			else if (bc == "fermion_periodic") {

			}
			else if (bc == "open") {
				switchOpenBc(lattice);
			}
			else {
				static int count = 0;
				if (count == 0 && isOutputProcess()) {
					std::cout << "Warning, boundary conditions not or badly set, using antiperiodic!" << std::endl;
					count = count + 1;
				}
				switchAntiperiodicBc(lattice);
			}
		}
		catch (NotFoundOption& ex) {
			static int count = 0;
			if (count == 0 && isOutputProcess()) {
				std::cout << "Warning, boundary conditions not or badly set, using antiperiodic!" << std::endl;
				count = count + 1;
			}
			switchAntiperiodicBc(lattice);
		}
	}

	//Get the gauge field interacting with the fermions
	extended_fermion_lattice_t& getFermionLattice() {
		return fermionicLinkConfiguration;
	}

	const extended_fermion_lattice_t& getFermionLattice() const {
		return fermionicLinkConfiguration;
	}

	//Get the gauge field in the adjoint representation of the gauge group
	extended_adjoint_lattice_t& getAdjointLattice() {
                return adjointLinkConfiguration;
        }

        const extended_adjoint_lattice_t& getAdjointLattice() const {
                return adjointLinkConfiguration;
        }

	//Get the gauge field in the fundamental representation of the gauge group
	extended_gauge_lattice_t& getFundamentalLattice() {
		return gaugeLinkConfiguration;
	}

	const extended_gauge_lattice_t& getFundamentalLattice() const {
		return gaugeLinkConfiguration;
	}

	//Gauge link configuration
	extended_gauge_lattice_t gaugeLinkConfiguration;
	//Gauge link configuration interacting with fermions
	extended_fermion_lattice_t fermionicLinkConfiguration;
	//Gauge link configuration in the adjoint representation
	extended_adjoint_lattice_t adjointLinkConfiguration;

	std::vector<extended_adjoint_color_vector_t> adjoint_scalar_fields;
	std::vector<extended_color_vector_t> fundamental_scalar_fields;

	//The configurations of the main program
	StorageParameters configurations;
	
	unsigned int sweep;
	unsigned int iteration;

	bool measurement;
};

typedef Environment environment_t;

}

#endif /* ENVIROMENT_H_ */
