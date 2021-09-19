#ifndef LANDAUGHOSTPROPAGATOR_H
#define LANDAUGHOSTPROPAGATOR_H

#include "LatticeSweep.h"
#include "gauge_fixing/LandauGaugeFixing.h"

namespace Update {

class LandauGhostPropagator : public Update::LandauGaugeFixing {
public:
	LandauGhostPropagator();
	LandauGhostPropagator(const LandauGhostPropagator&);
	~LandauGhostPropagator();

	void execute(environment_t& environment);

	static void registerParameters(std::map<std::string, Option>&);
protected:
	void ghostMatrix(extended_adjoint_complex_color_vector_t& output, const extended_adjoint_complex_color_vector_t& input, const extended_adjoint_lattice_t& A, const extended_adjoint_lattice_t& B, const extended_adjoint_lattice_t& C) const;
	bool ghostPropagatorCG(std::complex<real_t>& result, const extended_gauge_lattice_t& lattice, const std::vector<real_t>& momentum, int c, real_t epsilon, unsigned int max_steps) const;
	bool ghostPropagatorBiCGStab(std::complex<real_t>& result, const extended_gauge_lattice_t& lattice, const std::vector<real_t>& momentum, int c, real_t epsilon, unsigned int max_steps) const;
};

}

#endif
