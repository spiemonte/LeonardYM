/*
 * Eigenvalues.cpp
 *
 *  Created on: Jul 16, 2012
 *      Author: spiem_01
 */

#include "Eigenvalues.h"
#include "DiracEigenSolver.h"
#include "io/GlobalOutput.h"
#include "algebra_utils/AlgebraUtils.h"
#include "dirac_functions/Polynomial.h"
#include "utils/ToString.h"

namespace Update {

Eigenvalues::Eigenvalues() : LatticeSweep(), diracEigenSolver(0) { }

Eigenvalues::Eigenvalues(const Eigenvalues& toCopy) : LatticeSweep(toCopy), diracEigenSolver(0) { }

Eigenvalues::~Eigenvalues() {
	if (diracEigenSolver != 0) delete diracEigenSolver;
}

void Eigenvalues::execute(environment_t& environment) {
	typedef reduced_dirac_vector_t::Layout Layout;
	//Take the Dirac Operator
	DiracOperator* squareDiracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
	squareDiracOperator->setLattice(environment.getFermionLattice());

	if (diracEigenSolver == 0)  diracEigenSolver = new DiracEigenSolver();
	diracEigenSolver->setPrecision(environment.configurations.get<double>("generic_inverter_precision"));
	diracEigenSolver->setExtraSteps(environment.configurations.get<unsigned int>("number_extra_vectors_eigensolver"));

	std::vector< std::complex<real_t> > computed_eigenvalues;
	std::vector< reduced_dirac_vector_t > computed_eigenvectors;

	diracEigenSolver->maximumEigenvalues(squareDiracOperator, computed_eigenvalues, computed_eigenvectors, environment.configurations.get<unsigned int>("number_eigenvalues"), LargestReal);
	if (isOutputProcess()) std::cout << "Eigenvalues::Maximal Eigenvalue of square hermitian: " << computed_eigenvalues.front() << std::endl;

	if (isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("maximal_eigenvalues_square_dirac_operator");

		for (unsigned int i = 0; i < computed_eigenvalues.size(); ++i) {
			output->write("maximal_eigenvalues_square_dirac_operator", real(computed_eigenvalues[i]));
		}

		output->pop("maximal_eigenvalues_square_dirac_operator");
	}
/*
	std::vector< std::complex<real_t> > coeff = environment.configurations.get< std::vector< std::complex<real_t> > >("eigenvalues_map");
	Polynomial map;
	map.setScaling(coeff.front());
	coeff.erase(coeff.begin(),coeff.begin()+1);
	map.setRoots(coeff);

	diracEigenSolver->minimumEigenvalues(diracOperator, computed_eigenvalues, computed_eigenvectors, map, 10, 0);
	if (isOutputProcess()) std::cout << "Eigenvalues::Minimal Eigenvalue of square hermitian: " << computed_eigenvalues.front() << std::endl;
	if (isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("minimal_eigenvalues");

		for (int i = 0; i < computed_eigenvalues.size(); ++i) {
			output->write("minimal_eigenvalues", computed_eigenvalues[i]);
		}

		output->pop("minimal_eigenvalues");
	}
	*/
	//std::vector< std::complex<real_t> > computed_eigenvalues;
	//std::vector< reduced_dirac_vector_t > computed_eigenvectors;
	diracEigenSolver->setExtraSteps(30);//TODO TODO TODO
	diracEigenSolver->minimumEigenvalues(squareDiracOperator, computed_eigenvalues, computed_eigenvectors, environment.configurations.get<unsigned int>("number_eigenvalues"));
	if (isOutputProcess()) std::cout << "Minimal Eigenvalue of square hermitian: " << computed_eigenvalues.front() << std::endl;

	if (isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("minimal_eigenvalues_square_dirac_operator");

		std::vector< std::complex<real_t> >::iterator it;

		for (it = computed_eigenvalues.begin(); it != computed_eigenvalues.end(); ++it) {
			output->write("minimal_eigenvalues_square_dirac_operator", real(*it));
		}

		output->pop("minimal_eigenvalues_square_dirac_operator");
	}

	DiracOperator* diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	diracOperator->setLattice(environment.getFermionLattice());
	
	if (isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("minimal_eigenvalues_dirac_operator");
	}

	reduced_dirac_vector_t tmp;
	for (unsigned int i = 0; i < computed_eigenvectors.size(); ++i) {
		diracOperator->multiply(tmp, computed_eigenvectors[i]);
		std::complex<real_t> eigenvalue = static_cast< std::complex<real_t> >(AlgebraUtils::dot(computed_eigenvectors[i],tmp));

		if (isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();
			output->write("minimal_eigenvalues_dirac_operator", real(eigenvalue));
		}
	}

	if (isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->pop("minimal_eigenvalues_dirac_operator");
	}

	delete diracOperator;
	delete squareDiracOperator;

}

} /* namespace Update */
