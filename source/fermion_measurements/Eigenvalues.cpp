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

	if (diracEigenSolver == 0)  diracEigenSolver = new DiracEigenSolver();
	diracEigenSolver->setInverterPrecision(environment.configurations.get<double>("Eigenvalues::inverter_precision"));
	diracEigenSolver->setTolerance(environment.configurations.get<double>("Eigenvalues::eigensolver_precision"));

	int restarts = environment.configurations.get<unsigned int>("Eigenvalues::maximal_number_restarts_eigensolver");
	diracEigenSolver->setMaximalNumberOfRestarts(restarts);

	diracEigenSolver->setExtraSteps(environment.configurations.get<unsigned int>("Eigenvalues::number_extra_vectors"));	
	
	std::string hermitian = environment.configurations.get<std::string>("Eigenvalues::hermitian");

	std::string mode = environment.configurations.get<std::string>("Eigenvalues::mode");
	
	std::vector< std::complex<real_t> > computed_eigenvalues;
	std::vector< reduced_dirac_vector_t > computed_eigenvectors;

	if (environment.configurations.get<std::string>("Eigenvalues::use_chebyshev") == "true" && mode == "LR") {
		unsigned int chebyshevOrder = environment.configurations.get<unsigned int>("Eigenvalues::chebyshev_order");
		double chebyshevLeft = environment.configurations.get<double>("Eigenvalues::chebyshev_left");
		double chebyshevRight = environment.configurations.get<double>("Eigenvalues::chebyshev_right");

		diracEigenSolver->setUseChebyshev(true);
		diracEigenSolver->getChebyshevRecursion()->setParameters(chebyshevLeft, chebyshevRight, chebyshevOrder);
	}


	if (hermitian == "true") {
		//Take the Dirac Operator
		DiracOperator* squareDiracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
		squareDiracOperator->setLattice(environment.getFermionLattice());
		DiracOperator* diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
		diracOperator->setLattice(environment.getFermionLattice());
		
		if (mode == "LR") {
			diracEigenSolver->maximumEigenvalues(squareDiracOperator, computed_eigenvalues, computed_eigenvectors, environment.configurations.get<unsigned int>("Eigenvalues::number_eigenvalues"), LargestReal);
			if (isOutputProcess()) std::cout << "Eigenvalues::Maximal Eigenvalue of square hermitian: " << computed_eigenvalues.front() << std::endl;

			if (isOutputProcess()) {
				GlobalOutput* output = GlobalOutput::getInstance();
				output->push("maximal_eigenvalues_hermitian_square_dirac_operator");
	
				for (unsigned int i = 0; i < computed_eigenvalues.size(); ++i) {
					output->write("maximal_eigenvalues_hermitian_square_dirac_operator", real(computed_eigenvalues[i]));
				}

				output->pop("maximal_eigenvalues_hermitian_square_dirac_operator");
			}

			DiracOperator* diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
			diracOperator->setLattice(environment.getFermionLattice());
	
			if (isOutputProcess()) {
				GlobalOutput* output = GlobalOutput::getInstance();
				output->push("maximal_eigenvalues_hermitian_dirac_operator");
				output->push("maximal_eigenvalues_hermitian_chirality");
			}	
	
			reduced_dirac_vector_t tmp;
			for (unsigned int i = 0; i < computed_eigenvectors.size(); ++i) {
				diracOperator->multiply(tmp, computed_eigenvectors[i]);
				std::complex<real_t> eigenvalue = static_cast< std::complex<real_t> >(AlgebraUtils::dot(computed_eigenvectors[i],tmp));
				std::complex<long_real_t> chirality = AlgebraUtils::gamma5dot(computed_eigenvectors[i],computed_eigenvectors[i]);
			


				if (isOutputProcess()) {
					GlobalOutput* output = GlobalOutput::getInstance();
					output->write("maximal_eigenvalues_hermitian_dirac_operator", real(eigenvalue));
					output->write("maximal_eigenvalues_hermitian_chirality", real(chirality));
					output->write("maximal_eigenvalues_hermitian_chirality", imag(chirality));
				}
			}
			if (isOutputProcess()) {
				GlobalOutput* output = GlobalOutput::getInstance();
				output->pop("maximal_eigenvalues_hermitian_dirac_operator");
				output->pop("maximal_eigenvalues_hermitian_chirality");
			}
		}
		else if (mode == "SR") {
			if (environment.configurations.get<std::string>("Eigenvalues::use_inverse_power_method") == "true") {
				diracEigenSolver->minimumEigenvalues(squareDiracOperator, computed_eigenvalues, computed_eigenvectors, environment.configurations.get<unsigned int>("Eigenvalues::number_eigenvalues"));
			}
			else diracEigenSolver->maximumEigenvalues(squareDiracOperator, computed_eigenvalues, computed_eigenvectors, environment.configurations.get<unsigned int>("Eigenvalues::number_eigenvalues"), SmallestReal);

			if (isOutputProcess()) std::cout << "Minimal Eigenvalue of square hermitian: " << computed_eigenvalues.front() << std::endl;

			if (isOutputProcess()) {
				GlobalOutput* output = GlobalOutput::getInstance();
				output->push("minimal_eigenvalues_hermitian_square_dirac_operator");

				std::vector< std::complex<real_t> >::iterator it;

				for (it = computed_eigenvalues.begin(); it != computed_eigenvalues.end(); ++it) {
					output->write("minimal_eigenvalues_hermitian_square_dirac_operator", real(*it));
				}

				output->pop("minimal_eigenvalues_hermitian_square_dirac_operator");
			}
	
			if (isOutputProcess()) {
				GlobalOutput* output = GlobalOutput::getInstance();
				output->push("minimal_eigenvalues_hermitian_dirac_operator");
				output->push("minimal_eigenvalues_hermitian_chirality");
			}	
	
			reduced_dirac_vector_t tmp;
			for (unsigned int i = 0; i < computed_eigenvectors.size(); ++i) {
				diracOperator->multiply(tmp, computed_eigenvectors[i]);
				std::complex<real_t> eigenvalue = static_cast< std::complex<real_t> >(AlgebraUtils::dot(computed_eigenvectors[i],tmp));
				std::complex<long_real_t> chirality = AlgebraUtils::gamma5dot(computed_eigenvectors[i],computed_eigenvectors[i]);
			


				if (isOutputProcess()) {
					GlobalOutput* output = GlobalOutput::getInstance();
					output->write("minimal_eigenvalues_hermitian_dirac_operator", real(eigenvalue));
					output->write("minimal_eigenvalues_hermitian_chirality", real(chirality));
					output->write("minimal_eigenvalues_hermitian_chirality", imag(chirality));
				}
			}
			if (isOutputProcess()) {
				GlobalOutput* output = GlobalOutput::getInstance();
				output->pop("minimal_eigenvalues_hermitian_dirac_operator");
				output->pop("minimal_eigenvalues_hermitian_chirality");
			}
		}

		delete diracOperator;
		delete squareDiracOperator;
	}
	else {
		DiracOperator* squareHermitianDiracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
		squareHermitianDiracOperator->setLattice(environment.getFermionLattice());
		squareHermitianDiracOperator->setGamma5(true);
		DiracOperator* diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
		diracOperator->setLattice(environment.getFermionLattice());
		diracOperator->setGamma5(false);

		if (mode == "LR") {
			diracEigenSolver->maximumEigenvalues(diracOperator, computed_eigenvalues, computed_eigenvectors, environment.configurations.get<unsigned int>("Eigenvalues::number_eigenvalues"), LargestReal);
			if (isOutputProcess()) std::cout << "Eigenvalues::Maximal Eigenvalue of non-hermitian: " << computed_eigenvalues.front() << std::endl;

			if (isOutputProcess()) {
				GlobalOutput* output = GlobalOutput::getInstance();
				output->push("maximal_eigenvalues_non_hermitian_dirac_operator");

				for (unsigned int i = 0; i < computed_eigenvalues.size(); ++i) {
					output->push("maximal_eigenvalues_non_hermitian_dirac_operator");
					output->write("maximal_eigenvalues_non_hermitian_dirac_operator", real(computed_eigenvalues[i]));
					output->write("maximal_eigenvalues_non_hermitian_dirac_operator", imag(computed_eigenvalues[i]));
					output->pop("maximal_eigenvalues_non_hermitian_dirac_operator");
				}

				output->pop("maximal_eigenvalues_non_hermitian_dirac_operator");
			}
		}
		else if (mode == "SR") {
			if (environment.configurations.get<std::string>("Eigenvalues::use_inverse_power_method") == "true") {
				diracEigenSolver->minimumNonHermitianEigenvalues(diracOperator, squareHermitianDiracOperator, computed_eigenvalues, computed_eigenvectors, environment.configurations.get<unsigned int>("Eigenvalues::number_eigenvalues"));
			}
			else diracEigenSolver->maximumEigenvalues(diracOperator, computed_eigenvalues, computed_eigenvectors, environment.configurations.get<unsigned int>("Eigenvalues::number_eigenvalues"), SmallestReal);

			if (isOutputProcess()) std::cout << "Minimal Eigenvalue of non-hermitian: " << computed_eigenvalues.front() << std::endl;

			if (isOutputProcess()) {
				GlobalOutput* output = GlobalOutput::getInstance();
				output->push("minimal_eigenvalues_non_hermitian_dirac_operator");

				std::vector< std::complex<real_t> >::iterator it;

				for (it = computed_eigenvalues.begin(); it != computed_eigenvalues.end(); ++it) {
					output->push("minimal_eigenvalues_non_hermitian_dirac_operator");
					output->write("minimal_eigenvalues_non_hermitian_dirac_operator", real(*it));
					output->write("minimal_eigenvalues_non_hermitian_dirac_operator", imag(*it));
					output->pop("minimal_eigenvalues_non_hermitian_dirac_operator");
				}

				output->pop("minimal_eigenvalues_non_hermitian_dirac_operator");
			}

			std::vector< std::complex<real_t> > chirality(computed_eigenvectors.size());
			for (unsigned int i = 0; i < computed_eigenvectors.size(); ++i) {
				chirality[i] = AlgebraUtils::gamma5dot(computed_eigenvectors[i], computed_eigenvectors[i]);
			}
		

			if (isOutputProcess()) {
				GlobalOutput* output = GlobalOutput::getInstance();
				output->push("minimal_eigenvalues_non_hermitian_chirality");

				std::vector< std::complex<real_t> >::iterator it;

				for (it = chirality.begin(); it != chirality.end(); ++it) {
					output->push("minimal_eigenvalues_non_hermitian_chirality");
					output->write("minimal_eigenvalues_non_hermitian_chirality", real(*it));
					output->write("minimal_eigenvalues_non_hermitian_chirality", imag(*it));
					output->pop("minimal_eigenvalues_non_hermitian_chirality");
				}

				output->pop("minimal_eigenvalues_non_hermitian_chirality");
			}
		}

		delete diracOperator;
		delete squareHermitianDiracOperator;
	}

}

void Eigenvalues::registerParameters(std::map<std::string, Option>& desc) {
	desc["Eigenvalues::mode"] = Option("Eigenvalues::mode", "LR", "Eigenvalues computed from LargestReal (LR) or SmallestReal (SR)?");
	desc["Eigenvalues::use_chebyshev"] = Option("Eigenvalues::use_chebyshev", "false", "Use Chebyshev acceleration? Option valid only in the LR mode. (true/false)");
	desc["Eigenvalues::chebyshev_left"] = Option("Eigenvalues::chebyshev_left", 0.2, "Left interval of the Chebyshev polynomial.");
	desc["Eigenvalues::chebyshev_right"] = Option("Eigenvalues::chebyshev_right", 7.0, "Right interval of the Chebyshev polynomial.");
	desc["Eigenvalues::chebyshev_order"] = Option("Eigenvalues::chebyshev_order", 15, "Order of the Chebyshev acceleration. It must be an odd number");
	desc["Eigenvalues::hermitian"] = Option("Eigenvalues::hermitian", "true", "Eigenvalues of the Hermitian or Non-Hermitian Dirac operator? (true/false)");
	desc["Eigenvalues::inverter_precision"] = Option("Eigenvalues::inverter_precision", 1e-10, "set the precision used by the inverter");
	desc["Eigenvalues::eigensolver_precision"] = Option("Eigenvalues::eigensolver_precision", 1e-7, "set the precision used by the eigensolver");
	desc["Eigenvalues::number_extra_vectors"] = Option("Eigenvalues::number_extra_vectors", 30, "Number of extra vectors for the Arnoldi algorithm used in the computation of the eigenvectors, increase this number to increase precision");
	desc["Eigenvalues::maximal_number_restarts_eigensolver"] = Option("Eigenvalues::maximal_number_restarts_eigensolver", 50, "Maximal number of restarts for the implicitly restarted Arnoldi algorithm");
	desc["Eigenvalues::number_eigenvalues"] = Option("Eigenvalues::number_eigenvalues", 20, "Number of eigenvalues of the dirac wilson operator to be computed");
	desc["Eigenvalues::use_inverse_power_method"] = Option("Eigenvalues::use_inverse_power_method", "true", "Should we use the inverse power method to compute the minimal eigenvalues");
}

} /* namespace Update */
