/*
 * Eigenvalues.cpp
 *
 *  Created on: Jul 16, 2012
 *      Author: spiem_01
 */

#include "Eigenvalues.h"
#include "DiracEigenSolver.h"
#include "GlobalOutput.h"
#include "AlgebraUtils.h"
#include "Polynomial.h"
#ifdef HAVE_ARPACK
#include "../math/diracvect/arnoldidiracoperator.h"
#include "../math/diracvect/eigenspacediracoperator.h"
#include "../math/linalg/bindings/arpack/basicarnoldi.h"
#endif

namespace Update {

#ifdef HAVE_ARPACK
class DiracOperatorWrapper : public ::Math::DiracVect::VectorOperator< Update::dirac_vector_t > {
public:
	DiracOperatorWrapper(DiracOperator* _dirac): dirac(_dirac) { }

	virtual void multiply(Vector& In, Vector& Out) {
		dirac->multiply(Out,In);
	}

	virtual void multiplyAdd(dirac_vector_t& In, dirac_vector_t& Out, const dirac_vector_t& Add, const VectorElement& fakt) {
		dirac->multiplyAdd(Out,In,Add,fakt);
	}

	virtual size_t offset() const {
		return 0;
	}

	virtual size_t vectorSize() const {
		return dirac_vector_t::completesize;
	}

private:
	DiracOperator* dirac;
};
#endif

Eigenvalues::Eigenvalues() : LatticeSweep(), diracEigenSolver(0) { }

Eigenvalues::Eigenvalues(const Eigenvalues& toCopy) : LatticeSweep(toCopy), diracEigenSolver(0) { }

Eigenvalues::~Eigenvalues() {
	if (diracEigenSolver != 0) delete diracEigenSolver;
}

void Eigenvalues::execute(environment_t& environment) {
	//Take the Dirac Operator
	DiracOperator* diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
	diracOperator->setLattice(environment.getFermionLattice());

	if (diracEigenSolver == 0)  diracEigenSolver = new DiracEigenSolver();
	diracEigenSolver->setPrecision(environment.configurations.get<double>("generic_inverter_precision"));
	diracEigenSolver->setExtraSteps(50);

	std::vector< std::complex<real_t> > computed_eigenvalues;
	std::vector< reduced_dirac_vector_t > computed_eigenvectors;

	diracEigenSolver->maximumEigenvalues(diracOperator, computed_eigenvalues, computed_eigenvectors, 10);
	if (isOutputProcess()) std::cout << "Eigenvalues::Maximal Eigenvalue of square hermitian: " << computed_eigenvalues.front() << std::endl;

	if (isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("maximal_eigenvalues");

		for (int i = 0; i < computed_eigenvalues.size(); ++i) {
			output->write("maximal_eigenvalues", computed_eigenvalues[i]);
		}

		output->pop("maximal_eigenvalues");
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
	/*std::vector< std::complex<real_t> > min_eigenvalues = diracEigenSolver->minimumEigenvalues(diracOperator, 100);
	if (isOutputProcess()) std::cout << "Minimal Eigenvalue of square hermitian: " << min_eigenvalues.front() << std::endl;*/

	/*if (isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("minimal_eigenvalues");

		std::vector< std::complex<real_t> >::iterator it;

		for (it = min_eigenvalues.begin(); it != min_eigenvalues.end(); ++it) {
			output->write("minimal_eigenvalues", *it);
		}

		output->pop("minimal_eigenvalues");
	}*/

#ifdef HAVE_ARPACK
	Math::LinAlg::Bind::Arpack::ArnoldiParameter< real_t > arnoldiparameters; // take the default Arnoldi paramters

	arnoldiparameters.eigenvalueRegion = Math::LinAlg::Bind::Arpack::SR;
	//SR : smallest real
	//LM : largest magnitude
	//SM : smallest magnitude
	arnoldiparameters.numberEigenvalues = 20;
	arnoldiparameters.numberDummyEigenvalues = 10;
	arnoldiparameters.tolerance = 0.000002; // to see the agreement of the methods you have to use a high precision.

	Math::LinAlg::Bind::Arpack::BasicArnoldi<std::complex<real_t>,std::complex<real_t> > barno(arnoldiparameters); // Basic arnoldi method

	DiracOperatorWrapper diracop(diracOperator);

	/*dirac_vector_t test1, test2, test3, test4;
	AlgebraUtils::generateRandomVector(test1);
	diracop.multiply(test1,test2);
	diracOperator->multiply(test3,test1);*/

	//std::cout << "Test: " << AlgebraUtils::differenceNorm(test2,test3) << std::endl;

	Math::DiracVect::ArnoldiDiracInterface< dirac_vector_t > arnoldioperator(&diracop,false);

	//arnoldioperator.multiply(reinterpret_cast<std::complex<real_t>*>(test1.getRaw()),reinterpret_cast<std::complex<real_t>*>(test4.getRaw()));

	//std::cout << "Test: " << AlgebraUtils::differenceNorm(test2,test4) << std::endl;

	Math::DiracVect::VectorEigenspace< dirac_vector_t > eigen1(arnoldiparameters.numberEigenvalues);

	barno(arnoldioperator, eigen1); // determine the smallest Eigenvalues of the Hermitian Dirac operator

	if (isOutputProcess()) std::cout << "Minimal eigenvalue: " << eigen1.getEigenvalue(0) << std::endl;

	if (isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("minimal_eigenvalues");

		for (int i = 0; i < 20; ++i) {
			output->write("minimal_eigenvalues", eigen1.getEigenvalue(i));
		}

		output->pop("minimal_eigenvalues");
	}

#endif

	delete diracOperator;

}

} /* namespace Update */
