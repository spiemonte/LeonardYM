#include "MeanScalarField.h"
#include <iostream>
#include "io/GlobalOutput.h"

namespace Update {

	MeanScalarField::MeanScalarField() : LatticeSweep() { }

	MeanScalarField::~MeanScalarField() { }

	void MeanScalarField::execute(environment_t& environment) {
		unsigned int aNf = environment.configurations.get<unsigned int>("adjoint_nf_scalars");

		if (aNf > 0) {
			long_real_t mean_field_re = 0., mean_field_squared = 0., mean_field_abs_re = 0.;

			for (unsigned int f = 0; f < aNf; ++f) {
#pragma omp parallel for reduction(+:mean_field_re, mean_field_squared, mean_field_abs_re)
				for (int site = 0; site < environment.adjoint_scalar_fields[f].localsize; ++site) {
					for (int c = 0; c < numberColors*numberColors - 1; ++c) {
						mean_field_squared += (environment.adjoint_scalar_fields[f][site][c]*environment.adjoint_scalar_fields[f][site][c]);

						mean_field_re += environment.adjoint_scalar_fields[f][site][c];

						mean_field_abs_re += fabs(environment.adjoint_scalar_fields[f][site][c]);
					}
				}
			}

			reduceAllSum(mean_field_squared);
			reduceAllSum(mean_field_re);
			reduceAllSum(mean_field_abs_re);

			if (environment.measurement && isOutputProcess()) {
				GlobalOutput* output = GlobalOutput::getInstance();
				output->push("adjoint_scalar_field");

				typedef extended_gauge_lattice_t::Layout Layout;
				std::cout << "MeanScalarField::Mean adjoint scalar field expectation value: " << mean_field_re/(aNf*(numberColors*numberColors - 1)*Layout::globalVolume) << std::endl;
				std::cout << "MeanScalarField::Squared mean adjoint scalar field expectation value: " << mean_field_squared/(aNf*(numberColors*numberColors - 1)*Layout::globalVolume) << std::endl;
				output->write("adjoint_scalar_field", mean_field_re/(aNf*(numberColors*numberColors - 1)*Layout::globalVolume));
				output->write("adjoint_scalar_field", mean_field_squared/(aNf*(numberColors*numberColors - 1)*Layout::globalVolume));
				output->write("adjoint_scalar_field", mean_field_abs_re/(aNf*(numberColors*numberColors - 1)*Layout::globalVolume));

				output->pop("adjoint_scalar_field");
			}
		}

		unsigned int fNf = environment.configurations.get<unsigned int>("fundamental_nf_scalars");

		if (fNf > 0) {
			long_real_t mean_field_re = 0., mean_field_squared = 0., mean_field_abs_re = 0.;

			for (unsigned int f = 0; f < fNf; ++f) {
#pragma omp parallel for reduction(+:mean_field_re, mean_field_squared, mean_field_abs_re)
				for (int site = 0; site < environment.fundamental_scalar_fields[f].localsize; ++site) {
					for (int c = 0; c < numberColors; ++c) {
						mean_field_squared += real(conj(environment.fundamental_scalar_fields[f][site][c])*environment.fundamental_scalar_fields[f][site][c]);

						mean_field_re += real(environment.fundamental_scalar_fields[f][site][c]);

						mean_field_abs_re += fabs(environment.fundamental_scalar_fields[f][site][c]);
					}
				}
			}

			if (environment.measurement && isOutputProcess()) {
				GlobalOutput* output = GlobalOutput::getInstance();
				output->push("fundamental_scalar_field");

				typedef extended_gauge_lattice_t::Layout Layout;
				std::cout << "MeanScalarField::Mean fundamental scalar field expectation value: " << mean_field_re/(fNf*numberColors*Layout::globalVolume) << std::endl;
				std::cout << "MeanScalarField::Squared mean fundamental scalar field expectation value: " << mean_field_squared/(fNf*numberColors*Layout::globalVolume) << std::endl;
				output->write("fundamental_scalar_field", mean_field_re/(fNf*numberColors*Layout::globalVolume));
				output->write("fundamental_scalar_field", mean_field_squared/(fNf*numberColors*Layout::globalVolume));
				output->write("fundamental_scalar_field", mean_field_abs_re/(fNf*numberColors*Layout::globalVolume));

				output->pop("fundamental_scalar_field");
			}
		}

	}

	long_real_t MeanScalarField::meanValueSquared(const extended_adjoint_real_color_vector_t& field) const {
		long_real_t mean_field_squared = 0.;

#pragma omp parallel for reduction(+:mean_field_squared)
		for (int site = 0; site < field.localsize; ++site) {
			for (int c = 0; c < numberColors*numberColors - 1; ++c) {
				mean_field_squared += (field[site][c]*field[site][c]);
			}
		}

		reduceAllSum(mean_field_squared);

		return mean_field_squared;
	}

} /* namespace Update */
