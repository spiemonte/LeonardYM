/*
 * PolyakovLoop.cpp
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

#include "PolyakovLoop.h"
#include "io/GlobalOutput.h"

namespace Update {

PolyakovLoop::PolyakovLoop() { }

PolyakovLoop::~PolyakovLoop() { }

void PolyakovLoop::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t::Layout Layout;
	typedef extended_gauge_lattice_t LT;

	std::string write_3D_config = environment.configurations.get<std::string>("PolyakovLoop::write_polyakov_loop_config");

	if (environment.measurement && isOutputProcess()) {
                GlobalOutput* output = GlobalOutput::getInstance();
                output->push("polyakov");
        }

	for (int nu = 3; nu >= 0; --nu) {
		long_real_t polyakovLoopRe = 0;
		long_real_t polyakovLoopIm = 0;

		extended_gauge_lattice_t tmp = environment.gaugeLinkConfiguration;
		extended_gauge_lattice_t swap;
		extended_gauge_lattice_t polyakov;

#pragma omp parallel for
		for (int site = 0; site < Layout::localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				set_to_identity(polyakov[site][mu]);
			}
		}

		for (int t = 0; t < Layout::glob[nu]; ++t) {
#pragma omp parallel for
			for (int site = 0; site < Layout::localsize; ++site) {
				if (Layout::globalIndex(site,nu) == 0) {
					polyakov[site][nu] = polyakov[site][nu]*tmp[site][nu];
				}
			}

			//Antialias
			swap = tmp;

#pragma omp parallel for
			for (int site = 0; site < Layout::localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					tmp[site][mu] = swap[LT::sup(site,nu)][mu];
				}
			}
			tmp.updateHalo();
		}

#pragma omp parallel for reduction(+:polyakovLoopRe,polyakovLoopIm)
		for (int site = 0; site < Layout::localsize; ++site) {
			if (Layout::globalIndex(site,nu) == 0) {
				std::complex<real_t> polyakovLoop = trace(polyakov[site][nu]);
				polyakovLoopRe += real(polyakovLoop);
				polyakovLoopIm += imag(polyakovLoop);
			}
		}

		reduceAllSum(polyakovLoopRe);
		reduceAllSum(polyakovLoopIm);

		unsigned int spatialVolume = Layout::glob_spatial_volume;

		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();

			std::cout << "Polyakov Loop in the direction " << nu << " is (re) " << polyakovLoopRe/(numberColors*spatialVolume) << std::endl;
			std::cout << "Polyakov Loop in the direction " << nu << " is (im) " << polyakovLoopIm/(numberColors*spatialVolume) << std::endl;

			output->write("polyakov", polyakovLoopRe/(numberColors*spatialVolume));
			output->write("polyakov", polyakovLoopIm/(numberColors*spatialVolume));

		}

		if (write_3D_config == "true" && nu == 3) {
			std::string config_name = environment.configurations.get<std::string>("PolyakovLoop::output_configuration_name");
			write_polyakov_loop_config(polyakov, config_name, std::complex<real_t>(polyakovLoopRe/(numberColors*spatialVolume), polyakovLoopIm/(numberColors*spatialVolume)));
		}
	}

	if (environment.measurement && isOutputProcess()) {
                GlobalOutput* output = GlobalOutput::getInstance();
                output->pop("polyakov");
	}
}

void PolyakovLoop::write_polyakov_loop_config(const extended_gauge_lattice_t& polyakov, const std::string& output_name, const std::complex<real_t>& average_polyakov_loop) const {
	typedef extended_gauge_lattice_t::Layout Layout;
	typedef extended_gauge_lattice_t LT;

	std::string output_directory = environment.configurations.get<std::string>("output_directory_configurations");
	int offset = environment.configurations.get<unsigned int>("output_offset");


	FILE* fout(NULL);
	if (isOutputProcess()) {
		std::ostringstream filenamestream;
		filenamestream << output_directory << output_name << "-polyakov_config-" << environment.sweep+offset;
		if (isOutputProcess()) std::cout << "OutputSweep::Starting field write on file " << filenamestream.str() << std::endl;
		fout = fopen(filenamestream.str().c_str(), "w");

		if (!fout) {
			if (isOutputProcess()) std::cout << "OutputSweep::File not writable!" << std::endl;
			exit(19);
			return;
		}
	}

	XDR xout;

	if (isOutputProcess()) {
		xdrstdio_create(&xout, fout, XDR_ENCODE);
	}

	// Write fundamental polyakov loop
#ifdef ENABLE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	float *write_field = new float[2*Layout::glob_x*numberColors*numberColors];
#endif

	int wrt(0);

	for (int z = 0; z < Layout::glob_z; ++z) {
		for (int y = 0; y < Layout::glob_y; ++y) {
#ifndef ENABLE_MPI
			for (int x = 0; x < Layout::glob_x; ++x) {
				int globsite = Layout::getGlobalCoordinate(x,y,z,0);

				GaugeGroup tmp = environment.gaugeLinkConfiguration[globsite][t_dir];
				for (unsigned int c1 = 0; c1 < numberColors; ++c1) {
					for (unsigned int c2 = 0; c2 < numberColors; ++c2) {
						float tmp1 = static_cast<float>(real(tmp(c1,c2)));
						float tmp2 = static_cast<float>(imag(tmp(c1,c2)));
						wrt += xdr_float(&xout, &tmp1);
						wrt += xdr_float(&xout, &tmp2);
					}
				}
			}
#endif
#ifdef ENABLE_MPI
			int count = 0;
			int rank = Layout::rankTable(0,y,z,0);
			for (int x = 0; x < Layout::glob_x; ) {
				while (x < Layout::glob_x && rank == Layout::rankTable(x,y,z,t)) {
					int globsite = Layout::getGlobalCoordinate(x,y,z,t);
					int localsite = Layout::localIndex[globsite];
							
					if (localsite != -1 && localsite < Layout::localsize) {
						GaugeGroup tmp = environment.gaugeLinkConfiguration[localsite][t_dir];

						for (unsigned int c1 = 0; c1 < numberColors; ++c1) {
							for (unsigned int c2 = 0; c2 < numberColors; ++c2) {
								write_field[count++] = static_cast<float>(real(tmp(c1,c2)));
								write_field[count++] = static_cast<float>(imag(tmp(c1,c2)));
							}
						}
			
					}
					else {
						count += 2*numberColors*numberColors;
					}

					++x;
				}
						
				if (isOutputProcess()) {
					if (rank == Layout::this_processor) {
						for (int ii = 0; ii < count; ++ii) {
							wrt += xdr_float(&xout, &write_field[ii]);
						}
					}
					else {
						MPI_Recv(write_field, count, MPI_FLOAT, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						for (int ii = 0; ii < count; ++ii) {
							wrt += xdr_float(&xout, &write_field[ii]);
						}
					}
				}
				else {
					if (rank == Layout::this_processor) {
						MPI_Send(write_field, count, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
					}
				}

				if (x < Layout::glob_x) {
						rank = Layout::rankTable(x,y,z,0);
						count = 0;
					}
				}
			}
#endif
		}
	}

#ifdef ENABLE_MPI
	delete[] write_field;
#endif

	if (isOutputProcess()) {
		if (wrt != 2 * numberColors*numberColors * Layout::glob_spatial_volume){
			std::cout << "OutputSweep::write error for gauge fields on lattice "	<< 0 << ": " << (wrt / (2 * numberColors*numberColors)) << " : " << static_cast<size_t>(2 * numberColors*numberColors * Layout::glob_spatial_volume) << std::endl;
		}
		else {
			std::cout << "OutputSweep::number of gauge fields written: " << wrt / (2 * numberColors*numberColors) << " using XDR-format " << std::endl;
		}
	}

	float tmp1 = static_cast<float>(real(polyakov_loop_average));
	float tmp2 = static_cast<float>(imag(polyakov_loop_average));
	wrt += xdr_float(&xout, &tmp1);
	wrt += xdr_float(&xout, &tmp2);
}

void PolyakovLoop::registerParameters(po::options_description& desc) {
	desc.add_options()
		("PolyakovLoop::write_polyakov_loop_config", po::value<std::string>()->default_value("false"), "Write the whole 3D configuration of the Polyakov loop")
		("PolyakovLoop::output_configuration_name", po::value<std::string>()->default_value("polyakov_config_"), "Output name to write the 3D config");
}

} /* namespace Update */
