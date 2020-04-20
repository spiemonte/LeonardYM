#include "OutputSweep.h"
#include "Environment.h"
#include "wilson_loops/Plaquette.h"
#include "utils/ToString.h"
#include <fstream>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#include <sys/time.h>
#include <iomanip>

namespace Update {

OutputSweep::OutputSweep() : LatticeSweep() { }

OutputSweep::~OutputSweep() { }

void OutputSweep::execute(environment_t& environment) {

	timeval start, stop, result;
	gettimeofday(&start,NULL);	
	
	std::string format_name;
	try {
		format_name = environment.configurations.get<std::string>("output_format_name");
	}
	catch (NotFoundOption& ex) {
		format_name = environment.configurations.get<std::string>("format_name");
	}
	
	typedef extended_gauge_lattice_t::Layout LT;
	
	if (format_name == "leonard_format") {
		std::string output_name = environment.configurations.get<std::string>("output_configuration_name");
		std::string output_directory = environment.configurations.get<std::string>("output_directory_configurations");
		int offset = environment.configurations.get<unsigned int>("output_offset");

		if (isOutputProcess()) std::cout << "OutputSweep::Writing configuration to file " << output_directory << output_name << "_" << toString(environment.sweep+offset) << std::endl;		

		FILE* fout(NULL);
		std::string output_file = output_directory+output_name+"_"+toString(environment.sweep+offset)+"_"+toString(LT::this_processor)+".txt";

		fout = fopen(output_file.c_str(), "w");

		if (!fout) {
			std::cout << "File not writeble!" << std::endl;
			return;
		}

		XDR xout;
		xdrstdio_create(&xout, fout, XDR_ENCODE);

		for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
			xdr_int(&xout, &LT::globalCoordinate[site].x);
			xdr_int(&xout, &LT::globalCoordinate[site].y);
			xdr_int(&xout, &LT::globalCoordinate[site].z);
			xdr_int(&xout, &LT::globalCoordinate[site].t);
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (int i = 0; i < numberColors; ++i) {
					for (int j = 0; j < numberColors; ++j) {
						double re = environment.gaugeLinkConfiguration[site][mu].at(i,j).real();
						double im = environment.gaugeLinkConfiguration[site][mu].at(i,j).imag();
						xdr_double(&xout,&re);
						xdr_double(&xout,&im);
					}
				}
			}
		}

		xdr_destroy( &xout);
		fclose(fout);

		double plaquette = Plaquette::temporalPlaquette(environment.gaugeLinkConfiguration);

		if (isOutputProcess()) {
			std::string descriptor_name = output_directory+output_name+"_"+toString(environment.sweep+offset)+".descriptor.txt";
			std::fstream descriptor;
			descriptor.open(descriptor_name.c_str(), std::fstream::out);
			descriptor << LT::glob_x << " " << LT::glob_y << " " << LT::glob_z << " " << LT::glob_t << std::endl;
			descriptor << LT::pgrid_x << " " << LT::pgrid_y << " " << LT::pgrid_z << " " << LT::pgrid_t << std::endl;
			descriptor << LT::numberProcessors << std::endl;
			descriptor << LT::localsize << std::endl;
			descriptor << std::setprecision(13) << plaquette << std::endl;
			descriptor.close();
		}

	}
	else if (format_name == "muenster_format") {
#if NUMCOLORS == 2
		typedef extended_gauge_lattice_t::Layout Layout;
		std::string output_name = environment.configurations.get<std::string>("output_configuration_name");
		std::string output_directory = environment.configurations.get<std::string>("output_directory_configurations");
		int offset = environment.configurations.get<unsigned int>("output_offset");


		FILE* fout(NULL);
		if (isOutputProcess()) {
			std::ostringstream filenamestream;
			filenamestream << output_directory << 0 << ":" << environment.sweep+offset << ":" << output_name;
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

		// Write fundamental gauge field
#ifdef ENABLE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		float *write_field = new float[Layout::glob_x*16]; // 16=4*2*2
#endif
		int wrt(0);

		for (int t = 0; t < Layout::glob_t; ++t) {
			for (int z = 0; z < Layout::glob_z; ++z) {
				for (int y = 0; y < Layout::glob_y; ++y) {
#ifndef ENABLE_MPI
					for (int x = 0; x < Layout::glob_x; ++x) {
						int globsite = Layout::getGlobalCoordinate(x,y,z,t);

						for (unsigned int mu = 0; mu < 4; ++mu) {
							// Assuming now Istvans format with U^dag saved insted of U
							GaugeGroup tmpl = environment.gaugeLinkConfiguration[globsite][mu];
							GaugeGroup tmp = htrans(tmpl);
							float tmp1 = static_cast<float>(real(tmp(0,0)));
							float tmp2 = static_cast<float>(imag(tmp(0,0)));
							wrt += xdr_float(&xout, &tmp1);
							wrt += xdr_float(&xout, &tmp2);
							tmp1 = static_cast<float>(real(tmp(0,1)));
							tmp2 = static_cast<float>(imag(tmp(0,1)));
							wrt += xdr_float(&xout, &tmp1);
							wrt += xdr_float(&xout, &tmp2);
						}
					}
#endif
#ifdef ENABLE_MPI
					int count = 0;
					int rank = Layout::rankTable(0,y,z,t);
					for (int x = 0; x < Layout::glob_x; ) {
						while (x < Layout::glob_x && rank == Layout::rankTable(x,y,z,t)) {
							int globsite = Layout::getGlobalCoordinate(x,y,z,t);
							int localsite = Layout::localIndex[globsite];
							
							if (localsite != -1 && localsite < Layout::localsize) {
								for (unsigned int mu = 0; mu < 4; ++mu) {
									//Assuming now Istvans format with U^dag saved insted of U
									GaugeGroup tmpl = environment.gaugeLinkConfiguration[localsite][mu];
									GaugeGroup tmp = htrans(tmpl);
									write_field[count++] = static_cast<float>(real(tmp(0,0)));
									write_field[count++] = static_cast<float>(imag(tmp(0,0)));
									write_field[count++] = static_cast<float>(real(tmp(0,1)));
									write_field[count++] = static_cast<float>(imag(tmp(0,1)));
								}
							}
							else {
								count += 16;
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
							rank = Layout::rankTable(x,y,z,t);
							count = 0;
						}
					}
#endif
				}
			}
		}

#ifdef ENABLE_MPI
		delete[] write_field;
#endif

		if (isOutputProcess()) {
			if (wrt != 16 * Layout::globalVolume){
				std::cout << "OutputSweep::write error for gauge fields on lattice "	<< 0 << ": " << (wrt / 4) << " : " << static_cast<size_t>(4 * Layout::globalVolume) << std::endl;
			}
			else {
				std::cout << "OutputSweep::number of gauge fields written: " << wrt / 4 << " using XDR-format " << std::endl;
			}
		}

		double tplaq = Plaquette::temporalPlaquette(environment.gaugeLinkConfiguration);
		
		if (isOutputProcess()) {
			// Write Plaquette value, correction factor and smallest eigenvalue
			float tmp = static_cast<float >(tplaq);
			xdr_float(&xout, &tmp);
			tmp = static_cast<float >(0.);
			xdr_float(&xout, &tmp);
			tmp = static_cast<float >(0.);
			xdr_float(&xout, &tmp);
			//Write parameters
			tmp = static_cast<float >(environment.configurations.get<double>("beta"));
			xdr_float(&xout, &tmp);

			std::vector<real_t> kappas; //< in the update several kappas are considered
			std::vector<real_t> alphas; //< information about the spectrum and the polynoms (update)
			std::vector<real_t> epsilons; //< information about the spectrum and the polynoms (update)
			std::vector<real_t> lambdas; //< information about the spectrum and the polynoms (update)
			std::vector<int> polyorder;

			kappas.push_back(environment.configurations.get<double>("kappa"));//TODO
			alphas.push_back(0);//TODO
			epsilons.push_back(0);//TODO
			lambdas.push_back(0);//TODO

			polyorder.push_back(0);//TODO
			polyorder.push_back(0);//TODO
			polyorder.push_back(0);//TODO
			polyorder.push_back(0);//TODO
			polyorder.push_back(0);//TODO

			for (size_t ii(0); ii < kappas.size() - 1; ++ii) { // This way
				tmp = -static_cast<float>(kappas[ii]);
				xdr_float(&xout, &tmp);
			}

			tmp = static_cast<float>(kappas[kappas.size() - 1]); // ensures

			xdr_float(&xout, &tmp);
			// compatib.

			for (size_t ii(0); ii < alphas.size(); ++ii) {
				tmp = static_cast<float >(alphas[ii]);
				xdr_float(&xout, &tmp);
			}

			for (size_t ii(0); ii < alphas.size(); ++ii) {
				tmp = static_cast<float >(epsilons[ii]);
				xdr_float(&xout, &tmp);
				tmp = static_cast<float >(lambdas[ii]);
				xdr_float(&xout, &tmp);
			}

			{
				int tmpint(offset + environment.sweep);
				xdr_int(&xout, &tmpint);
			}

			for (size_t kk(0); kk < polyorder.size(); ++kk) {
				int ord(polyorder[kk]);
				xdr_int(&xout, &ord);
			}

			int lats[4] = { Layout::glob_x, Layout::glob_y, Layout::glob_z, Layout::glob_t };

			for (int ii(0); ii < 4; ++ii) xdr_int(&xout, &lats[ii]);

			int lvrc = 0;//TODO

			xdr_int(&xout, &lvrc);
		}

		if (isOutputProcess()) {
			xdr_destroy(&xout);
			fclose(fout);
		}
#endif
#if NUMCOLORS == 3
		typedef extended_gauge_lattice_t::Layout Layout;
		std::string output_name = environment.configurations.get<std::string>("output_configuration_name");
		std::string output_directory = environment.configurations.get<std::string>("output_directory_configurations");
		int offset = environment.configurations.get<unsigned int>("output_offset");


		FILE* fout(NULL);
		if (isOutputProcess()) {
			std::ostringstream filenamestream;
			filenamestream << output_directory << output_name << ".config-"<< environment.sweep+offset;
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

		// Write fundamental gauge field
#ifdef ENABLE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		float *write_field = new float[Layout::glob_x*48]; // 48=4*6*2
#endif
		int wrt(0);

		for (int t = 0; t < Layout::glob_t; ++t) {
			for (int z = 0; z < Layout::glob_z; ++z) {
				for (int y = 0; y < Layout::glob_y; ++y) {
#ifndef ENABLE_MPI
					for (int x = 0; x < Layout::glob_x; ++x) {
						int globsite = Layout::getGlobalCoordinate(x,y,z,t);

						
						for (unsigned int mu = 0; mu < 4; ++mu) {
							GaugeGroup mt = environment.gaugeLinkConfiguration[globsite][mu];
							for (size_t ii = 0; ii < 3; ++ii) {
								for (size_t jj = 0; jj < 2; ++jj) {
									/** Assuming now Istvans format with U^* saved insted of U*/
									float tmp = mt.at(ii,jj).real(), tmp2 = -mt.at(ii,jj).imag();
									wrt += xdr_float(&xout, &tmp);
									wrt += xdr_float(&xout, &tmp2);
								}
							}
						}
					}
#endif
#ifdef ENABLE_MPI
					int count = 0;
					int rank = Layout::rankTable(0,y,z,t);
					for (int x = 0; x < Layout::glob_x; ) {
						while (x < Layout::glob_x && rank == Layout::rankTable(x,y,z,t)) {
							int globsite = Layout::getGlobalCoordinate(x,y,z,t);
							int localsite = Layout::localIndex[globsite];
							
							if (localsite != -1 && localsite < Layout::localsize) {
								for (unsigned int mu = 0; mu < 4; ++mu) {
									GaugeGroup mt = environment.gaugeLinkConfiguration[localsite][mu];
									for (size_t ii = 0; ii < 3; ++ii) {
										for (size_t jj = 0; jj < 2; ++jj) {
											/** Assuming now Istvans format with U^* saved insted of U*/
											float tmp = mt.at(ii,jj).real(), tmp2 = -mt.at(ii,jj).imag();
											write_field[count++] = tmp;
											write_field[count++] = tmp2;
										}
									}
								}
							}
							else {
								count += 48;
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
							rank = Layout::rankTable(x,y,z,t);
							count = 0;
						}
					}
#endif
				}
			}
		}

#ifdef ENABLE_MPI
		delete[] write_field;
#endif

		double t_plaq = Plaquette::temporalPlaquette(environment.gaugeLinkConfiguration);

		if (isOutputProcess()) {
			if (wrt != 48 * Layout::globalVolume){
				std::cout << "OutputSweep::write error for gauge fields on lattice "	<< 0 << ": " << (wrt / 4) << " : " << static_cast<size_t>(4 * Layout::globalVolume) << std::endl;
			}
			else {
				std::cout << "OutputSweep::number of gauge fields written: " << wrt / 4 << " using XDR-format " << std::endl;
			}
		}
		
		if (isOutputProcess()) {
			if (wrt != 48 * Layout::globalVolume) {
				std::cout << " write error for gauge fields on lattice " << ": " << wrt << std::endl;
			}
			else {
				std::cout << "OutputSweep::Finished writing " << wrt / 4 << " gauge fields" << std::endl;
			}

			/// - read the global lattice extend and lvrco
          		int lats[4] = { Layout::glob_x, Layout::glob_y, Layout::glob_z, Layout::glob_t };

			for (int ii(0); ii < 4; ++ii) xdr_int(&xout, &lats[ii]);

          		double beto = static_cast<float >(environment.configurations.get<double>("beta"));
          		xdr_double(&xout, &beto);

     			double c1(-1./12.);
          		xdr_double(&xout, &c1);

			double fplq = t_plaq;	
			xdr_double(&xout, &fplq);
		}

		if (isOutputProcess()) {
			xdr_destroy(&xout);
			fclose(fout);
		}
		
#endif
#if NUMCOLORS > 3
		if (isOutputProcess()) {
			std::cout << "OuputSweep::munster_format not implemented for NC>3!" << std::endl;
		}
#endif
	}
	else if (format_name == "mathematica_format") {
		std::string output_name = environment.configurations.get<std::string>("output_configuration_name");
		std::string output_directory = environment.configurations.get<std::string>("output_directory_configurations");
		int offset = environment.configurations.get<unsigned int>("output_offset");

		

		//FILE* fout(NULL);
		std::string output_file = output_directory+output_name+"_"+toString(environment.sweep+offset)+"_"+toString(LT::this_processor)+".txt";
		std::fstream fout(output_file.c_str(),std::fstream::out);

		//fout.open(output_file.c_str(), "w");

		/*if (!fout) {
			std::cout << "File not writeble!" << std::endl;
			return;
		}*/
		
		fout << "{";
		for (int site = 0; site < environment.gaugeLinkConfiguration.localsize-1; ++site) {
			fout << "{{" << LT::globalCoordinate[site].x << ",";
			fout << LT::globalCoordinate[site].y << ",";
			fout << LT::globalCoordinate[site].z << ",";
			fout << LT::globalCoordinate[site].t << "},";
			fout << toString(environment.gaugeLinkConfiguration[site][0]) << ",";
			fout << toString(environment.gaugeLinkConfiguration[site][1]) << ",";
			fout << toString(environment.gaugeLinkConfiguration[site][2]) << ",";
			fout << toString(environment.gaugeLinkConfiguration[site][3]) << "},";
		}
		{
			int site = environment.gaugeLinkConfiguration.localsize-1;
			fout << "{{" << LT::globalCoordinate[site].x << ",";
			fout << LT::globalCoordinate[site].y << ",";
			fout << LT::globalCoordinate[site].z << ",";
			fout << LT::globalCoordinate[site].t << "},";
			fout << toString(environment.gaugeLinkConfiguration[site][0]) << ",";
			fout << toString(environment.gaugeLinkConfiguration[site][1]) << ",";
			fout << toString(environment.gaugeLinkConfiguration[site][2]) << ",";
			fout << toString(environment.gaugeLinkConfiguration[site][3]) << "}";
		}
		fout << "}";
		fout.close();
	}
	
	gettimeofday(&stop,NULL);
	timersub(&stop,&start,&result);
	if (isOutputProcess()) std::cout << "OutputSweep::Configuration written in: " << (double)result.tv_sec + result.tv_usec/1000000.0 << " sec" << std::endl;
	
	/*
	
	

	int blocksize = 4*sizeof(int) + 4*numberColors*numberColors*2*8;

	MPI_File fh;

	

	if (isOutputProcess()) std::cout << output_file << std::endl;

	char s[4096];
	for (unsigned int i = 0; i < 4096 && i < output_file.size(); ++i) {
		s[i] = output_file[i];
		s[i + 1] = 0;
	}
	
	MPI_File_open(MPI_COMM_WORLD, s, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

	MPI_File_set_size(fh, LT::globalVolume*blocksize);
	
	for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
		int coord = LT::getGlobalCoordinate(LT::globalCoordinate[site]);
		MPI_File_seek(fh, coord*blocksize, MPI_SEEK_SET);
		int buf[4] = {LT::globalCoordinate[site].x, LT::globalCoordinate[site].y, LT::globalCoordinate[site].z, LT::globalCoordinate[site].t};
		MPI_File_write(fh, buf, 4, MPI_INT, MPI_STATUS_IGNORE);
		MPI_File_write(fh, &environment.gaugeLinkConfiguration[site], 4*numberColors*numberColors*2, MPI_DOUBLE, MPI_STATUS_IGNORE);
	}
	
	MPI_File_close(&fh);*/
}


} /* namespace Update */
