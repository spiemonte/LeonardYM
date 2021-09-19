#include "ReadStartGaugeConfiguration.h"
#include "wilson_loops/Plaquette.h"
#include <string>
#include <fstream>
#ifndef ENABLE_MPI
#ifdef __APPLE__
typedef bool bool_t;
typedef int enum_t;
typedef long long quad_t;
typedef unsigned long long u_quad_t;
typedef int caddr_t;
#define TRUE (true)
#define FALSE (false)
#endif
#endif
#include <rpc/xdr.h>
#include "utils/ToString.h"

namespace Update {

real_t absq(const std::complex<real_t>& num) { return real(num)*real(num) + imag(num)*imag(num); }

ReadStartGaugeConfiguration::ReadStartGaugeConfiguration() : StartGaugeConfiguration() { }

ReadStartGaugeConfiguration::~ReadStartGaugeConfiguration() { }

void ReadStartGaugeConfiguration::execute(environment_t& environment) {
	int numberfile = environment.configurations.get<unsigned int>("input_number");
	bool success = this->readConfiguration(environment,numberfile);
	if (!success) {
		if (isOutputProcess()) std::cout << "ReadStartGaugeConfiguration::Reading failed!" << std::endl;
		exit(49);
	}
}

bool ReadStartGaugeConfiguration::readConfiguration(environment_t& environment, int numberfile) {
	std::string format_name;
	try {
		format_name = environment.configurations.get<std::string>("input_format_name");
	}
	catch (NotFoundOption& ex) {
		format_name = environment.configurations.get<std::string>("format_name");
	}
	
	typedef extended_gauge_lattice_t::Layout LT;
	typedef extended_gauge_lattice_t::Layout Layout;
	double read_plaquette;

	bool success = true;

	if (format_name == "leonard_format") {
		std::string directory = environment.configurations.get<std::string>("input_directory_configurations");
		std::string input_name = environment.configurations.get<std::string>("input_name");
		
		if (isOutputProcess()) std::cout << "ReadStartGaugeConfiguration::Reading configuration from file " << directory << input_name << "_" << toString(numberfile) << std::endl;

		int read_glob_x, read_glob_y, read_glob_z, read_glob_t;
		int read_pgrid_x, read_pgrid_y, read_pgrid_z, read_pgrid_t;
		int numberProcessors;
		int read_localsize;

		std::string descriptor_name = directory+input_name+"_"+toString(numberfile)+".descriptor.txt";
		std::fstream descriptor;
		if (isOutputProcess())  std::cout << "ReadStartGaugeConfiguration::Reading configuration descriptor " << descriptor_name << std::endl;
		descriptor.open(descriptor_name.c_str(), std::fstream::in);
		descriptor >> read_glob_x >> read_glob_y >> read_glob_z >> read_glob_t;
		descriptor >> read_pgrid_x >> read_pgrid_y >> read_pgrid_z >> read_pgrid_t;
		descriptor >> numberProcessors;
		descriptor >> read_localsize;
		descriptor >> read_plaquette;
		descriptor.close();

		if (read_glob_x != LT::glob_x || read_glob_y != LT::glob_y || read_glob_z != LT::glob_z || read_glob_t != LT::glob_t) {
			std::cout << "Different lattice size in reading configuration!" << std::endl;
			std::cout << "Configured: " << LT::glob_x << " " << LT::glob_y << " " << LT::glob_z << " " << LT::glob_t << std::endl;
			std::cout << "Readed: " << read_glob_x << " " << read_glob_y << " " << read_glob_z << " " << read_glob_t << std::endl;
			return false;
		}

		for (int processor = 0; processor < numberProcessors; ++processor) {
			std::string input_file = directory+input_name+"_"+toString(numberfile) + "_" + toString(processor) + ".txt";

			FILE* fin(NULL);
			fin = fopen(input_file.c_str(), "r");

			if (!fin) {
				std::cout << "ReadStartGaugeConfiguration::File not readble!" << std::endl;
				return false;
			}

			XDR xin;
			xdrstdio_create(&xin, fin, XDR_DECODE);

			for (int index = 0; index < read_localsize; ++index) {
				int x, y, z, t;
				xdr_int(&xin, &x);
				xdr_int(&xin, &y);
				xdr_int(&xin, &z);
				xdr_int(&xin, &t);
				int site = LT::localIndex[LT::getGlobalCoordinate(x,y,z,t)];
				if (site != -1) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						for (int i = 0; i < numberColors; ++i) {
							for (int j = 0; j < numberColors; ++j) {
								double re;
								double im;
								xdr_double(&xin,&re);
								xdr_double(&xin,&im);
								environment.gaugeLinkConfiguration[site][mu].at(i,j) = std::complex<real_t>(re,im);
							}
						}
					}
				}
				else {
					//TODO
					for (unsigned int mu = 0; mu < 4; ++mu) {
						for (int i = 0; i < numberColors; ++i) {
							for (int j = 0; j < numberColors; ++j) {
								double tmp;
								xdr_double(&xin,&tmp);
								xdr_double(&xin,&tmp);
							}
						}
					}
				}
			}

			xdr_destroy(&xin);
			fclose(fin);
		}
	}
	else if (format_name == "muenster_format") {
#if NUMCOLORS == 2
		std::string directory = environment.configurations.get<std::string>("input_directory_configurations");
		std::string input_name = environment.configurations.get<std::string>("input_name");

		std::string filename = directory + "0:" + toString(numberfile) + ":" + input_name;

		FILE* fin(NULL);


		fin = fopen(filename.c_str(), "r");

		if (!fin) {
			if (isOutputProcess()) std::cout << "ReadStartGaugeConfiguration::File " << filename << " impossible to read!" << std::endl;
			return false;
		}


		//initialisation of xdr
		XDR xin;
		xdrstdio_create(&xin, fin, XDR_DECODE);

		//read fundamental gauge field
		int red = 0;

		for (int t = 0; t < Layout::glob_t; ++t) {
			for (int z = 0; z < Layout::glob_z; ++z) {
				for (int y = 0; y < Layout::glob_y; ++y) {
					for (int x = 0; x < Layout::glob_x; ++x) {
						int globsite = Layout::getGlobalCoordinate(x,y,z,t);
#ifndef ENABLE_MPI
						for (size_t mu = 0; mu < 4; ++mu) {
							for (size_t ii = 0; ii < 2; ++ii) {
								float tmp, tmp2;
								red += xdr_float(&xin, &tmp);
								red += xdr_float(&xin, &tmp2);
								/** Assuming now Istvans format with U^dag saved insted of U*/
								environment.gaugeLinkConfiguration[globsite][mu].at(0,ii) = (ii == 0 ? conj(complex(tmp, tmp2)) : -complex(tmp, tmp2) );
								environment.gaugeLinkConfiguration[globsite][mu].at(1,ii==0 ? 1:0) = (ii == 0 ? complex(tmp, tmp2) : conj(complex(tmp, tmp2)) );
							}
						}
#endif
#ifdef ENABLE_MPI
						int localsite = Layout::localIndex[globsite];
						if (localsite != -1) {
							for (size_t mu = 0; mu < 4; ++mu) {
								for (size_t ii = 0; ii < 2; ++ii) {
									float tmp, tmp2;
									red += xdr_float(&xin, &tmp);
									red += xdr_float(&xin, &tmp2);
									/** Assuming now Istvans format with U^dag saved insted of U*/
									environment.gaugeLinkConfiguration[localsite][mu].at(0,ii) = (ii == 0 ? conj(complex(tmp, tmp2)) : -complex(tmp, tmp2) );
									environment.gaugeLinkConfiguration[localsite][mu].at(1,ii==0 ? 1:0) = (ii == 0 ? complex(tmp, tmp2) : conj(complex(tmp, tmp2)) );
								}
							}
						}
						else {
							for (size_t mu = 0; mu < 4; ++mu) {
								for (size_t ii = 0; ii < 2; ++ii) {
									float tmp, tmp2;
									red += xdr_float(&xin, &tmp);
									red += xdr_float(&xin, &tmp2);
								}
							}
						}
#endif
					}
				}
			}
		}




		// - check if the number of read fields is correct - reading failure if not
		if (isOutputProcess()) {
			if (red != 16 * Layout::globalVolume) {
				std::cout << " read error for gauge fields on lattice " << ": " << red << std::endl;
				success = false;
			}
			else {
				std::cout << "ReadStartGaugeConfiguration::Finished reading " << red / 4 << " gauge fields" << std::endl;
				std::cout << "ReadStartGaugeConfiguration::From file " << filename << std::endl;
			}

			// - read Plaquette value, correction factor and smallest eigenvalue, beta
			float plaq(0.0);
			xdr_float(&xin, &plaq);
			read_plaquette = plaq;

			float cf(0.0);
			xdr_float(&xin, &cf);

			float evs(0.0);
			xdr_float(&xin, &evs);

			float beto(0.0);
			xdr_float(&xin, &beto);

			// - read several kappas, alphas and epsilons the number of these paramters (fermiono) is the number of input kappas with a negative sign+1. (Of course only the fabs of these is used in the update.)
			// Note that in the measurements only the last kappa value (with positive sign) counts.
			float kappo(-1.0);

			int fermiono(0);
			// stop reading if positive sign appears.

			while (kappo < 0.0) {
				xdr_float(&xin, &kappo);
				fermiono++;
			}

			float alpho(0.0);

			for (int num1(0); num1 < fermiono; ++num1)	xdr_float(&xin, &alpho);

			float epsilo(0.0), lambdo(0.0);

			for (int num1(0); num1 < fermiono; ++num1) {
				xdr_float(&xin, &epsilo);
				xdr_float(&xin, &lambdo);
			}

			// - read sweep number (relevant for the output) and polynomial order (not used here)
			int sweepo(0);

			xdr_int(&xin, &sweepo);

			for (int num2(0); num2 < fermiono; ++num2)
				for (int num1(0); num1 < 5; ++num1) {
					int polo;
					xdr_int(&xin, &polo);
				}

			// - read the global lattice extend
			int lato[4];

			for (int num1 = 0; num1 < 4; ++num1) xdr_int(&xin, &lato[num1]);

			int lvrco(0);
			xdr_int(&xin, &lvrco);

			// - compare read parameters to the ones in measureconfig.dat / Stop program if they are not the same (This is first done directly on the master node without communication of all the paramters)
			// - test for lattice size - new included: reading failure if given incorrect
			if (static_cast<int>(lato[0]) != Layout::glob_x
					|| static_cast<int>(lato[1]) != Layout::glob_y
					|| static_cast<int>(lato[2]) != Layout::glob_z
					|| static_cast<int>(lato[3]) != Layout::glob_t) {
				if (isOutputProcess()) std::cout << "Lattice size is different:\n" << lato[0] << " " << lato[1] << " "
						<< lato[2] << " " << lato[3] << "\n" << Layout::glob_x << " "
						<< Layout::glob_y << " " << Layout::glob_z << " "
						<< Layout::glob_t;
				success = false;
			}

			//Load the beta and kappa
			double beta = environment.configurations.get<double>("beta");
			double kappa = environment.configurations.get<double>("kappa");

			// - test for Beta and Kappa - new included: reading failure if given incorrect
			if (std::fabs(beta - beto) / (std::fabs(beta) + std::fabs(beto)) > 0.0000001
					|| std::fabs(kappa - kappo) / (std::fabs(kappa) + std::fabs(kappo)) > 0.0000001) {
				if (isOutputProcess()) std::cout << "Some parameter values are changed on lattice "
						<< ":\n" << " Beta, Kappa:\n" << beto << " " << kappo << "\n"
						<< beta << " " << kappa;
			}

		} // on of Io node only

		xdr_destroy( &xin);
		fclose(fin);
#endif
#if NUMCOLORS == 3
		std::string directory = environment.configurations.get<std::string>("input_directory_configurations");
		std::string input_name = environment.configurations.get<std::string>("input_name");

		std::string filename = directory + input_name + ".config-" + toString(numberfile);

		FILE* fin(NULL);


		fin = fopen(filename.c_str(), "r");

		if (!fin) {
			if (isOutputProcess()) std::cout << "ReadStartGaugeConfiguration::File " << filename << " impossible to read!" << std::endl;
			return false;
		}


		//initialisation of xdr
		XDR xin;
		xdrstdio_create(&xin, fin, XDR_DECODE);

		//read fundamental gauge field
		int red = 0;

		for (int t = 0; t < Layout::glob_t; ++t) {
			for (int z = 0; z < Layout::glob_z; ++z) {
				for (int y = 0; y < Layout::glob_y; ++y) {
					for (int x = 0; x < Layout::glob_x; ++x) {
						int globsite = Layout::getGlobalCoordinate(x,y,z,t);
#ifndef ENABLE_MPI
						GaugeGroup mt;
						for (size_t mu = 0; mu < 4; ++mu) {
							for (size_t ii = 0; ii < 3; ++ii) {
								for (size_t jj = 0; jj < 2; ++jj) {
									float tmp, tmp2;
									red += xdr_float(&xin, &tmp);
									red += xdr_float(&xin, &tmp2);
									/** Assuming now Istvans format with U^* saved insted of U*/
									mt.at(ii,jj) = std::complex<real_t>(tmp, -tmp2);
								}
							}
							
							{
         							const std::complex<real_t> temp = 1.0 / sqrt(absq(mt(0,0)) + absq(mt(1,0)) + absq(mt(2,0)));
								mt(0,0) *= temp;
								mt(1,0) *= temp;
								mt(2,0) *= temp;
							}

							{
								const std::complex<real_t> temp = mt(0,1) * conj(mt(0,0)) + mt(1,1) * conj(mt(1,0)) + mt(2,1) * conj(mt(2,0));
								mt(0,1) -= temp * mt(0,0);
								mt(1,1) -= temp * mt(1,0);
								mt(2,1) -= temp * mt(2,0);
							}

							{
								const std::complex<real_t> temp = 1.0 / sqrt(absq(mt(0,1)) + absq(mt(1,1)) + absq(mt(2,1)));
								mt(0,1) *= temp;
								mt(1,1) *= temp;
 								mt(2,1) *= temp;
							}

							mt(0,2) = conj(mt(1,0)) * conj(mt(2,1)) - conj(mt(2,0)) * conj(mt(1,1));
							mt(1,2) = conj(mt(2,0)) * conj(mt(0,1)) - conj(mt(0,0)) * conj(mt(2,1));
							mt(2,2) = conj(mt(0,0)) * conj(mt(1,1)) - conj(mt(1,0)) * conj(mt(0,1));

							environment.gaugeLinkConfiguration[globsite][mu] = mt;
						}
#endif
#ifdef ENABLE_MPI
						int localsite = Layout::localIndex[globsite];
						if (localsite != -1) {
							GaugeGroup mt;
							for (size_t mu = 0; mu < 4; ++mu) {
								for (size_t ii = 0; ii < 3; ++ii) {
									for (size_t jj = 0; jj < 2; ++jj) {
										float tmp, tmp2;
										red += xdr_float(&xin, &tmp);
										red += xdr_float(&xin, &tmp2);
										/** Assuming now Istvans format with U^* saved insted of U*/
										mt.at(ii,jj) = std::complex<real_t>(tmp, -tmp2);
									}
								}
							
								{
         								const std::complex<real_t> temp = 1.0 / sqrt(absq(mt(0,0)) + absq(mt(1,0)) + absq(mt(2,0)));
									mt(0,0) *= temp;
									mt(1,0) *= temp;
									mt(2,0) *= temp;
								}

								{
									const std::complex<real_t> temp = mt(0,1) * conj(mt(0,0)) + mt(1,1) * conj(mt(1,0)) + mt(2,1) * conj(mt(2,0));
									mt(0,1) -= temp * mt(0,0);
									mt(1,1) -= temp * mt(1,0);
									mt(2,1) -= temp * mt(2,0);
								}

								{
									const std::complex<real_t> temp = 1.0 / sqrt(absq(mt(0,1)) + absq(mt(1,1)) + absq(mt(2,1)));
									mt(0,1) *= temp;
									mt(1,1) *= temp;
 									mt(2,1) *= temp;
								}

								mt(0,2) = conj(mt(1,0)) * conj(mt(2,1)) - conj(mt(2,0)) * conj(mt(1,1));
								mt(1,2) = conj(mt(2,0)) * conj(mt(0,1)) - conj(mt(0,0)) * conj(mt(2,1));
								mt(2,2) = conj(mt(0,0)) * conj(mt(1,1)) - conj(mt(1,0)) * conj(mt(0,1));
								environment.gaugeLinkConfiguration[localsite][mu] = mt;
							}
						}
						else {
							for (size_t mu = 0; mu < 4; ++mu) {
								for (size_t ii = 0; ii < 6; ++ii) {
									float tmp, tmp2;
									red += xdr_float(&xin, &tmp);
									red += xdr_float(&xin, &tmp2);
								}
							}
						}
#endif
					}
				}
			}
		}




		// - check if the number of read fields is correct - reading failure if not
		if (isOutputProcess()) {
			if (red != 48 * Layout::globalVolume) {
				std::cout << " read error for gauge fields on lattice " << ": " << red << std::endl;
				success = false;
			}
			else {
				std::cout << "ReadStartGaugeConfiguration::Finished reading " << red / 4 << " gauge fields" << std::endl;
				std::cout << "ReadStartGaugeConfiguration::From file " << filename << std::endl;
			}

			/// - read the global lattice extend and lvrco
          		int lato[4];

         		for (int ii(0); ii < 4; ++ii) xdr_int(&xin, &lato[ii]);

          		double beto(0.0);
          		xdr_double(&xin, &beto);

     			double c1(0.0);
          		xdr_double(&xin, &c1);

         		std::cout << "ReadStartGaugeConfiguration::Beta and c1 values (not checked) beta=" << beto << " c1=" << c1 << std::endl;

			// - compare read parameters to the ones in measureconfig.dat / Stop program if they are not the same (This is first done directly on the master node without communication of all the paramters)
			// - test for lattice size - new included: reading failure if given incorrect
			if (static_cast<int>(lato[0]) != Layout::glob_x
					|| static_cast<int>(lato[1]) != Layout::glob_y
					|| static_cast<int>(lato[2]) != Layout::glob_z
					|| static_cast<int>(lato[3]) != Layout::glob_t) {
				if (isOutputProcess()) std::cout << "Lattice size is different:\n" << lato[0] << " " << lato[1] << " "
						<< lato[2] << " " << lato[3] << "\n" << Layout::glob_x << " "
						<< Layout::glob_y << " " << Layout::glob_z << " "
						<< Layout::glob_t;
				success = false;
			}

			double plaq(0.0);
			if (xdr_double(&xin, &plaq)) {
				read_plaquette = plaq;
			}
			else {
				std::cout << "ReadStartGaugeConfiguration::Warning, no temporal plaquette value found!" << std::endl;
			}
			

		} // on of Io node only

		xdr_destroy( &xin);
		fclose(fin);
		
#endif
	}
	environment.gaugeLinkConfiguration.updateHalo();
	environment.synchronize();

	double plaquette = Plaquette::temporalPlaquette(environment.gaugeLinkConfiguration);

	if (isOutputProcess()) std::cout << "ReadStartGaugeConfiguration::Plaquette difference: " << (plaquette - read_plaquette) << std::endl;

	return success;


}

} /* namespace Update */
