/*
 * main.cpp
 *
 *  Created on: Feb 21, 2012
 *      Author: spiem_01
 */

#include "./source/Environment.h"
#include "./source/io/StorageParameters.h"
#include "./source/MatrixTypedef.h"
#include "./source/utils/RandomSeed.h"
#include "./source/utils/ToString.h"
#include "./source/Simulation.h"
#include "./source/LatticeSweep.h"
#include "./source/io/GlobalOutput.h"
#include "./source/MPILattice/ReducedStencil.h"
#include "./source/MPILattice/StandardStencil.h"
#include "./source/MPILattice/ExtendedStencil.h"
#include "./source/MPILattice/LocalLayout.h"
#include "./source/utils/LieGenerators.h"
#include "./source/utils/ToString.h"
#include <iostream>
#include <fenv.h>

int Update::RandomSeed::counter = -1;
boost::mt19937 Update::RandomSeed::rng;
boost::uniform_int<> Update::RandomSeed::dist = boost::uniform_int<>(-1000000,1000000);//TODO

#ifdef ENABLE_MPI
MPI_Datatype MpiType<short int>::type = MPI_SHORT;
MPI_Datatype MpiType<int>::type = MPI_INT;
MPI_Datatype MpiType<int[4]>::type = MPI_INT;
MPI_Datatype MpiType<Update::real_t>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::real_t[4]>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::FundamentalGroup[4]>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::AdjointGroup[4]>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::FundamentalGroup[6]>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::AdjointGroup[6]>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::FundamentalVector[4]>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::AdjointVector[4]>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::AdjointVector>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::FundamentalVector>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::FundamentalGroup>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::AdjointGroup>::type = MPI_DOUBLE;
#ifdef ADJOINT
MPI_Datatype MpiType<Update::FermionicForceMatrix[4]>::type = MPI_DOUBLE;
#endif
#endif


namespace po = boost::program_options;

template<int N>
inline bool testBasicMatrixFunctions(){
    typedef Eigen::Matrix<Update::complex, N, N> TMatrix;
    TMatrix tmp1(TMatrix::Zero());
    TMatrix tmp2(TMatrix::Zero());
    TMatrix tmp3(TMatrix::Zero());
    std::cout<< trace(tmp1*tmp2) <<std::endl;
    return true;
}


int main(int ac, char* av[]) {
	//feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
	/*Update::LieGenerator<Update::GaugeGroup> lieGenerator;
	Update::LieGenerator<Update::AdjointGroup> adjointLieGenerator;
	for (unsigned int t = 0; t < Update::numberColors*Update::numberColors - 1; ++t) std::cout << lieGenerator.get(t) << std::endl << std::endl;
	std::cout << std::endl << std::endl;
	for (unsigned int t = 0; t < Update::numberColors*Update::numberColors - 1; ++t) std::cout << adjointLieGenerator.get(t) << std::endl << std::endl;*/
	/*Update::AdjointGroup* test = new Update::AdjointGroup[3];
	test[0] = Update::AdjointGroup::Identity();
	set_to_zero(test[1]);
	test[2] = Update::AdjointGroup::Identity();

	typedef Update::GaugeVector vt4[4];

	vt4* vect = new vt4[11];
	for (unsigned int i = 0; i < 11; ++i) {
		for (unsigned int k = 0; k < 4; ++k)
		for (unsigned int j = 0; j < 3; ++j) {
			vect[i][k][j] = i+ j + 2*k;
		}
	}
	std::cout << "Qui ci arrivo!" << std::endl;
	std::complex<Update::real_t>* test2 = reinterpret_cast<std::complex<Update::real_t>*>(vect);
	//((Update::real_t (*)[5])test)[1][2]
	std::cout << ((std::complex<Update::real_t> (*)[4][3])test2)[7][2][2] << " " << vect[7][2][2] << std::endl;

	exit(0);
	Update::GaugeGroup test1, test2;
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			test1.at(i,j) = (static_cast<Update::real_t>(i+ 2*j +1 ))/7.;
			test2.at(i,j) = (static_cast<Update::real_t>(3*i- 2*j +1) )/2.;
		}
	}
	Update::GaugeGroup test3(test1*test2);
	Update::complex res1 = trace(test1*test2);
	Update::complex res2 = trace(test3);
	std::cout << "Test: " << res1 << " " << res2 << std::endl;
	std::cout << Update::toString(test1) << std::endl;
	std::cout << Update::toString(test2) << std::endl;
	std::cout << Update::toString(test3) << std::endl;

	testBasicMatrixFunctions<5>();*/

#ifdef ENABLE_MPI
	MPI_Init(&ac, &av);
#endif

	//The descriptor for the options
	po::options_description desc("Allowed options");
	desc.add_options()
	    ("help", "produce help message")
	    ("list_sweeps", "produce the list of all the available sweeps")
	    ("configfile", po::value<std::string>(), "set the configuration file")
	    ("beta", po::value<Update::real_t>(), "set the \\beta parameter of the simulations")
	    ("kappa", po::value<Update::real_t>(), "set the \\kappa parameter of the simulations")
	    ("mass", po::value<Update::real_t>(), "set the mass parameter of the overlap operator")
		("glob_x", po::value<unsigned int>(), "The x lattice size")
		("glob_y", po::value<unsigned int>(), "The y lattice size")
		("glob_z", po::value<unsigned int>(), "The z lattice size")
		("glob_t", po::value<unsigned int>(), "The t lattice size")
		("pgrid_x", po::value<unsigned int>(), "The grid subdivision in the x direction")
		("pgrid_y", po::value<unsigned int>(), "The grid subdivision in the y direction")
		("pgrid_z", po::value<unsigned int>(), "The grid subdivision in the z direction")
		("pgrid_t", po::value<unsigned int>(), "The grid subdivision in the t direction")
	    ("start", po::value<std::string>(), "the start method of the gauge configuration for the simulations (hotstart/coldstart/readstart)")
	    ("number_warm_up_sweeps", po::value<unsigned int>(), "the number of warm-up sweeps")
	    ("number_measurement_sweeps", po::value<unsigned int>(), "the number of measurement sweeps")
	    ("warm_up_sweeps", po::value< std::string >(), "the vector of the warm up sweeps to do (example: {{PureGaugeCM,1,1},{Plaquette,1,1}} )")
	    ("measurement_sweeps", po::value< std::string >(), "the vector of the measurement sweeps to do (example: {{PureGaugeCM,1,1},{Plaquette,1,1}} )")
	    ("name_action", po::value<std::string>(), "the name of the gauge part of the action (\"StandardWilson\")")
	    ("name_integrator", po::value<std::string>(), "the name of the type of integrator (\"second_order\")")
		("adjoint_nf_scalars", po::value<unsigned int>()->default_value(0), "set the number of the adjoint scalar fields")
		("fundamental_nf_scalars", po::value<unsigned int>()->default_value(0), "set the number of the fundamental scalar fields")
	    ("hmc_t_length", po::value<Update::real_t>(), "the length of a single HMC step (examples: 0.1, 0.05 ...)")
	    ("number_hmc_steps", po::value<std::string>(), "the vector of the numbers of HMC steps for a single trajectory (examples: 2, 7 ...)")
	    ("dirac_operator", po::value<std::string>(), "the name of the dirac wilson operator that should be used (supported: DiracWilson)")
	    ("number_pseudofermions", po::value<unsigned int>(), "the number of pseudofermions used")
	    ("number_force_levels", po::value<unsigned int>(), "the number of levels used for the force")
		("check_rational_approximations", po::value< std::string >(), "Set to true for checking the approximations in the beginning of the simulation")
		("twisted_mass_squared", po::value<Update::real_t>(), "The twisted mass squared used by twisted updater to regularize to RHMC")
		("number_twisted_correction_noise_vectors", po::value<unsigned int>(), "Number of noise vectors used to estimate the correction factor")
		("number_block_correction_noise_vectors", po::value<unsigned int>(), "Number of noise vectors used to estimate the correction factor")
		("theory_power_factor", po::value<Update::real_t>(), "The half of the number of fermion (1/4 for SUSY, 1 for two flavor, etc)")
		("twisted_breakup", po::value<unsigned int>(), "The value of the determinant breakup used in the correction step")
		("deflation_block_size", po::value<std::string>(), "The vector of the block size {lx,ly,lz,lt}")
		("t_source_origin", po::value<unsigned int>(), "The t-origin for the source of the pion and gluino correlators")
		("preconditioner_recursions", po::value<unsigned int>(), "Number of preconditioner recursions used for evaluating the rational approximations")
		("preconditioner_precision", po::value<Update::real_t>(), "The precision for the inversion of the preconditioner used for evaluating the rational approximations")
		("boundary_conditions", po::value<std::string>(), "Boundary conditions to use: periodic (fermions), antiperiodic (fermions), spatialantiperiodic (fermion), open")
		("measure_condensate_connected", po::value<bool>(), "Should I measure the connected part of the condensate susceptibility?")		
		("eigenvalues_map", po::value<std::string>(), "the map used for computing the lowest eigenvalues (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})")
	    ("output_directory_configurations", po::value<std::string>(), "The directory for the output of the configurations")
		("output_directory_measurements", po::value<std::string>(), "The directory for the output of the measurements")
	    ("output_name", po::value<std::string>(), "The name for the beginning part of the file of the output")
		("input_directory_configurations", po::value<std::string>(), "The directory for the input of the configurations")
	    ("input_name", po::value<std::string>(), "The name for the beginning part of the file of the input")
		("input_number", po::value<unsigned int>(), "The number of the file of the input")
	    ("number_threads", po::value<unsigned int>(), "The number of threads for openmp")
	    ("number_stochastic_estimators", po::value<unsigned int>(), "The number of stochastic estimators")
	    ("generic_inverter_precision", po::value<Update::real_t>(), "The precision for the inverter")
	    ("force_inverter_precision", po::value<Update::real_t>(), "The precision for the inverter in the force step")
	    ("metropolis_inverter_precision", po::value<Update::real_t>(), "The precision for the inverter in the metropolis step")
	    ("start_gauge_configuration_file", po::value<std::string>(), "The name of the beginning file with the gauge configuration for reading")
	    ("start_configuration_number", po::value<unsigned int>(), "The beginning number of the output configuration written (zero as default)")
	    ("output_configuration_name", po::value<std::string>(), "The name of the output of the field")
	    ("csw", po::value<Update::real_t>(), "The clover term coefficient")
	    ("version", "Version message")
	    ("momentum_operator", po::value<std::string>(),"\"true\" if the operators will be measured also with momentum")
	    ("maximum_momentum", po::value<unsigned int>(),"maximum momentum measured")
	    ("generic_inverter_max_steps", po::value<unsigned int>(),"maximum level of steps used by the inverters")
	    ("metropolis_inverter_max_steps", po::value<unsigned int>(),"maximum level of steps used by the inverters for computing the energy of the metropolis step")
	    ("force_inverter_max_steps", po::value<unsigned int>(),"maximum level of steps used by the inverters for computing the force")
	    ("max_r_wilsonloop", po::value<unsigned int>(),"maximum R for the wilson loops")
	    ("min_r_wilsonloop", po::value<unsigned int>(),"minimum R for the wilson loops")
	    ("max_t_wilsonloop", po::value<unsigned int>(),"maximum T for the wilson loops")
	    ("min_t_wilsonloop", po::value<unsigned int>(),"minimum T for the wilson loops")
		("number_subsweeps_luescher", po::value<unsigned int>(),"number of subsweeps of the luescher algorithm")
		("size_slice_luescher", po::value<unsigned int>(),"the size of a slice of the luescher algorithm")
	    ("level_stout_smearing_meson", po::value<unsigned int>(), "The number of levels for stout smearing in measuring meson masses")
		("level_stout_smearing_glueball", po::value<unsigned int>(), "The number of levels for stout smearing in measuring glueball masses")
		("level_stout_smearing_gluinoglue", po::value<unsigned int>(), "The number of levels for stout smearing in measuring gluinoglue masses")
		("level_stout_smearing_wilson_loop", po::value<unsigned int>(), "The number of levels for stout smearing in measuring the wilson loop")
		("level_stout_smearing_polyakov_correlator", po::value<unsigned int>(), "The number of levels for stout smearing in measuring the Polyakov loop correlator")
	    ("rho_stout_smearing", po::value<Update::real_t>(), "The rho parameter for stout smearing in measuring masses")
	    ("calculate_disconnected_contributions", po::value<bool>(), "Should the program calculate disconnected contributions?")
	    ("number_multiplication_test_speed", po::value<unsigned int>(), "How many multiplications should I use in the tests?")
		("output_offset", po::value<unsigned int>(), "The offset for the number of the output configurations")
		("load_layout", "If the layout should be loaded from the disk")
		("print_report_layout", "If the full report of the layout should be printed")
		("format_name", po::value<std::string>(), "leonard_format/muenster_format for reading and writing configurations")
		("input_format_name", po::value<std::string>(), "leonard_format/muenster_format only for reading configurations")
		("output_format_name", po::value<std::string>(), "leonard_format/muenster_format only for writing configurations")
		("OverlapOperator::squareRootApproximation", po::value<std::string>(), "Approximation of x^(1/2) used for the Overlap fermion sign function (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})")
		("ExactOverlapOperator::squareRootApproximation", po::value<std::string>()->default_value(""), "Approximation of x^(1/2) used for the Overlap fermion sign function (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})")
		("ExactOverlapOperator::eigensolver::use_chebyshev", po::value<std::string>()->default_value("false"), "Use Chebyshev acceleration? (true/false)")
		("ExactOverlapOperator::eigensolver::chebyshev_left", po::value<double>()->default_value(0.2), "Left interval of the Chebyshev polynomial.")
		("ExactOverlapOperator::eigensolver::chebyshev_right", po::value<double>()->default_value(7.), "Right interval of the Chebyshev polynomial.")
		("ExactOverlapOperator::eigensolver::chebyshev_order", po::value<unsigned int>()->default_value(15), "Order of the Chebyshev acceleration. It must be an odd number")
		("ExactOverlapOperator::eigensolver::eigensolver_precision", po::value<double>()->default_value(0.000000001), "set the precision used by the eigensolver")
		("ExactOverlapOperator::eigensolver::number_extra_vectors", po::value<unsigned int>(), "Number of extra vectors for the Arnoldi algorithm used in the computation of the eigenvectors, increase this number to increase precision")
		("ExactOverlapOperator::eigensolver::maximal_number_restarts_eigensolver", po::value<unsigned int>()->default_value(50), "Number of restarts for the implicitly restarted Arnoldi algorithm")
		("ExactOverlapOperator::eigensolver::number_eigenvalues", po::value<unsigned int>(), "Number of eigenvalues of the dirac wilson operator to be computed")
		("correction_step_breakup_level", po::value<unsigned int>(), "The level of breakup \"lb\" for the correction factors")
		("correction_approximation_direct_1", po::value<std::string>(), "Approximation of x^(nf/(2lb)) used for the correction step (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})")
		("correction_approximation_inverse_1", po::value<std::string>(), "Approximation of x^(-nf/(2lb)) used for the correction step (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})")
		("correction_approximation_direct_2", po::value<std::string>(), "Approximation of x^(nf/2) used for the correction step (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})")
		("correction_approximation_inverse_2", po::value<std::string>(), "Approximation of x^(-nf/(2lb)) used for the correction step (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})")
		("correction_approximation_direct_3", po::value<std::string>(), "Approximation of x^(nf/(2lb)) used for the correction step (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})")
		("correction_approximation_inverse_3", po::value<std::string>(), "Approximation of x^(-nf/(2lb)) used for the correction step (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})")
		("correction_approximation_direct_4", po::value<std::string>(), "Approximation of x^(nf/(2lb)) used for the correction step (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})")
		("correction_approximation_inverse_4", po::value<std::string>(), "Approximation of x^(-nf/(2lb)) used for the correction step (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})")
		("read_start_number", po::value<unsigned int>(), "From which configurations start to read again for the analysis")
		("read_step", po::value<unsigned int>(), "The step in the reading analysis")
		("flow_type", po::value<std::string>(), "The type of flow equations (Wilson/Symanzik)")
		("flow_step", po::value<Update::real_t>(), "The integration step of the flow")
		("flow_integration_intervals", po::value<unsigned int>(), "The number of intervals for each flow integration step")
		("flow_time", po::value<Update::real_t>(), "The total time of the flow")
		("stout_smearing_levels", po::value<int>(), "The levels for the stout smearing of the dirac operator")
		("stout_smearing_rho", po::value<Update::real_t>(), "The rho for the stout smearing of the dirac operator")
	;

	for (int level = 1; level < 4; ++level) {
		for (int i = 1; i < 33; ++i) {
			desc.add_options()((std::string("force_rational_fraction_")+Update::toString(i)+"_level_"+Update::toString(level)).c_str(),po::value<std::string>(),"Rational approximation of the force used by MultiStep integrator");
		}
		desc.add_options()(std::string("force_inverter_precision_level_"+Update::toString(level)).c_str(),po::value<Update::real_t>(),"Precision of the inverter of the force used by MultiStep integrator");
	}
	for (int i = 1; i < 33; ++i) {
		desc.add_options()((std::string("heatbath_rational_fraction_")+Update::toString(i)).c_str(),po::value<std::string>(),"the rational fraction approximation that should be used for the heatbath (syntax: {alpha_1,..,alpha_n,beta_1, ..., beta_n})");
		desc.add_options()((std::string("metropolis_rational_fraction_")+Update::toString(i)).c_str(),po::value<std::string>(),"the rational fraction approximation that should be used for the metropolis (syntax: {alpha_1,..,alpha_n,beta_1, ..., beta_n})");
		desc.add_options()((std::string("force_rational_fraction_")+Update::toString(i)).c_str(),po::value<std::string>(),"the rational fraction approximation that should be used for the force (syntax: {alpha_1,..,alpha_n,beta_1, ..., beta_n})");
	}

	Update::LatticeSweep::addParameters(desc);

	//Now we can parse the command line
	po::variables_map vm;
	po::store(po::parse_command_line(ac, av, desc), vm);
	po::notify(vm);

	//If there is a configfile we parse also it
	if (vm.count("configfile")) {
		std::ifstream in(vm["configfile"].as<std::string>().c_str());
		po::store(po::parse_config_file(in, desc, true), vm);
		po::notify(vm);
	} else if (vm.count("help")) {
		//Print the help description of the parameters of the program
		std::cout << desc << std::endl;
		std::cout << "type --list_sweeps for having a list of all the possible sweeps" << std::endl;
		return 0;
	} else if (vm.count("version")) {
		std::cout << "Updater Muenster version " << 1932 << std::endl;
#ifdef ADJOINT
		std::cout << " with adjoint fermions" << std::endl;
#endif
#ifndef ADJOINT
		std::cout << " with fundamental fermions" << std::endl;
#endif
		std::cout << " with number of colors: " << Update::numberColors << std::endl;
#ifdef ENABLE_MPI
		std::cout << " with mpi enabled" << std::endl;
#endif
#ifdef TESTLATTICE
		std::cout << " with test lattice";
#endif
#ifndef TESTLATTICE
		std::cout << " with standard lattice";
#endif
		std::cout << " of size: " << 5 << "^3 x " << 3 << std::endl;
		return 0;
	} else if (vm.count("list_sweeps")) {
		Update::LatticeSweep::printSweepsName();
		return 0;
	}

#ifndef TESTLATTICE
	//Initialize the layout of the program
	//Math::StandardLayout::startup_initialize(ac, av);
	//Math::StandardLayout::Env::gcor().
#endif

	Lattice::StandardStencil::initializeNeighbourSites();
	Lattice::ExtendedStencil::initializeNeighbourSites();
	Lattice::ReducedStencil::initializeNeighbourSites();

#ifndef ENABLE_MPI
	Lattice::LocalLayout::pgrid_t = 1;
	Lattice::LocalLayout::pgrid_x = 1;
	Lattice::LocalLayout::pgrid_y = 1;
	Lattice::LocalLayout::pgrid_z = 1;

	Lattice::LocalLayout::glob_t = vm["glob_t"].as<unsigned int>();
	Lattice::LocalLayout::glob_x = vm["glob_x"].as<unsigned int>();
	Lattice::LocalLayout::glob_y = vm["glob_y"].as<unsigned int>();
	Lattice::LocalLayout::glob_z = vm["glob_z"].as<unsigned int>();

	Lattice::LocalLayout::initialize();

	if (Update::isOutputProcess()) std::cout << "Lattice size (x,y,z,t): (" << Lattice::LocalLayout::glob_x << "," << Lattice::LocalLayout::glob_y << "," << Lattice::LocalLayout::glob_z << "," << Lattice::LocalLayout::glob_t << ")" << std::endl;
#endif
#ifdef ENABLE_MPI
	Lattice::MpiLayout<Lattice::ExtendedStencil>::pgrid_t = vm["pgrid_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::pgrid_x = vm["pgrid_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::pgrid_y = vm["pgrid_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::pgrid_z = vm["pgrid_z"].as<unsigned int>();
	
	Lattice::MpiLayout<Lattice::StandardStencil>::pgrid_t = vm["pgrid_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::pgrid_x = vm["pgrid_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::pgrid_y = vm["pgrid_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::pgrid_z = vm["pgrid_z"].as<unsigned int>();

	Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_t = vm["pgrid_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_x = vm["pgrid_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_y = vm["pgrid_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_z = vm["pgrid_z"].as<unsigned int>();


	Lattice::MpiLayout<Lattice::ExtendedStencil>::glob_t = vm["glob_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::glob_x = vm["glob_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::glob_y = vm["glob_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::glob_z = vm["glob_z"].as<unsigned int>();
	
	Lattice::MpiLayout<Lattice::StandardStencil>::glob_t = vm["glob_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::glob_x = vm["glob_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::glob_y = vm["glob_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::glob_z = vm["glob_z"].as<unsigned int>();

	Lattice::MpiLayout<Lattice::ReducedStencil>::glob_t = vm["glob_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::glob_x = vm["glob_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::glob_y = vm["glob_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::glob_z = vm["glob_z"].as<unsigned int>();

	if (vm.count("load_layout")) {
		std::string basename = vm["output_directory_configurations"].as<std::string>();
		Lattice::MpiLayout<Lattice::ExtendedStencil>::load(basename+"extended_layout");
		Lattice::MpiLayout<Lattice::StandardStencil>::load(basename+"standard_layout");
		Lattice::MpiLayout<Lattice::ReducedStencil>::load(basename+"reduced_layout");
		//Lattice::MpiLayout<Lattice::ReducedUpStencil>::load(basename+"reduced_up_layout");
		//Lattice::MpiLayout<Lattice::ReducedDownStencil>::load(basename+"reduced_down_layout");
	} else {
		Lattice::MpiLayout<Lattice::ExtendedStencil>::initialize();
		Lattice::MpiLayout<Lattice::StandardStencil>::initialize();
		Lattice::MpiLayout<Lattice::ReducedStencil>::initialize();
		//Lattice::MpiLayout<Lattice::ReducedUpStencil>::initialize();
		//Lattice::MpiLayout<Lattice::ReducedDownStencil>::initialize();
		
		std::string basename = vm["output_directory_configurations"].as<std::string>();
		Lattice::MpiLayout<Lattice::ExtendedStencil>::save(basename+"extended_layout");
		Lattice::MpiLayout<Lattice::StandardStencil>::save(basename+"standard_layout");
		Lattice::MpiLayout<Lattice::ReducedStencil>::save(basename+"reduced_layout");
		//Lattice::MpiLayout<Lattice::ReducedUpStencil>::save(basename+"reduced_up_layout");
		//Lattice::MpiLayout<Lattice::ReducedDownStencil>::save(basename+"reduced_down_layout");
	}

	if (vm.count("print_report_layout")) {
		if (Update::isOutputProcess()) Lattice::MpiLayout<Lattice::ExtendedStencil>::printReport();
		if (Update::isOutputProcess()) Lattice::MpiLayout<Lattice::ReducedStencil>::printReport();
	}

	if (Update::isOutputProcess()) std::cout << "Lattice size (x,y,z,t): (" << Lattice::MpiLayout<Lattice::ReducedStencil>::glob_x << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::glob_y << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::glob_z << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::glob_t << ")" << std::endl;
	if (Update::isOutputProcess()) std::cout << "Mpi grid (px,py,pz,pt): (" << Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_x << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_y << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_z << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_t << ")" << std::endl;
#endif

	//Initialize the enviroment of the program
	Update::environment_t* environment = new Update::environment_t(vm);

	//Set the number of threads
	if (vm.count("number_threads")) {
#ifdef MULTITHREADING
		omp_set_num_threads(vm["number_threads"].as<unsigned int>());
		int num_threads = omp_get_max_threads();
		if (Update::isOutputProcess()) std::cout << "Number of threads: " << num_threads << std::endl;
#endif
	}

	Update::Simulation simulation(*environment);
	simulation.starter();
	simulation.warmUp();
	simulation.measurement();

	Update::GlobalOutput* output = Update::GlobalOutput::getInstance();
	output->print();
	output->destroy();

	delete environment;

#ifdef ENABLE_MPI
	Lattice::MpiLayout<Lattice::ExtendedStencil>::destroy();
	Lattice::MpiLayout<Lattice::StandardStencil>::destroy();
	Lattice::MpiLayout<Lattice::ReducedStencil>::destroy();
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

	return 0;
}


