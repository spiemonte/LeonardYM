#include "Environment.h"
#include "program_options/StorageParameters.h"
#include "MatrixTypedef.h"
#include "utils/RandomSeed.h"
#include "utils/ToString.h"
#include "Simulation.h"
#include "LatticeSweep.h"
#include "io/GlobalOutput.h"
#include "MPILattice/ReducedStencil.h"
#include "MPILattice/StandardStencil.h"
#include "MPILattice/ExtendedStencil.h"
#include "MPILattice/LocalLayout.h"
#include "utils/LieGenerators.h"
#include "utils/ToString.h"
#include <iostream>
#include <fenv.h>
#include <regex>

//MPI Datatype initilialization
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
MPI_Datatype MpiType<Update::AdjointRealVector[4]>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::AdjointRealVector>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::AdjointComplexVector[4]>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::AdjointComplexVector>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::FundamentalVector>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::FundamentalGroup>::type = MPI_DOUBLE;
MPI_Datatype MpiType<Update::AdjointGroup>::type = MPI_DOUBLE;
#ifdef ADJOINT
MPI_Datatype MpiType<Update::FermionicForceMatrix[4]>::type = MPI_DOUBLE;
#endif
#endif

std::string clean_line(const std::string& input) {
	std::string result = input;
	if (result.rfind("#",0) != std::string::npos) {
		result = result.substr(0, result.rfind("#",0));
	}
	result.erase(std::remove(result.begin(), result.end(), ' '), result.end());
	return result;
}

std::pair<std::string, std::string> parse_parameter(const std::string& input) {
    std::regex re("[a-zA-Z0-9\\.:,{}\\ _\\-/]+");
    std::sregex_iterator
        first(input.begin(), input.end(), re),
        last;
    
    std::vector<std::string> matches;
    for (std::sregex_iterator i = first; i != last; ++i) {
    	std::string m = i->str();
    	m.erase(std::remove(m.begin(), m.end(), ' '), m.end());
    	matches.push_back(m);
    }

    if (matches.size() == 2) {
    	return std::pair<std::string, std::string>(matches[0], matches[1]);
    }
    if (matches.size() == 1) {
    	return std::pair<std::string, std::string>(matches[0], "");
    }

    std::cout << "Wrongly specified option " << input << matches.size() << std::endl;
   	exit(77);
}


int main(int ac, char* av[]) {
#ifdef ENABLE_MPI
	//Initialize MPI if in MPI mode
	MPI_Init(&ac, &av);
#endif

	//The map descriptor for the global options
	std::map<std::string, Update::Option> desc;

	//Action specifications
	desc["name_action"]	= Update::Option("name_action", "Improved", "the name of the gauge part of the action (\"StandardWilson/Improved\")");
	desc["dirac_operator"] = Update::Option("dirac_operator", "DiracWilson", "the name of the dirac wilson operator that should be used (supported: DiracWilson/Improved/Overlap/ExactOverlap)");
		
	//Action and dirac operator options
	desc["beta"] = Update::Option("beta", 2.0, "set the \\beta parameter of the simulations");
	desc["kappa"] = Update::Option("kappa", 0.125, "set the \\kappa parameter of the simulations");
	desc["mass"] = Update::Option("mass", 0.0, "set the mass parameter of the overlap operator");
	desc["csw"] =	Update::Option("csw", 1.0, "The clover term coefficient");
	desc["stout_smearing_levels"] = Update::Option("stout_smearing_levels", 0, "The levels for the stout smearing of the dirac operator");
	desc["stout_smearing_rho"] = Update::Option("stout_smearing_rho", 0.15, "The rho for the stout smearing of the dirac operator");
		
	//GRID parallelization and lattice options
	desc["glob_x"] = Update::Option("glob_x", 4, "The x lattice size");
	desc["glob_y"] = Update::Option("glob_y", 4, "The y lattice size");
	desc["glob_z"] = Update::Option("glob_z", 4, "The z lattice size");
	desc["glob_t"] = Update::Option("glob_t", 4, "The t lattice size");
	desc["pgrid_x"] = Update::Option("pgrid_x", 1, "The grid subdivision in the x direction");
	desc["pgrid_y"] = Update::Option("pgrid_y", 1, "The grid subdivision in the y direction");
	desc["pgrid_z"] = Update::Option("pgrid_z", 1, "The grid subdivision in the z direction");
	desc["pgrid_t"] = Update::Option("pgrid_t", 1, "The grid subdivision in the t direction");
	desc["number_threads"] = Update::Option("number_threads", 2, "The number of threads for openmp");
	desc["load_layout"] = Update::Option("load_layout", "false", "If the MPI layout should be loaded from the disk");
	desc["print_report_layout"]	= Update::Option("print_report_layout", "false", "If the full report of the MPI layout should be printed");

	//Boundary conditions
	desc["boundary_conditions"] = Update::Option("boundary_conditions", "fermion_periodic", "Boundary conditions to use: periodic (fermions), antiperiodic (fermions), spatialantiperiodic (fermion), open");

	//Input-output options
	desc["output_directory_configurations"] = Update::Option("output_directory_configurations", "./", "The directory for the output of the configurations");
	desc["output_directory_measurements"] = Update::Option("output_directory_measurements", "./", "The directory for the output of the measurements");
	desc["output_name"] = Update::Option("output_name", "measurements", "The name for the beginning part of the file of the output");
	desc["input_directory_configurations"] = Update::Option("input_directory_configurations", "./", "The directory for the input of the configurations");
	desc["input_name"] = Update::Option("input_name", "gauge_fields", "The name for the beginning part of the file of the input");
	desc["input_number"] = Update::Option("input_number", 1, "The number of the file of the input");
	desc["output_configuration_name"] = Update::Option("output_configuration_name", "gauge_fields", "The name of the output of the field");
	desc["output_offset"] = Update::Option("output_offset", 1, "The offset for the number of the output configurations");
	//desc["format_name"]	= Update::Option("format_name", "muester_format", "leonard_format/muenster_format for reading and writing configurations");
	desc["input_format_name"] = Update::Option("input_format_name", "muester_format", "leonard_format/muenster_format only for reading configurations");
	desc["output_format_name"] = Update::Option("output_format_name", "muenster_format", "leonard_format/muenster_format only for writing configurations");
	desc["measurement_output_format"] = Update::Option("measurement_output_format", "txt", "output format for the measurements (xml/txt)");
		
	//Start, warm up and measurement specifications
	desc["start"] = Update::Option("start", "hotstart", "the start gauge configuration for the simulations (hotstart/coldstart/readstart)");
	desc["start_gauge_configuration_file"] = Update::Option("start_gauge_configuration_file", "gauge_fields", "The name of the beginning file with the gauge configuration for reading");
	desc["start_configuration_number"] = Update::Option("start_configuration_number", 1, "The beginning number of the output configuration written (zero as default)");
	desc["number_warm_up_sweeps"] = Update::Option("number_warm_up_sweeps", 1, "the number of warm-up sweeps");
	desc["number_measurement_sweeps"] = Update::Option("number_measurement_sweeps", 1, "the number of measurement sweeps");
	desc["warm_up_sweeps"] = Update::Option("warm_up_sweeps", "{{TestCommunications,1,1},{PureGaugeCM,1,1},{Plaquette,1,1}}", "the vector of the warm up sweeps to do (example: {{PureGaugeCM,1,1},{Plaquette,1,1}} )");
	desc["measurement_sweeps"] = Update::Option("measurement_sweeps", "{{PureGaugeCM,1,1},{Plaquette,1,1}}", "the vector of the measurement sweeps to do (example: {{PureGaugeCM,1,1},{Plaquette,1,1}} )");
		
	//Scalar field options
	desc["adjoint_nf_scalars"] = Update::Option("adjoint_nf_scalars", 0, "set the number of the adjoint scalar fields");
	desc["fundamental_nf_scalars"] = Update::Option("fundamental_nf_scalars", 0, "set the number of the fundamental scalar fields");
		
	//HMC options
	desc["name_integrator"] = Update::Option("name_integrator", "omelyan", "the name of the type of integrator (\"second_order, omelyan, fourth_order, fourth_omelyan\")");
	desc["hmc_t_length"] = Update::Option("hmc_t_length", 1.0, "the length of a single HMC step (examples: 0.1, 0.05 ...)");
	desc["number_hmc_steps"] = Update::Option("number_hmc_steps", "{1,1,1}", "the vector of the numbers of HMC steps for a single trajectory (examples: 2, 7 ...)");
		
	//RHMC options
	desc["force_inverter_precision"] = Update::Option("force_inverter_precision", 1e-11, "The precision for the inverter in the force step");
	desc["metropolis_inverter_precision"] = Update::Option("metropolis_inverter_precision", 1e-13, "The precision for the inverter in the metropolis step");
	desc["metropolis_inverter_max_steps"] = Update::Option("metropolis_inverter_max_steps", 50000, "maximum number of steps used by the inverters for computing the energy of the metropolis step");
	desc["force_inverter_max_steps"] = Update::Option("force_inverter_max_steps", 20000, "maximum number of steps used by the inverters for computing the force");
	desc["number_pseudofermions"] = Update::Option("number_pseudofermions", 2, "the number of pseudofermions used");
	desc["number_force_levels"] = Update::Option("number_force_levels", 2, "the number of levels used for the force");
	desc["check_rational_approximations"] = Update::Option("check_rational_approximations", "true", "Set to true for checking the approximations in the beginning of the simulation");
	desc["twisted_mass_squared"] = Update::Option("twisted_mass_squared", 0.0, "The twisted mass squared used by twisted updater to regularize to RHMC");
	desc["number_twisted_correction_noise_vectors"] = Update::Option("number_twisted_correction_noise_vectors", 20, "Number of noise vectors used to estimate the correction factor");
	desc["number_block_correction_noise_vectors"] = Update::Option("number_block_correction_noise_vectors", 20, "Number of noise vectors used to estimate the correction factor");
	desc["theory_power_factor"] = Update::Option("theory_power_factor", 0.25, "The half of the number of fermion (1/4 for SUSY, 1 for two flavor, etc)");
	desc["twisted_breakup"] = Update::Option("twisted_breakup", 2, "The value of the determinant breakup used in the correction step");
	desc["deflation_block_size"] = Update::Option("deflation_block_size", "{4,4,4,4}", "The vector of the block size {lx,ly,lz,lt}");
	
	//Options for the ReadGaugeConfiguration sweep
	desc["read_start_number"] = Update::Option("read_start_number", 1, "From which configurations start to read again for the analysis");
	desc["read_step"] = Update::Option("read_step", 1, "The step in the reading analysis");

	//Overlap operator options
	desc["OverlapOperator::squareRootApproximation"] = Update::Option("OverlapOperator::squareRootApproximation", "", "Approximation of x^(1/2) used for the Overlap fermion sign function (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})");
	desc["ExactOverlapOperator::squareRootApproximation"] = Update::Option("ExactOverlapOperator::squareRootApproximation", "", "Approximation of x^(1/2) used for the Overlap fermion sign function (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})");
	desc["ExactOverlapOperator::eigensolver::use_chebyshev"] = Update::Option("ExactOverlapOperator::eigensolver::use_chebyshev", "false", "Use Chebyshev acceleration? (true/false)");
	desc["ExactOverlapOperator::eigensolver::chebyshev_left"] = Update::Option("ExactOverlapOperator::eigensolver::chebyshev_left", 0.2, "Left interval of the Chebyshev polynomial.");
	desc["ExactOverlapOperator::eigensolver::chebyshev_right"] = Update::Option("ExactOverlapOperator::eigensolver::chebyshev_right", 7.0, "Right interval of the Chebyshev polynomial.");
	desc["ExactOverlapOperator::eigensolver::chebyshev_order"] = Update::Option("ExactOverlapOperator::eigensolver::chebyshev_order", 15, "Order of the Chebyshev acceleration. It must be an odd number");
	desc["ExactOverlapOperator::eigensolver::eigensolver_precision"] = Update::Option("ExactOverlapOperator::eigensolver::eigensolver_precision", 1e-9, "set the precision used by the eigensolver");
	desc["ExactOverlapOperator::eigensolver::number_extra_vectors"] = Update::Option("ExactOverlapOperator::eigensolver::number_extra_vectors", 30, "Number of extra vectors for the Arnoldi algorithm used in the computation of the eigenvectors, increase this number to increase precision");
	desc["ExactOverlapOperator::eigensolver::maximal_number_restarts_eigensolver"] = Update::Option("ExactOverlapOperator::eigensolver::maximal_number_restarts_eigensolver", 50, "Number of restarts for the implicitly restarted Arnoldi algorithm");
	desc["ExactOverlapOperator::eigensolver::number_eigenvalues"] = Update::Option("ExactOverlapOperator::eigensolver::number_eigenvalues", 30, "Number of eigenvalues of the Dirac-Wilson operator to be computed");

	//Options for the inverter
	desc["generic_inverter_max_steps"] = Update::Option("generic_inverter_max_steps", 20000,"maximum level of steps used by the inverters");
	desc["generic_inverter_precision"] = Update::Option("generic_inverter_precision", 1e-11, "The precision for the inverter");

	for (int level = 1; level < 4; ++level) {
		for (int i = 1; i < 33; ++i) {
			std::string option_name = std::string("force_rational_fraction_")+Update::toString(i)+"_level_"+Update::toString(level);
			desc[option_name] = Update::Option(option_name, "","Rational approximation of the force used by MultiStep integrator");
		}
		std::string option_name = std::string("force_inverter_precision_level_")+Update::toString(level);
		desc[option_name] = Update::Option(option_name, "","Precision of the inverter of the force used by MultiStep integrator");
	}
	for (int i = 1; i < 33; ++i) {
		std::string option_name = std::string("heatbath_rational_fraction_")+Update::toString(i);
		desc[option_name] = Update::Option(option_name, "", "the rational fraction approximation that should be used for the heatbath (syntax: {alpha_1,..,alpha_n,beta_1, ..., beta_n})");
		option_name = std::string("metropolis_rational_fraction_")+Update::toString(i);
		desc[option_name] = Update::Option(option_name, "", "the rational fraction approximation that should be used for the metropolis (syntax: {alpha_1,..,alpha_n,beta_1, ..., beta_n})");
		option_name = std::string("force_rational_fraction_")+Update::toString(i);
		desc[option_name] = Update::Option(option_name, "", "the rational fraction approximation that should be used for the force (syntax: {alpha_1,..,alpha_n,beta_1, ..., beta_n})");
	}

	Update::LatticeSweep::addParameters(desc);

	//Now we can parse the command line
	for (int i = 0; i < ac; ++i) {
		//If there is a configfile we parse also it
		if (std::string(av[i]).rfind("--configfile", 0) == 0) {
			std::string filename = parse_parameter(av[i]).second;
			std::ifstream file(filename);
			if (file.is_open()) {
  				std::string line;
    			while (std::getline(file, line)) {
    				line = clean_line(line);
    				if (line.size() != 0) {
        				std::pair<std::string, std::string> parameter = parse_parameter(line);
        				if (desc.count(parameter.first)) {
        					desc[parameter.first] = parameter.second;
        				}
        				else {
        					std::cout << "Warning, illegal option " << parameter.first << "!" << std::endl;
        				}
        			}
    			}
    			file.close();

    			std::cout << "Options set to:" << std::endl;
    			for (auto elem: desc) {
					std::cout << elem.second.toString() << std::endl;
				}
			}
		} else if (std::string(av[i]) == "--help" || std::string(av[i]) == "help") {
			//Print the help description of the parameters of the program
			for (auto elem: desc) {
				std::cout << elem.second.getHelp() << std::endl;
			}
			
			std::cout << "type --list_sweeps for having a list of all the possible sweeps" << std::endl;
			return 0;
		} else if (std::string(av[i]) == "version" || std::string(av[i]) == "--version") {
			std::cout << "LeonardQCD version: 1.1" << std::endl;
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
			return 0;
		} else if (std::string(av[i]) == "list_sweeps" || std::string(av[i]) == "--list_sweeps") {
			Update::LatticeSweep::printSweepsName();
			return 0;
		}
	}

	//Initialized up/down stencil
	Lattice::StandardStencil::initializeNeighbourSites();
	Lattice::ExtendedStencil::initializeNeighbourSites();
	Lattice::ReducedStencil::initializeNeighbourSites();

	//Initialize lattice layout
#ifndef ENABLE_MPI
	Lattice::LocalLayout::pgrid_t = 1;
	Lattice::LocalLayout::pgrid_x = 1;
	Lattice::LocalLayout::pgrid_y = 1;
	Lattice::LocalLayout::pgrid_z = 1;

	Lattice::LocalLayout::glob_t = desc["glob_t"].as<unsigned int>();
	Lattice::LocalLayout::glob_x = desc["glob_x"].as<unsigned int>();
	Lattice::LocalLayout::glob_y = desc["glob_y"].as<unsigned int>();
	Lattice::LocalLayout::glob_z = desc["glob_z"].as<unsigned int>();

	Lattice::LocalLayout::initialize();

	if (isOutputProcess()) std::cout << "Lattice size (x,y,z,t): (" << Lattice::LocalLayout::glob_x << "," << Lattice::LocalLayout::glob_y << "," << Lattice::LocalLayout::glob_z << "," << Lattice::LocalLayout::glob_t << ")" << std::endl;
#endif
#ifdef ENABLE_MPI
	Lattice::MpiLayout<Lattice::ExtendedStencil>::pgrid_t = desc["pgrid_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::pgrid_x = desc["pgrid_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::pgrid_y = desc["pgrid_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::pgrid_z = desc["pgrid_z"].as<unsigned int>();
	
	Lattice::MpiLayout<Lattice::StandardStencil>::pgrid_t = desc["pgrid_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::pgrid_x = desc["pgrid_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::pgrid_y = desc["pgrid_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::pgrid_z = desc["pgrid_z"].as<unsigned int>();

	Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_t = desc["pgrid_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_x = desc["pgrid_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_y = desc["pgrid_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_z = desc["pgrid_z"].as<unsigned int>();


	Lattice::MpiLayout<Lattice::ExtendedStencil>::glob_t = desc["glob_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::glob_x = desc["glob_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::glob_y = desc["glob_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ExtendedStencil>::glob_z = desc["glob_z"].as<unsigned int>();
	
	Lattice::MpiLayout<Lattice::StandardStencil>::glob_t = desc["glob_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::glob_x = desc["glob_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::glob_y = desc["glob_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::StandardStencil>::glob_z = desc["glob_z"].as<unsigned int>();

	Lattice::MpiLayout<Lattice::ReducedStencil>::glob_t = desc["glob_t"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::glob_x = desc["glob_x"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::glob_y = desc["glob_y"].as<unsigned int>();
	Lattice::MpiLayout<Lattice::ReducedStencil>::glob_z = desc["glob_z"].as<unsigned int>();

	if (desc["load_layout"].as<std::string>() == "true") {
		std::string basename = desc["output_directory_configurations"].as<std::string>();
		Lattice::MpiLayout<Lattice::ExtendedStencil>::load(basename+"extended_layout");
		Lattice::MpiLayout<Lattice::StandardStencil>::load(basename+"standard_layout");
		Lattice::MpiLayout<Lattice::ReducedStencil>::load(basename+"reduced_layout");
	} else {
		Lattice::MpiLayout<Lattice::ExtendedStencil>::initialize();
		Lattice::MpiLayout<Lattice::StandardStencil>::initialize();
		Lattice::MpiLayout<Lattice::ReducedStencil>::initialize();
		
		std::string basename = desc["output_directory_configurations"].as<std::string>();
		Lattice::MpiLayout<Lattice::ExtendedStencil>::save(basename+"extended_layout");
		Lattice::MpiLayout<Lattice::StandardStencil>::save(basename+"standard_layout");
		Lattice::MpiLayout<Lattice::ReducedStencil>::save(basename+"reduced_layout");
	}

	if (desc["print_report_layout"].as<std::string>() == "true") {
		if (isOutputProcess()) Lattice::MpiLayout<Lattice::ExtendedStencil>::printReport();
		if (isOutputProcess()) Lattice::MpiLayout<Lattice::ReducedStencil>::printReport();
	}

	if (isOutputProcess()) std::cout << "Lattice size (x,y,z,t): (" << Lattice::MpiLayout<Lattice::ReducedStencil>::glob_x << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::glob_y << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::glob_z << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::glob_t << ")" << std::endl;
	if (isOutputProcess()) std::cout << "Mpi grid (px,py,pz,pt): (" << Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_x << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_y << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_z << "," << Lattice::MpiLayout<Lattice::ReducedStencil>::pgrid_t << ")" << std::endl;
#endif

	//Initialize the enviroment of the program
	Update::environment_t* environment = new Update::environment_t(desc);

	//Set the number of threads
	if (desc["number_threads"].as<unsigned int>() > 0) {
#ifdef MULTITHREADING
		omp_set_num_threads(desc["number_threads"].as<unsigned int>());
		int num_threads = omp_get_max_threads();
		if (isOutputProcess()) std::cout << "Number of threads: " << num_threads << std::endl;
#endif
	}

	//Set the output to format
	Update::GlobalOutput* output = Update::GlobalOutput::getInstance();
	output->setFormat(desc["measurement_output_format"].as<std::string>());

	//Finally create the simulation
	Update::Simulation simulation(*environment);
	//Start and run warmup
	simulation.starter();
	simulation.warmUp();
	//Perform the measurements
	simulation.measurement();

	//Print and finalize the output to file
	output->print();
	output->destroy();

	//destroy the environment
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


