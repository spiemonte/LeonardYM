#ifndef STORAGEPARAMETERS_H_
#define STORAGEPARAMETERS_H_
#include <string>
#include <boost/program_options.hpp>
#include <iostream>
#include <vector>
#include <complex>
#include <sstream>
#include <exception>
#include <stdexcept>
#include "utils/FromString.h"

namespace Update {

class NotFoundOption : public std::runtime_error {
public:
	NotFoundOption(const std::string& _nameOption) throw () : std::runtime_error((_nameOption+": option not found!").c_str()) { }
	~NotFoundOption() throw () { };
};

namespace implement {

template<typename T> T get(const boost::program_options::variables_map& vm, const std::string& nameOption) {
	if (vm.count(nameOption)) {
		try {
			return vm[nameOption].as<T>();
		}
		catch (const boost::bad_any_cast& e) {
			std::cout << "Error, option " << nameOption << " wrongly specified! (bad type?)" << std::endl;
			exit(133);
		}
		catch (const Update::NotFoundOption& e) {
			std::cout << "Error, option " << nameOption << " not found!" << std::endl;
			exit(113);
		}
	}
	else {
		throw NotFoundOption(nameOption);
	}
}

template<> std::vector< std::complex<double> > get< std::vector< std::complex<double> > >(const boost::program_options::variables_map& vm, const std::string& nameOption);

template<> std::vector< double > get< std::vector< double > >(const boost::program_options::variables_map& vm, const std::string& nameOption);

template<> std::vector< unsigned int > get< std::vector< unsigned int > >(const boost::program_options::variables_map& vm, const std::string& nameOption);

template<> std::vector<std::string> get< std::vector<std::string> >(const boost::program_options::variables_map& vm, const std::string& nameOption);

}

class StorageParameters {
public:
	StorageParameters();
	StorageParameters(const boost::program_options::variables_map& _vm);
	~StorageParameters();

	/**
	 * This function returns back the wanted option with name "nameOption"
	 *  (usage:
	 *   StorageParameters sp(vm);
	 *   double val = sp.get<double>("my_value_name");
	 *  )
	 * @param nameOption
	 * @return
	 */
	template<typename T> T get(const std::string& nameOption) const {
		return implement::get<T>(vm, nameOption);
	}

private:
	boost::program_options::variables_map vm;
};

} /* namespace Update */
#endif /* STORAGEPARAMETERS_H_ */
