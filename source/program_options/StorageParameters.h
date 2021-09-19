#ifndef STORAGEPARAMETERS_H_
#define STORAGEPARAMETERS_H_
#include <string>
#include <any>
#include <iostream>
#include <vector>
#include <complex>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <map>
#include "utils/FromString.h"
#include "program_options/Option.h"

namespace Update {

class NotFoundOption : public std::runtime_error {
public:
	NotFoundOption(const std::string& _nameOption) throw () : std::runtime_error((_nameOption+": option not found!").c_str()) { }
	~NotFoundOption() throw () { };
};

namespace implement {

template<typename T> T get(const std::map<std::string, Option>& vm, const std::string& nameOption) {
	if (vm.count(nameOption)) {
		try {
			return vm.at(nameOption).as<T>();
		}
		catch (...) {
			std::cout << "Error, option " << nameOption << " wrongly specified! (bad type?)" << std::endl;
			exit(133);
		}
	}
	else {
		throw NotFoundOption(nameOption);
	}
}

template<> std::vector< std::complex<double> > get< std::vector< std::complex<double> > >(const std::map<std::string, Option>& vm, const std::string& nameOption);

template<> std::vector< double > get< std::vector< double > >(const std::map<std::string, Option>& vm, const std::string& nameOption);

template<> std::vector< unsigned int > get< std::vector< unsigned int > >(const std::map<std::string, Option>& vm, const std::string& nameOption);

template<> std::vector<std::string> get< std::vector<std::string> >(const std::map<std::string, Option>& vm, const std::string& nameOption);

}

class StorageParameters {
public:
	StorageParameters();
	StorageParameters(const std::map<std::string, Option>& _vm);
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
	std::map<std::string, Option> vm;
};

} /* namespace Update */
#endif /* STORAGEPARAMETERS_H_ */
