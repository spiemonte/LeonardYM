/*
 * StorageParameters.h
 *
 *  Created on: Feb 21, 2012
 *      Author: spiem_01
 */

#ifndef STORAGEPARAMETERS_H_
#define STORAGEPARAMETERS_H_
#include <string>
#include <boost/program_options.hpp>
#include <iostream>
#include <vector>
#include <complex>
#include <sstream>
#include <exception>
#include "utils/FromString.h"

namespace Update {

class NotFoundOption : public std::exception {
public:
	NotFoundOption(const std::string& _nameOption) throw () : nameOption(_nameOption) { }
	~NotFoundOption() throw () { };
private:
	virtual const char* what() const throw() {
		std::ostringstream os;
		os << "No parameter " << nameOption << " found!" << std::endl;
		return os.str().c_str();
	}
	std::string nameOption;
};

namespace implement {

template<typename T> T get(const boost::program_options::variables_map& vm, const std::string& nameOption) throw(NotFoundOption) {
	if (vm.count(nameOption)) {
		return vm[nameOption].as<T>();
	}
	else {
		throw NotFoundOption(nameOption);
	}
}

template<> std::vector< std::complex<double> > get< std::vector< std::complex<double> > >(const boost::program_options::variables_map& vm, const std::string& nameOption) throw(NotFoundOption);

template<> std::vector< double > get< std::vector< double > >(const boost::program_options::variables_map& vm, const std::string& nameOption) throw(NotFoundOption);

template<> std::vector< unsigned int > get< std::vector< unsigned int > >(const boost::program_options::variables_map& vm, const std::string& nameOption) throw(NotFoundOption);

template<> std::vector<std::string> get< std::vector<std::string> >(const boost::program_options::variables_map& vm, const std::string& nameOption) throw(NotFoundOption);

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

	//template<> std::vector<std::string> get< std::vector<std::string> >(const std::string& nameOption) const;

private:
	boost::program_options::variables_map vm;
};

} /* namespace Update */
#endif /* STORAGEPARAMETERS_H_ */
