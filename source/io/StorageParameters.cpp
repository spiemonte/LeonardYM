/*
 * StorageParameters.cpp
 *
 *  Created on: Feb 21, 2012
 *      Author: spiem_01
 */

#include "StorageParameters.h"

namespace Update {

StorageParameters::StorageParameters() { }

StorageParameters::StorageParameters(const boost::program_options::variables_map& _vm) : vm(_vm) { }

StorageParameters::~StorageParameters() { }

namespace implement {

template<> std::vector<std::string> get< std::vector<std::string> >(const boost::program_options::variables_map& vm, const std::string& nameOption) {
	if (vm.count(nameOption)) {
		std::string toRead = vm[nameOption].as<std::string>();
		std::vector<std::string> ris;
		std::string::const_iterator index = toRead.begin();
		while (index != toRead.end() && *index != '{') { ++index; }
		if (index == toRead.end()) {
			std::cout << "Empty or non-valid string given to parser of vector line" << std::endl;
			exit(1);
		}
		else {
			++index;
			unsigned int state = 1;
			std::string tmp;
			while (state != 0) {
				switch (*index) {
				case '{':
					++state;
					tmp.append(1,*index);
					++index;
					break;
				case '}':
					--state;
					if (state != 0) {
						tmp.append(1,*index);
						++index;
					}
					else {
						ris.push_back(tmp);
						tmp.clear();
						++index;
					}
					break;
				case ',':
					if (state != 1) {
						tmp.append(1,*index);
						++index;
					}
					else {
						ris.push_back(tmp);
						tmp.clear();
						++index;
					}
					break;
				case '(':
					++state;
					tmp.append(1,*index);
					++index;
					break;
				case ')':
					--state;
					tmp.append(1,*index);
					++index;
					break;
				default:
					tmp.append(1,*index);
					++index;
					break;
				}
				if (index == toRead.end()) {
					return ris;
				}
			}
		}
		return ris;
	}
	else {
		throw NotFoundOption(nameOption);
	}
}

template<> std::vector< std::complex<double> > get< std::vector< std::complex<double> > >(const boost::program_options::variables_map& vm, const std::string& nameOption) {
	std::vector<std::string> tmp = get< std::vector<std::string> >(vm, nameOption);
	std::vector< std::complex<double> > ris;
	std::vector<std::string>::iterator i;
	for (i = tmp.begin(); i != tmp.end(); ++i) {
		ris.push_back(fromString< std::complex<double> >(*i));
	}
	return ris;
}

template<> std::vector< double > get< std::vector< double > >(const boost::program_options::variables_map& vm, const std::string& nameOption) {
	std::vector<std::string> tmp = get< std::vector<std::string> >(vm, nameOption);
	std::vector< double > ris;
	std::vector<std::string>::iterator i;
	for (i = tmp.begin(); i != tmp.end(); ++i) {
		ris.push_back(fromString< double >(*i));
	}
	return ris;
}

template<> std::vector< unsigned int > get< std::vector< unsigned int > >(const boost::program_options::variables_map& vm, const std::string& nameOption) {
	std::vector<std::string> tmp = get< std::vector<std::string> >(vm, nameOption);
	std::vector< unsigned int > ris;
	std::vector<std::string>::iterator i;
	for (i = tmp.begin(); i != tmp.end(); ++i) {
		ris.push_back(fromString< unsigned int >(*i));
	}
	return ris;
}

} /* namespace implement */

} /* namespace Update */
