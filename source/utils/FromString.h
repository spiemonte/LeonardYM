#ifndef FROMSTRING_H_
#define FROMSTRING_H_
#include <sstream>

namespace Update {

template<typename T> T fromString(const char* toPrint) {
	std::istringstream iss(toPrint);
	T ris;
	iss >> ris;
	return ris;
}

template<typename T> T fromString(const std::string& toPrint) {
	std::istringstream iss(toPrint);
	T ris;
	iss >> ris;
	return ris;
}

}


#endif /* FROMSTRING_H_ */
