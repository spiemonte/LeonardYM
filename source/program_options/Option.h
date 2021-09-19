#ifndef OPTION_H
#define OPTION_H
#include <any>
#include <string>
#include "utils/ToString.h"
#include "utils/FromString.h"

namespace Update {

class Option {
public:
	template<typename T> Option(const std::string& _name, const T& _default_value, const std::string& _help) : name(_name), help(_help) {
		value = Update::toString(_default_value);
	}

	Option(const std::string& _name) : name(_name) { }

	Option(const std::string& _name, const std::string& _help) : name(_name), help(_help) { }

	Option() {}

	template<typename T> T as() const {
		return fromString<T>(value);
	}

	Option& operator=(const Option& copy) {
		value = copy.value;
		name = copy.name;
		help = copy.help;
		return *this;
	}

	Option& operator=(const std::string& _value) {
		value = _value;
		return *this;
	}

	Option(const Option& copy) : value(copy.value), name(copy.name), help(copy.help) { }

	std::string getHelp() const {
		return name + ": " + help;
	}

	std::string toString() const {
		return name + " = " + value + " # " + help;
	}
private:
	std::string value;
	std::string name;
	std::string help;
};

}

#endif