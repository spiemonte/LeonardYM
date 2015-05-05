/*
 * GlobalOutput.h
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

#ifndef GLOBALOUTPUT_H_
#define GLOBALOUTPUT_H_
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <map>
#include <iomanip>

namespace Update {

class GlobalOutput {
	GlobalOutput();

	static GlobalOutput* ptr;
public:
	~GlobalOutput();

	static GlobalOutput* getInstance() {
		if (ptr == 0) {
			ptr = new GlobalOutput();
		}
		return ptr;
	}

	void destroy() {
		delete ptr;
		ptr = 0;
	}

	void push(const std::string& name) {
		std::map<std::string, std::string>::iterator it = output_streams.find(name);
		if (it != output_streams.end()) {
			output_status[name] = false;
			it->second.append("\t{");
		}
		else {
			output_streams[name] = "{\n\t{";
			output_status[name] = false;
		}
	}

	void pop(const std::string& name) {
		std::map<std::string, std::string>::iterator it = output_streams.find(name);
		if (it != output_streams.end()) {
			output_status[name] = false;
			it->second.append("},\n");
		}
		else {
			std::cout << "Fatal error in output" << std::endl;
			exit(1);
		}
	}

	template<typename T> void write(const std::string& name, const T& what) {
		std::map<std::string, std::string>::iterator it = output_streams.find(name);
		if (it != output_streams.end()) {
			if (output_status[name] == true) {
				it->second.append(", ");
			}
			else {
				output_status[name] = true;
			}
			std::ostringstream oss;
			oss << std::setprecision(23) << what;
			it->second.append(oss.str());
		}
		else {
			std::ostringstream oss;
			oss << "{" << std::setprecision(23) << what;
			output_streams[name] = oss.str();
			output_status[name] = true;
		}
	}

	void print() {
		std::map<std::string, std::string>::iterator it;
		for (it=output_streams.begin(); it != output_streams.end(); ++it) {
			std::ofstream ofs;
			ofs.open((baseFolder+baseName+"_"+it->first+".txt").c_str(), std::fstream::out | std::fstream::app);
			ofs << it->second;
			it->second.clear();
			ofs.close();
		}
	}

	static void setBaseName(const std::string& _baseName) {
		baseName = _baseName;
	}

	static void setBaseFolder(const std::string& _baseFolder) {
		baseFolder = _baseFolder;
	}
private:
	std::map<std::string, std::string> output_streams;
	std::map<std::string, bool> output_status;
	static std::string baseName;
	static std::string baseFolder;
};

} /* namespace Update */
#endif /* GLOBALOUTPUT_H_ */
