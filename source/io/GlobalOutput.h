#ifndef GLOBALOUTPUT_H_
#define GLOBALOUTPUT_H_
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <map>
#include <iomanip>

namespace Update {

inline bool file_exists(const char *fileName) {
	std::ifstream infile(fileName);
	return infile.good();
}

//Singleton class to manage to output of the measurements
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
		if (format == "txt") {
			if (it != output_streams.end()) {
				output_status[name] = false;
				it->second.append("\t{");
			}
			else {
				output_streams[name] = "{\n\t{";
				output_status[name] = false;
			}
		}
		else if (format == "xml")  {
			if (it != output_streams.end()) {
				output_status[name] = false;
				it->second.append(std::string("<")+name+">");
			}
			else {
				output_streams[name] = std::string("<")+name+">";
				output_status[name] = false;
			}
		}
	}

	void pop(const std::string& name) {
		std::map<std::string, std::string>::iterator it = output_streams.find(name);
		if (format == "txt") {
			if (it != output_streams.end()) {
				output_status[name] = false;
				it->second.append("},\n");
			}
			else {
				std::cout << "Fatal error in output, probably pop called without push" << std::endl;
				exit(1);
			}
		}
		else if (format == "xml") {
			if (it != output_streams.end()) {
				output_status[name] = false;
				it->second.append(std::string("</")+name+">");
			}
			else {
				std::cout << "Fatal error in output, probably pop called without push" << std::endl;
				exit(1);
			}
		}
	}

	template<typename T> void write(const std::string& name, const T& what) {
		std::map<std::string, std::string>::iterator it = output_streams.find(name);
		if (format == "txt") {
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
		else {
			if (it != output_streams.end()) {
				if (output_status[name] == true) {
					it->second.append(" ");
				}
				else {
					output_status[name] = true;
				}
				std::ostringstream oss;
				oss << std::setprecision(23) << what << " ";
				it->second.append(oss.str());
			}
			else {
				std::ostringstream oss;
				oss << std::setprecision(23) << what << " ";
				output_streams[name] = oss.str();
				output_status[name] = true;
			}
		}
	}

	void print() {
		std::map<std::string, std::string>::iterator it;
		for (it=output_streams.begin(); it != output_streams.end(); ++it) {
			if (format == "txt") {
				std::ofstream ofs;
				ofs.open((baseFolder+baseName+"_"+it->first+".txt").c_str(), std::fstream::out | std::fstream::app);
				ofs << it->second;
				it->second.clear();
				ofs.close();
			}
			else if (format == "xml") {
				if (file_exists((baseFolder+baseName+"_"+it->first+".xml").c_str())) {
					std::ifstream ifs;
					ifs.open((baseFolder+baseName+"_"+it->first+".xml").c_str(), std::fstream::in);
					std::stringstream buffer;
					buffer << ifs.rdbuf();
					std::string contents = buffer.str();
					contents.erase(contents.end()-9,contents.end());
					ifs.close();


					std::ofstream ofs;
					ofs.open((baseFolder+baseName+"_"+it->first+".xml").c_str(), std::fstream::trunc);
					ofs << contents;	
					ofs << "\n\n" << it->second;
					ofs << "\n</root>\n";
					it->second.clear();
					ofs.close();
				}
				else {
					std::ofstream ofs;
					ofs.open((baseFolder+baseName+"_"+it->first+".xml").c_str(), std::fstream::out | std::fstream::app);
					ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\n<root>\n";
					ofs << it->second;
					ofs << "\n</root>\n";
					it->second.clear();
					ofs.close();
				}
			}
		}
	}

	static void setBaseName(const std::string& _baseName) {
		baseName = _baseName;
	}

	static void setBaseFolder(const std::string& _baseFolder) {
		baseFolder = _baseFolder;
	}

	static void setFormat(const std::string& _format) {
		if (_format == "xml") {
			format = "xml";
		}
		else if (_format == "txt") {
			format = "txt";
		}
		else {
			std::cout << "Ouput format " << _format << " unsupported!" << std::endl;
		}
	}
private:
	std::map<std::string, std::string> output_streams;
	std::map<std::string, bool> output_status;
	static std::string baseName;
	static std::string baseFolder;
	static std::string format;
};

} /* namespace Update */
#endif /* GLOBALOUTPUT_H_ */
