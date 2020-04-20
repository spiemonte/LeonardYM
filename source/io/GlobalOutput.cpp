#include "GlobalOutput.h"

namespace Update {

GlobalOutput* GlobalOutput::ptr = 0;
std::string GlobalOutput::baseName = "";
std::string GlobalOutput::baseFolder = "";
std::string GlobalOutput::format = "";

GlobalOutput::GlobalOutput() { }

GlobalOutput::~GlobalOutput() { }

} /* namespace Update */
