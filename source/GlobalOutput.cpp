/*
 * GlobalOutput.cpp
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

#include "GlobalOutput.h"

namespace Update {

GlobalOutput* GlobalOutput::ptr = 0;
std::string GlobalOutput::baseName = "";
std::string GlobalOutput::baseFolder = "";

GlobalOutput::GlobalOutput() { }

GlobalOutput::~GlobalOutput() { }

} /* namespace Update */
