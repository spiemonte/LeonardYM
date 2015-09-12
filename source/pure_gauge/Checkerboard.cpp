/*
 * Checkerboard.cpp
 *
 *  Created on: Apr 11, 2012
 *      Author: spiem_01
 */

#include "Checkerboard.h"
#include <vector>
#include <algorithm>
#include <list>

namespace Update {

Checkerboard* Checkerboard::ptr = 0;

Checkerboard::Checkerboard() {
	typedef extended_gauge_lattice_t LT;
	typedef extended_gauge_lattice_t::Layout Layout;
	//Math::LinkConf::LinkConfiguration<unsigned int> test;
	std::list<int> (*checkerboardList)[4] = new std::list<int>[Layout::completesize][4];
	//Initalize the checkerboard
	for (int site = 0; site < Layout::localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			checkerboardList[site][mu].push_back(0);
			checkerboard[site][mu] = 0;
		}
	}
	//Now evaluate the checkerboard
	for (unsigned int j = 0; j < 3 && this->checkCorrectness() != 0; ++j) {
#ifdef ENABLE_MPI
		for (int processor = 0; processor < Layout::numberProcessors; ++processor) {
			if (processor == Layout::this_processor) {
#endif
				for (int site = 0; site < Layout::localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						for (unsigned int nu = 0; nu < 4; ++nu) {
							if (nu != mu) {
								//The links needed by the plaquette part of the action
								if (checkerboard[site][mu] == checkerboard[site][nu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[site][nu], checkerboard[site][mu], checkerboard[site][nu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[site][nu], checkerboard[site][mu], checkerboard[site][nu]);

								if (checkerboard[site][mu] == checkerboard[LT::sdn(site, mu)][mu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sdn(site, mu)][mu], checkerboard[site][mu], checkerboard[LT::sdn(site, mu)][mu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sdn(site, mu)][mu], checkerboard[site][mu], checkerboard[LT::sdn(site, mu)][mu]);

								if (checkerboard[site][mu] == checkerboard[LT::sdn(site, mu)][nu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sdn(site, mu)][nu], checkerboard[site][mu], checkerboard[LT::sdn(site, mu)][nu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sdn(site, mu)][nu], checkerboard[site][mu], checkerboard[LT::sdn(site, mu)][nu]);

								if (checkerboard[site][mu] == checkerboard[LT::sup(site, mu)][mu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sup(site, mu)][mu], checkerboard[site][mu], checkerboard[LT::sup(site, mu)][mu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sup(site, mu)][mu], checkerboard[site][mu], checkerboard[LT::sup(site, mu)][mu]);

								if (checkerboard[site][mu] == checkerboard[LT::sup(site, mu)][nu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sup(site, mu)][nu], checkerboard[site][mu], checkerboard[LT::sup(site, mu)][nu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sup(site, mu)][nu], checkerboard[site][mu], checkerboard[LT::sup(site, mu)][nu]);

								if (checkerboard[site][mu] == checkerboard[LT::sup(LT::sup(site, mu), mu)][nu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sup(LT::sup(site, mu), mu)][nu], checkerboard[site][mu], checkerboard[LT::sup(LT::sup(site, mu), mu)][nu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sup(LT::sup(site, mu), mu)][nu], checkerboard[site][mu], checkerboard[LT::sup(LT::sup(site, mu), mu)][nu]);

								if (checkerboard[site][mu] == checkerboard[LT::sdn(LT::sdn(site, nu), nu)][mu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sdn(LT::sdn(site, nu), nu)][mu], checkerboard[site][mu], checkerboard[LT::sdn(LT::sdn(site, nu), nu)][mu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sdn(LT::sdn(site, nu), nu)][mu], checkerboard[site][mu], checkerboard[LT::sdn(LT::sdn(site, nu), nu)][mu]);

								if (checkerboard[site][mu] == checkerboard[LT::sdn(LT::sdn(site, nu), nu)][nu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sdn(LT::sdn(site, nu), nu)][nu], checkerboard[site][mu], checkerboard[LT::sdn(LT::sdn(site, nu), nu)][nu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sdn(LT::sdn(site, nu), nu)][nu], checkerboard[site][mu], checkerboard[LT::sdn(LT::sdn(site, nu), nu)][nu]);

								if (checkerboard[site][mu] == checkerboard[LT::sup(LT::sdn(LT::sdn(site, nu), nu), mu)][nu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sup(LT::sdn(LT::sdn(site, nu), nu), mu)][nu], checkerboard[site][mu], checkerboard[LT::sup(LT::sdn(LT::sdn(site, nu), nu), mu)][nu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sup(LT::sdn(LT::sdn(site, nu), nu), mu)][nu], checkerboard[site][mu], checkerboard[LT::sup(LT::sdn(LT::sdn(site, nu), nu), mu)][nu]);

								if (checkerboard[site][mu] == checkerboard[LT::sdn(site, nu)][mu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sdn(site, nu)][mu], checkerboard[site][mu], checkerboard[LT::sdn(site, nu)][mu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sdn(site, nu)][mu], checkerboard[site][mu], checkerboard[LT::sdn(site, nu)][mu]);

								if (checkerboard[site][mu] == checkerboard[LT::sdn(site, nu)][nu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sdn(site, nu)][nu], checkerboard[site][mu], checkerboard[LT::sdn(site, nu)][nu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sdn(site, nu)][nu], checkerboard[site][mu], checkerboard[LT::sdn(site, nu)][nu]);

								if (checkerboard[site][mu] == checkerboard[LT::sdn(LT::sdn(site, mu), nu)][mu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sdn(LT::sdn(site, mu), nu)][mu], checkerboard[site][mu], checkerboard[LT::sdn(LT::sdn(site, mu), nu)][mu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sdn(LT::sdn(site, mu), nu)][mu], checkerboard[site][mu], checkerboard[LT::sdn(LT::sdn(site, mu), nu)][mu]);

								if (checkerboard[site][mu] == checkerboard[LT::sdn(LT::sdn(site, mu), nu)][nu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sdn(LT::sdn(site, mu), nu)][nu], checkerboard[site][mu], checkerboard[LT::sdn(LT::sdn(site, mu), nu)][nu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sdn(LT::sdn(site, mu), nu)][nu], checkerboard[site][mu], checkerboard[LT::sdn(LT::sdn(site, mu), nu)][nu]);

								if (checkerboard[site][mu] == checkerboard[LT::sup(LT::sdn(site, nu), mu)][mu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sup(LT::sdn(site, nu), mu)][mu], checkerboard[site][mu], checkerboard[LT::sup(LT::sdn(site, nu), mu)][mu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sup(LT::sdn(site, nu), mu)][mu], checkerboard[site][mu], checkerboard[LT::sup(LT::sdn(site, nu), mu)][mu]);

								if (checkerboard[site][mu] == checkerboard[LT::sup(LT::sdn(site, nu), mu)][nu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sup(LT::sdn(site, nu), mu)][nu], checkerboard[site][mu], checkerboard[LT::sup(LT::sdn(site, nu), mu)][nu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sup(LT::sdn(site, nu), mu)][nu], checkerboard[site][mu], checkerboard[LT::sup(LT::sdn(site, nu), mu)][nu]);

								if (checkerboard[site][mu] == checkerboard[LT::sdn(LT::sup(LT::sup(site, mu), mu), nu)][nu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sdn(LT::sup(LT::sup(site, mu), mu), nu)][nu], checkerboard[site][mu], checkerboard[LT::sdn(LT::sup(LT::sup(site, mu), mu), nu)][nu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sdn(LT::sup(LT::sup(site, mu), mu), nu)][nu], checkerboard[site][mu], checkerboard[LT::sdn(LT::sup(LT::sup(site, mu), mu), nu)][nu]);

								if (checkerboard[site][mu] == checkerboard[LT::sup(site, nu)][mu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sup(site, nu)][mu], checkerboard[site][mu], checkerboard[LT::sup(site, nu)][mu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sup(site, nu)][mu], checkerboard[site][mu], checkerboard[LT::sup(site, nu)][mu]);

								if (checkerboard[site][mu] == checkerboard[LT::sup(site, nu)][nu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sup(site, nu)][nu], checkerboard[site][mu], checkerboard[LT::sup(site, nu)][nu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sup(site, nu)][nu], checkerboard[site][mu], checkerboard[LT::sup(site, nu)][nu]);

								if (checkerboard[site][mu] == checkerboard[LT::sup(LT::sdn(site, mu), nu)][mu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sup(LT::sdn(site, mu), nu)][mu], checkerboard[site][mu], checkerboard[LT::sup(LT::sdn(site, mu), nu)][mu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sup(LT::sdn(site, mu), nu)][mu], checkerboard[site][mu], checkerboard[LT::sup(LT::sdn(site, mu), nu)][mu]);

								if (checkerboard[site][mu] == checkerboard[LT::sup(LT::sup(site, mu), nu)][mu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sup(LT::sup(site, mu), nu)][mu], checkerboard[site][mu], checkerboard[LT::sup(LT::sup(site, mu), nu)][mu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sup(LT::sup(site, mu), nu)][mu], checkerboard[site][mu], checkerboard[LT::sup(LT::sup(site, mu), nu)][mu]);

								if (checkerboard[site][mu] == checkerboard[LT::sup(LT::sup(site, mu), nu)][nu]) this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sup(LT::sup(site, mu), nu)][nu], checkerboard[site][mu], checkerboard[LT::sup(LT::sup(site, mu), nu)][nu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sup(LT::sup(site, mu), nu)][nu], checkerboard[site][mu], checkerboard[LT::sup(LT::sup(site, mu), nu)][nu]);

								if (checkerboard[site][mu] == checkerboard[LT::sup(LT::sup(site, nu), nu)][mu])	this->updateEqual(checkerboardList[site][mu],  checkerboardList[LT::sup(LT::sup(site, nu), nu)][mu], checkerboard[site][mu], checkerboard[LT::sup(LT::sup(site, nu), nu)][mu]);
								else this->updateDifferent(checkerboardList[site][mu],  checkerboardList[LT::sup(LT::sup(site, nu), nu)][mu], checkerboard[site][mu], checkerboard[LT::sup(LT::sup(site, nu), nu)][mu]);
							}
						}
					}
				}
#ifdef ENABLE_MPI
			}
			checkerboard.updateHalo();
		}
#endif
	}
	delete[] checkerboardList;
	this->statistics();

	int incorrect_sites = this->checkCorrectness();
	if (isOutputProcess() && incorrect_sites != 0) std::cout << "Checkerboard::WARNING, checkerboard incorrectly built, " << incorrect_sites << " cross dependent sites!" << std::endl;
}

Checkerboard::~Checkerboard() { }

void Checkerboard::statistics() {
	typedef extended_gauge_lattice_t::Layout Layout;
	numberLoops = 0;
	std::vector<int> datas;
	//std::map<>
	for (int site = 0; site < Layout::localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			std::vector<int>::iterator it = std::find(datas.begin(), datas.end(), checkerboard[site][mu]);
			if (it == datas.end()) {
				datas.push_back(checkerboard[site][mu]);
			}
		}
	}

#ifdef ENABLE_MPI
	//Now we take all the data from the other MPI nodes
	for (int processor = 0; processor < Layout::numberProcessors; ++processor) {
		int size_vector = -1;
		if (processor == Layout::this_processor) size_vector = datas.size();
		MPI_Bcast(&size_vector, 1, MPI_INT, processor, MPI_COMM_WORLD);
		int* tmp_vector = new int[size_vector];
		if (processor == Layout::this_processor) {
			for (int i = 0; i < size_vector; ++i) {
				tmp_vector[i] = datas[i];
			}
		}
		MPI_Bcast(tmp_vector, size_vector, MPI_INT, processor, MPI_COMM_WORLD);
		for (int i = 0; i < size_vector; ++i) {
			std::vector<int>::iterator it = std::find(datas.begin(), datas.end(), tmp_vector[i]);
			if (it == datas.end()) {
				datas.push_back(tmp_vector[i]);
			}
		}
		
	}
#endif
	
	if (isOutputProcess()) {
		std::cout << "Datas: {";
		for (unsigned int i = 0; i < datas.size() - 1; ++i) {
			std::cout << datas[i] << ", ";
		}
		std::cout << datas.back() << "}" << std::endl;
	}

	numberLoops = datas.size();

#ifdef ENABLE_MPI
	for (int processor = 0; processor < Layout::numberProcessors; ++processor) {
		if (processor == Layout::this_processor) {
			std::cout << "Number of loops for the mth environment for the processor "<< processor <<": " << numberLoops << std::endl;
#endif
#ifndef ENABLE_MPI
			std::cout << "Number of loops for the mth environment: " << numberLoops << std::endl;
#endif
			numberElementLoops = new int[numberLoops];
			for (int i = 0; i < numberLoops; ++i) {
				numberElementLoops[i] = 0;
			}
			siteColorList = new std::vector<Site>[numberLoops];
			for (int site = 0; site < Layout::localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					++numberElementLoops[checkerboard[site][mu]];
					siteColorList[checkerboard[site][mu]].push_back(Site(site, mu));
				}
			}

			std::cout << "Links in the loops:  {";
			for (int i = 0; i < numberLoops - 1; ++i) {
				std::cout << numberElementLoops[i] << ", ";
			}
			std::cout << numberElementLoops[numberLoops - 1] << "}" << std::endl;
#ifdef ENABLE_MPI
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif
}

Checkerboard* Checkerboard::getInstance() {
	if (ptr == 0) {
		ptr = new Checkerboard();
		return ptr;
	}
	else return ptr;
}

int Checkerboard::getNumberLoops() const {
	return numberLoops;
}

std::vector<Site>* Checkerboard::getSiteColorList() const {
	return siteColorList;
}

int Checkerboard::getColor(int site, int mu) const {
	return checkerboard[site][mu];
}

void Checkerboard::addValue(std::list<int>& l, int value) const {
	if (l.empty()) {
		l.push_back(value);
	} else {
		std::list<int>::iterator i;
		for (i = l.begin(); i != l.end(); ++i) {
			if (*i == value) return;
			else if (*i > value) {
				l.insert(i, value);
				return;
			}
		}
		l.push_back(value);
	}
}

int Checkerboard::getValue(std::list<int>& l) const {
	int tmp = 0;
	std::list<int>::iterator i;
	for (i = l.begin(); i != l.end(); ++i) {
		if (tmp != *i) return tmp;
		else ++tmp;
	}
	return tmp;
}

void Checkerboard::updateEqual(std::list<int>& l1, std::list<int>& l2, int& value1, int& value2) {
	/*//The second site cannot now take the value of the site 1
	this->addValue(l2, value1);
	//Get a free value from the list
	value2 = this->getValue(l2);
	//Now the first site cannot take the new value of the site 2
	this->addValue(l1, value2);*/
	//The first site cannot now take the value of the site 2
	this->addValue(l1, value2);
	//Get a free value from the list
	value1 = this->getValue(l1);
	//Now the second site cannot take the new value of the site 1
	this->addValue(l2, value1);
}

void Checkerboard::updateDifferent(std::list<int>& l1, std::list<int>& l2, int value1, int value2) {
	//The second site cannot now take the value of the site 1
	this->addValue(l2, value1);
	//The first site cannot now take the value of the site 2
	this->addValue(l1, value2);
}

int Checkerboard::checkCorrectness() const {
	typedef extended_gauge_lattice_t LT;
	typedef extended_gauge_lattice_t::Layout Layout;
	int incorrect = 0;
	
#pragma omp parallel for reduction(+:incorrect)
	for (int site = 0; site < Layout::localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = 0; nu < 4; ++nu) {
				if (nu != mu) {
					//The links needed by the plaquette part of the action
					if (checkerboard[site][mu] == checkerboard[site][nu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sdn(site, mu)][mu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sdn(site, mu)][nu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sup(site, mu)][mu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sup(site, mu)][nu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sup(LT::sup(site, mu), mu)][nu])	incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sdn(LT::sdn(site, nu), nu)][mu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sdn(LT::sdn(site, nu), nu)][nu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sup(LT::sdn(LT::sdn(site, nu), nu), mu)][nu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sdn(site, nu)][mu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sdn(site, nu)][nu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sdn(LT::sdn(site, mu), nu)][mu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sdn(LT::sdn(site, mu), nu)][nu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sup(LT::sdn(site, nu), mu)][mu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sup(LT::sdn(site, nu), mu)][nu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sdn(LT::sup(LT::sup(site, mu), mu), nu)][nu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sup(site, nu)][mu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sup(site, nu)][nu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sup(LT::sdn(site, mu), nu)][mu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sup(LT::sup(site, mu), nu)][mu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sup(LT::sup(site, mu), nu)][nu]) incorrect += 1;
					if (checkerboard[site][mu] == checkerboard[LT::sup(LT::sup(site, nu), nu)][mu]) incorrect += 1;
				}
			}
		}
	}

	reduceAllSum(incorrect);
	return incorrect;
}

} /* namespace Update */
