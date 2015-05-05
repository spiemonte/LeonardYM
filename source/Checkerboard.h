/*
 * Checkerbord.h
 *
 *  Created on: Apr 11, 2012
 *      Author: spiem_01
 */

#ifndef CHECKERBOARD_H_
#define CHECKERBOARD_H_
#include "Environment.h"
#include <list>

namespace Update {

struct Site {
	Site(int _site, unsigned int _mu) : site(site), mu(_mu) { }
	int site;
	int mu;
};

class Checkerboard {
	Checkerboard();
	void statistics();
	
	extended_index_lattice_t checkerboard;

	std::vector<Site>* siteColorList;

	static Checkerboard* ptr;
	int numberLoops;
	int* numberElementLoops;
public:
	~Checkerboard();

	static Checkerboard* getInstance();
	
	int getNumberLoops() const;
	int getColor(int site, int mu) const;

	std::vector<Site>* getSiteColorList() const;

private:
	void addValue(std::list<int>& l, int value) const;
	int getValue(std::list<int>& l) const;
	void updateEqual(std::list<int>& l1, std::list<int>& l2, int& value1, int& value2);
	void updateDifferent(std::list<int>& l1, std::list<int>& l2, int value1, int value2);

	int checkCorrectness() const;
};

} /* namespace Update */
#endif /* CHECKERBORD_H_ */
