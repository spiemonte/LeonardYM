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
	
	unsigned int (*checkerboard)[4];

	std::vector<Site>* siteColorList;

	static Checkerboard* ptr;
	unsigned int numberLoops;
	int* numberElementLoops;
public:
	~Checkerboard();

	static Checkerboard* getInstance();
	
	unsigned int getNumberLoops() const;
	unsigned int getColor(int site, int mu) const;

	std::vector<Site>* getSiteColorList() const;

private:
	void addValue(std::list<unsigned int>& l, unsigned int value) const;
	unsigned int getValue(std::list<unsigned int>& l) const;
	void updateEqual(std::list<unsigned int>& l1, std::list<unsigned int>& l2, unsigned int value1, unsigned int& value2);
	void updateDifferent(std::list<unsigned int>& l1, std::list<unsigned int>& l2, unsigned int value1, unsigned int value2);
};

} /* namespace Update */
#endif /* CHECKERBORD_H_ */
