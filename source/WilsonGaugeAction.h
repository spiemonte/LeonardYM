/*
 * WilsonGaugeAction.h
 *
 *  Created on: Feb 28, 2012
 *      Author: spiem_01
 */

#ifndef WILSONGAUGEACTION_H_
#define WILSONGAUGEACTION_H_

#include "GaugeAction.h"

namespace Update {

class WilsonGaugeAction : public Update::GaugeAction {
public:
	WilsonGaugeAction(real_t _beta);
	~WilsonGaugeAction();

#ifndef __IBMCPP__
	using Force::force;
#endif

	/**
	 * This function returns back the staple related to the link (site,mu) of the standard wilson action
	 * @param lattice
	 * @param site
	 * @param mu
	 * @return the staple of the link
	 */
	virtual GaugeGroup staple(const extended_gauge_lattice_t& lattice, int site, int mu) const;

	virtual GaugeGroup force(const extended_gauge_lattice_t& lattice, int site, int mu) const;

	virtual long_real_t energy(const environment_t& env);

	virtual real_t deltaAction(const extended_gauge_lattice_t& lattice, const GaugeGroup& trial, const GaugeGroup& staple, int site, int mu) const;
};

} /* namespace Update */
#endif /* WILSONGAUGEACTION_H_ */
