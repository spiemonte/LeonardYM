#ifndef IMPROVEDGAUGEACTION_H_
#define IMPROVEDGAUGEACTION_H_

#include "GaugeAction.h"

namespace Update {

class ImprovedGaugeAction : public GaugeAction {
public:
	ImprovedGaugeAction(real_t _beta, real_t _u0 = 1.);
	~ImprovedGaugeAction();

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
private:
	real_t u0;
};

} /* namespace Update */
#endif /* IMPROVEDGAUGEACTION_H_ */
