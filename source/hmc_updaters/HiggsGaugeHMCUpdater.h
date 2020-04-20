#ifndef HIGGSGAUGEHMCUPDATER_H_
#define HIGGSGAUGEHMCUPDATER_H_

#include "LatticeSweep.h"
#include "hmc_integrators/Integrate.h"
#include "HMCUpdater.h"
#include "actions/AdjointScalarAction.h"
#include "actions/FundamentalScalarAction.h"
#include "actions/MultiScalarAction.h"

namespace Update {

class HiggsGaugeHMCUpdater: public Update::LatticeSweep, public HMCUpdater {
public:
	HiggsGaugeHMCUpdater();
	HiggsGaugeHMCUpdater(const HiggsGaugeHMCUpdater& toCopy);
	~HiggsGaugeHMCUpdater();

	virtual void execute(environment_t& environment);

	static void registerParameters(po::options_description& desc);
protected:
        void initializeScalarAction(environment_t& environment);

private:
	//The new environment, provided by HMC
	environment_t environmentNew;
	//The conjugate momenta
	extended_gauge_lattice_t momenta;

	//The scalar action of the theory
        MultiScalarAction* scalarAction;
};

} /* namespace Update */
#endif /* PUREGAUGEHMCUPDATER_H_ */
