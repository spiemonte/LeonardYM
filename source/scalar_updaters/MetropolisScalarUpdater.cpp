#include "MetropolisScalarUpdater.h"
#include "utils/RandomSeed.h"
#include "actions/ScalarAction.h"
#include <omp.h>

namespace Update {

    MetropolisScalarUpdater::MetropolisScalarUpdater()
#ifndef MULTITHREADING
    : LatticeSweep(), randomGenerator(RandomSeed::randomSeed()), randomNormal(RandomSeed::getNormalNumberGenerator(randomGenerator,sqrt(0.5))), randomUniform(RandomSeed::getRandomNumberGenerator(randomGenerator)) { }
#endif
#ifdef MULTITHREADING
    : LatticeSweep()
    {
        randomGenerator = new random_generator_t*[omp_get_max_threads()];
        randomNormal = new random_normal_generator_t*[omp_get_max_threads()];
        randomUniform = new random_uniform_generator_t*[omp_get_max_threads()];
        for (int i = 0; i < omp_get_max_threads(); ++i) {
            randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
            randomNormal[i] = new random_normal_generator_t(RandomSeed::getNormalNumberGenerator(*randomGenerator[i],sqrt(0.5)));
            randomUniform[i] = new random_uniform_generator_t(RandomSeed::getRandomNumberGenerator(*randomGenerator[i]));
        }
    }
#endif

    MetropolisScalarUpdater::MetropolisScalarUpdater(const MetropolisScalarUpdater& toCopy)
#ifndef MULTITHREADING
    : LatticeSweep(toCopy), randomGenerator(RandomSeed::randomSeed()), randomNormal(RandomSeed::getNormalNumberGenerator(randomGenerator,sqrt(0.5))), randomUniform(RandomSeed::getRandomNumberGenerator(randomGenerator)) { }
#endif
#ifdef MULTITHREADING
    : LatticeSweep(toCopy) 
    {
        randomGenerator = new random_generator_t*[omp_get_max_threads()];
        randomNormal = new random_normal_generator_t*[omp_get_max_threads()];
        randomUniform = new random_uniform_generator_t*[omp_get_max_threads()];
        for (int i = 0; i < omp_get_max_threads(); ++i) {
            randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
            randomNormal[i] = new random_normal_generator_t(RandomSeed::getNormalNumberGenerator(*randomGenerator[i],sqrt(0.5)));
            randomUniform[i] = new random_uniform_generator_t(RandomSeed::getRandomNumberGenerator(*randomGenerator[i]));
        }
    }
#endif

    MetropolisScalarUpdater::~MetropolisScalarUpdater() {
#ifdef MULTITHREADING
        for (int i = 0; i < omp_get_max_threads(); ++i) {
            delete randomGenerator[i];
            delete randomNormal[i];
            delete randomUniform[i];
        }
        delete[] randomGenerator;
        delete[] randomNormal;
        delete[] randomUniform;
#endif
    }

    void MetropolisScalarUpdater::execute(environment_t& environment) {
        unsigned int aNf = environment.configurations.get<unsigned int>("adjoint_nf_scalars");
        if (aNf != environment.adjoint_scalar_fields.size()) {
            if (isOutputProcess()) std::cout << "MetropolisScalarUpdater::Number of adjoint scalars incorrectly set!" << std::endl;
            exit(109);
        } 

        unsigned int fNf = environment.configurations.get<unsigned int>("fundamental_nf_scalars");
        if (fNf != environment.fundamental_scalar_fields.size()) {
            if (isOutputProcess()) std::cout << "MetropolisScalarUpdater::Number of fundamental scalars incorrectly set!" << std::endl;
            exit(107);
        }

        double adjoint_scalar_mass = environment.configurations.get<double>("MetropolisScalarUpdater::adjoint_scalar_mass");
        double fundamental_scalar_mass = environment.configurations.get<double>("MetropolisScalarUpdater::fundamental_scalar_mass");
        double lambda_adjoint = environment.configurations.get<double>("MetropolisScalarUpdater::adjoint_quartic_coupling");
        double lambda_fundamental = environment.configurations.get<double>("MetropolisScalarUpdater::fundamental_quartic_coupling");
        double lambda_8 = environment.configurations.get<double>("MetropolisScalarUpdater::lambda_8");
        double lambda_mixed = environment.configurations.get<double>("MetropolisScalarUpdater::lambda_mixed");

        double epsilon = environment.configurations.get<double>("MetropolisScalarUpdater::epsilon");

        ScalarAction* action = new ScalarAction(adjoint_scalar_mass, fundamental_scalar_mass, lambda_fundamental, lambda_adjoint, lambda_mixed, lambda_8);

        unsigned int trials = environment.configurations.get<unsigned int>("MetropolisScalarUpdater::number_of_hits");

        typedef extended_adjoint_real_color_vector_t::Layout Layout;

        int acceptance = 0;

    /*{
        AdjointVector oldphi = environment.adjoint_scalar_fields[0][5];
        AdjointVector newphi = environment.adjoint_scalar_fields[0][5];
        newphi[0] += 0.22;
        newphi[1] -= 0.13;
        AdjointVector kinetic_coupling = action->getKineticCoupling(environment.getAdjointLattice(), environment.adjoint_scalar_fields[0], 5);
        real_t deltaMet = action->deltaEnergy(kinetic_coupling, oldphi, newphi);
        real_t deltaE = -action->energy(environment);
        environment.adjoint_scalar_fields[0][5] = newphi;
        deltaE += action->energy(environment);
        if (isOutputProcess()) std::cout << "Consistency test of the action " << deltaMet << " " << deltaE << std::endl; 
    }*/

        for (int block = 0; block < 2; ++block) {
#pragma omp parallel for reduction(+:acceptance)
            for (int site = 0; site < environment.getAdjointLattice().localsize; ++site) {
                //White/Black partitioning
                if ((Layout::globalIndexX(site) + Layout::globalIndexY(site) + Layout::globalIndexZ(site) + Layout::globalIndexT(site)) % 2 == block) {
                    for (unsigned int field_index = 0; field_index < aNf; ++field_index) {
                        AdjointRealVector adjoint_kinetic_coupling = action->getKineticCoupling(environment.getAdjointLattice(), environment.adjoint_scalar_fields[field_index], site);

                        std::vector<AdjointRealVector*> adjoint_scalar_fields_on_site = action->scalarFieldsAtSite(environment.adjoint_scalar_fields, site);
                        std::vector<FundamentalVector*> fundamental_scalar_fields_on_site = action->scalarFieldsAtSite(environment.fundamental_scalar_fields, site);

                        for (unsigned int trial = 0; trial < trials; ++trial) {
                            for (unsigned int c = 0; c < numberColors*numberColors - 1; ++c) {
                            
#ifndef MULTITHREADING
                                real_t proposal = epsilon*randomNormal();
#endif
#ifdef MULTITHREADING
                                real_t proposal = epsilon*(*randomNormal[omp_get_thread_num()])();
#endif
                                real_t oldAction = action->energy(adjoint_scalar_fields_on_site, fundamental_scalar_fields_on_site, adjoint_kinetic_coupling, field_index);
                                
                                (*adjoint_scalar_fields_on_site[field_index])[c] += proposal;

                                real_t newAction = action->energy(adjoint_scalar_fields_on_site, fundamental_scalar_fields_on_site, adjoint_kinetic_coupling, field_index);
                                
                                //Do the accept/reject metropolis
                                if (newAction - oldAction < 0.) {
                                    ++acceptance;
                                }
#ifndef MULTITHREADING
                                else if (randomUniform() < exp(-(newAction - oldAction))) {
#endif
#ifdef MULTITHREADING
                                else if ((*randomUniform[omp_get_thread_num()])() < exp(-(newAction - oldAction))) {
#endif
                                    ++acceptance;
                                }
                                else {
                                    //Proposal rejected, we restore the old field
                                    (*adjoint_scalar_fields_on_site[field_index])[c] -= proposal;
                                }
                            }
                        }

                    }

                    for (unsigned int field_index = 0; field_index < fNf; ++field_index) {
                        FundamentalVector fundamental_kinetic_coupling = action->getKineticCoupling(environment.getFundamentalLattice(), environment.fundamental_scalar_fields[field_index], site);

                        std::vector<AdjointRealVector*> adjoint_scalar_fields_on_site = action->scalarFieldsAtSite(environment.adjoint_scalar_fields, site);
                        std::vector<FundamentalVector*> fundamental_scalar_fields_on_site = action->scalarFieldsAtSite(environment.fundamental_scalar_fields, site);

                        for (unsigned int trial = 0; trial < trials; ++trial) {
                            for (unsigned int c = 0; c < numberColors; ++c) {
                            
#ifndef MULTITHREADING
                                std::complex<real_t> proposal = epsilon*std::complex<real_t>(randomNormal(),randomNormal());
#endif
#ifdef MULTITHREADING
                                std::complex<real_t> proposal = epsilon*std::complex<real_t>((*randomNormal[omp_get_thread_num()])(),(*randomNormal[omp_get_thread_num()])());
#endif
                                real_t oldAction = action->energy(adjoint_scalar_fields_on_site, fundamental_scalar_fields_on_site, fundamental_kinetic_coupling, field_index);
                                
                                (*fundamental_scalar_fields_on_site[field_index])[c] += proposal;

                                real_t newAction = action->energy(adjoint_scalar_fields_on_site, fundamental_scalar_fields_on_site, fundamental_kinetic_coupling, field_index);
                                
                                //Do the accept/reject metropolis
                                if (newAction - oldAction < 0.) {
                                    ++acceptance;
                                }
#ifndef MULTITHREADING
                                else if (randomUniform() < exp(-(newAction - oldAction))) {
#endif
#ifdef MULTITHREADING
                                else if ((*randomUniform[omp_get_thread_num()])() < exp(-(newAction - oldAction))) {
#endif
                                    ++acceptance;
                                }
                                else {
                                    //Proposal rejected, we restore the old field
                                    (*fundamental_scalar_fields_on_site[field_index])[c] -= proposal;
                                }
                            }
                        }

                    }
                }
            }
            for (unsigned int field_index = 0; field_index < aNf; ++field_index) {
                environment.adjoint_scalar_fields[field_index].updateHalo();
            }
            for (unsigned int field_index = 0; field_index < fNf; ++field_index) {
                environment.fundamental_scalar_fields[field_index].updateHalo();
            }
        }

        reduceAllSum(acceptance);

        if (isOutputProcess()) std::cout << "MetropolisScalarUpdater::Acceptance rate: " << acceptance/static_cast<double>(Layout::globalVolume*trials*(aNf*(numberColors*numberColors-1)+fNf*numberColors)) << std::endl;

        delete action;
    }

    void MetropolisScalarUpdater::registerParameters(std::map<std::string, Option>& desc) {
        desc["MetropolisScalarUpdater::adjoint_scalar_mass"] = Option("MetropolisScalarUpdater::adjoint_scalar_mass", 0.0, "set the value of the adjoint scalar mass m*m");
        desc["MetropolisScalarUpdater::fundamental_scalar_mass"] = Option("MetropolisScalarUpdater::fundamental_scalar_mass", 0.0, "set the value of the fundamental scalar mass m*m");

        desc["MetropolisScalarUpdater::adjoint_quartic_coupling"] = Option("MetropolisScalarUpdater::adjoint_quartic_coupling", 0.0, "set the value of the adjoint quartic coupling lambda");
        desc["MetropolisScalarUpdater::fundamental_quartic_coupling"] = Option("MetropolisScalarUpdater::fundamental_quartic_coupling", 0.0, "set the value of the fundamental quartic coupling lambda");

        desc["MetropolisScalarUpdater::lambda_8"] = Option("MetropolisScalarUpdater::lambda_8", 0.0, "set the value of the adjoint quartic coupling lambda");
        desc["MetropolisScalarUpdater::lambda_mixed"] = Option("MetropolisScalarUpdater::lambda_mixed", 0.0, "set the value of the mixed quartic coupling lambda");

        desc["MetropolisScalarUpdater::epsilon"] = Option("MetropolisScalarUpdater::epsilon", 0.1, "set the value of epsilon for the random update");
        desc["MetropolisScalarUpdater::number_of_hits"] = Option("MetropolisScalarUpdater::number_of_hits", 3, "set the number of the trials for the multihit metropolis");
    }

} /* namespace Update */
