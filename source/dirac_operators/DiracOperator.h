#ifndef DIRACOPERATOR_H_
#define DIRACOPERATOR_H_

#include "Environment.h"
#include "hmc_forces/FermionForce.h"

#include <string>

namespace Update {

class DiracOperator {
public:
	DiracOperator();
	DiracOperator(const extended_fermion_lattice_t& _lattice, double _kappa = 0., bool _gamma5 = true);
	virtual ~DiracOperator();
	static DiracOperator* getInstance(const std::string& name, unsigned int power, const StorageParameters& parameters, const std::string& basename = "");
	static DiracOperator* getSquareRoot(DiracOperator* dirac);
	static DiracOperator* getSquare(DiracOperator* dirac);

	static void registerParameters(std::map<std::string, Option>& desc, const std::string& basename);
	
#ifdef ENABLE_MPI
	void multiply(standard_dirac_vector_t& output, const standard_dirac_vector_t& input) {
		reduced_dirac_vector_t _output, _input = input;
		this->multiply(_output,_input);
		output = _output;
	}

	void multiply(extended_dirac_vector_t& output, const extended_dirac_vector_t& input) {
		reduced_dirac_vector_t _output, _input = input;
		this->multiply(_output,_input);
		output = _output;
	}

	virtual void multiplyAdd(standard_dirac_vector_t& output, const standard_dirac_vector_t& vector1, const standard_dirac_vector_t & vector2, const complex& alpha) {
		reduced_dirac_vector_t _output, _vector1 = vector1, _vector2 = vector2;
                this->multiplyAdd(_output, _vector1, _vector2, alpha);
                output = _output;
	}

	virtual void multiplyAdd(extended_dirac_vector_t& output, const extended_dirac_vector_t& vector1, const extended_dirac_vector_t & vector2, const complex& alpha) {
		reduced_dirac_vector_t _output, _vector1 = vector1, _vector2 = vector2;
		//if (isOutputProcess()) std::cout << "Vai a sapere" << std::endl;
		this->multiplyAdd(_output, _vector1, _vector2, alpha);
		//if (isOutputProcess()) std::cout << "Vai a sapere" << std::endl;
		output = _output;
        }
#endif

#ifdef ALIGNED_OPT
	//Dummy function to be overridden in the most relevant time-consuming cases
	virtual void multiply(reduced_soa_dirac_vector_t& output, const reduced_soa_dirac_vector_t& input) {
		reduced_dirac_vector_t d_output, d_input;
		input.copy_to(d_input);

		this->multiply(d_output, d_input);

		output = d_output;
	}

	//Dummy function to be overridden in the most relevant time-consuming cases
	virtual void multiplyAdd(reduced_soa_dirac_vector_t& output, const reduced_soa_dirac_vector_t& vector1, const reduced_soa_dirac_vector_t & vector2, const complex& alpha) {
		if (vector1 == vector2) {
			reduced_dirac_vector_t d_output, d_vector;
			vector1.copy_to(d_vector);
		

			this->multiplyAdd(d_output, d_vector, d_vector, alpha);

			output = d_output;
		} 
		else {
			reduced_dirac_vector_t d_output, d_vector1, d_vector2;
			vector1.copy_to(d_vector1);
			vector2.copy_to(d_vector2);
		

			this->multiplyAdd(d_output, d_vector1, d_vector2, alpha);

			output = d_output;
		}
	}
#endif
	
  	/**
	 * This routine multiplies the Dirac operator to input and stores the result in output
	 * @param output
	 * @param input
	 */
   	 virtual void multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) = 0;

    	/**
	 * This routine multiplies the Dirac operator to vector1 and stores the result in output adding alpha*vector2
	 * @param output
	 * @param vector1
	 * @param vector2
	 * @param alpha
	 */
    virtual void multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t & vector2, const complex& alpha) =0;

    virtual FermionForce* getForce() const = 0;

	/**
	 * This function sets the kappa parameter of the operator
	 * @param _kappa
	 */
	virtual void setKappa(real_t _kappa);

	/**
	 * This function returns back the kappa parameter of the operator
	 * @return kappa
	 */
	real_t getKappa() const;

	/**
	 * This function returns back the pointer to the lattice
	 * @return the pointer to the lattice
	 */
	const reduced_fermion_lattice_t *getLattice() const;

	/**
	 * This function set the lattice
	 * @param lattice
	 */
	virtual void setLattice(const extended_fermion_lattice_t& _lattice);

	void setGamma5(bool _gamma5);
	bool getGamma5() const;
	
	std::string getName() const;

protected:
	reduced_fermion_lattice_t lattice;

	real_t kappa;

	bool gamma5;

	std::string name;
};

} /* namespace Update */
#endif /* DIRACOPERATOR_H_ */
