#ifndef BLOCKIMPROVEDDIRACWILSONOPERATOR_H_
#define BLOCKIMPROVEDDIRACWILSONOPERATOR_H_
#include "BlockDiracOperator.h"
#include "utils/RandomSeed.h"

namespace Update {

class BlockImprovedDiracWilsonOperator : public BlockDiracOperator {
public:
	BlockImprovedDiracWilsonOperator(Color _color = Black);
	BlockImprovedDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa = 0., Color _color = Black);
	~BlockImprovedDiracWilsonOperator();

	virtual void multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input);

	virtual void multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const std::complex<real_t>& alpha);

	virtual FermionForce* getForce() const;

	virtual void setLattice(const extended_fermion_lattice_t& _lattice);

	real_t getCSW() const;
	void setCSW(real_t _csw);
private:
	Color color;

	//The clover term
	real_t csw;
	//The field strength
	reduced_field_strength_lattice_t F;

	void updateFieldStrength(const extended_fermion_lattice_t& _lattice);
};

} /* namespace Update */
#endif /* RANDOMDIRACWILSONOPERATOR_H_ */
