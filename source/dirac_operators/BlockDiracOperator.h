#ifndef BLOCKDIRACOPERATOR_H_
#define BLOCKDIRACOPERATOR_H_
#include "DiracOperator.h"
#include "utils/RandomSeed.h"

namespace Update {

enum Color {Black = 0, Red};

#ifdef ENABLE_MPI
typedef Lattice::Lattice<short int, Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_index_lattice_t;
#endif
#ifndef ENABLE_MPI
typedef Lattice::Lattice<short int, Lattice::LocalLayout > reduced_index_lattice_t;
#endif

class BlockDiracOperator : public DiracOperator {
public:
	BlockDiracOperator(Color _color = Black, real_t _twist = 0.);
	BlockDiracOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa, Color _color = Black, real_t _twist = 0.);
	~BlockDiracOperator();

	static BlockDiracOperator* getInstance(const std::string& name, unsigned int power, const StorageParameters& parameters, Color _color = Black);

	virtual void setBlockSize(const std::vector<unsigned int>& _blockSize);
	std::vector<unsigned int> getBlockSize() const;

	void project(reduced_dirac_vector_t& output);

	void setTwist(real_t _twist);
protected:
	int xBlockSize;
	int yBlockSize;
	int zBlockSize;
	int tBlockSize;
	Color color;
	real_t twist;
	
	reduced_index_lattice_t index_lattice;
	void initialize_index_lattice();
};

} /* namespace Update */
#endif /* RANDOMDIRACWILSONOPERATOR_H_ */
