#ifdef LAPLACIANCONJUGATEGRADIENT_H
#define LAPLACIANCONJUGATEGRADIENT_H

namespace Update {

class LaplacianConjugateGradient {
public:
	LaplacianConjugateGradient();
	~LaplacianConjugateGradient();	

	bool solve(const reduced_fermion_lattice_t& lattice, const reduced_color_vector_t& source, reduced_color_vector_t& solution, const block& s, int j_decay, const real_t& shift);

	void setPrecision(real_t _epsilon);
	real_t getPrecision() const;

	real_t getLastError() const;
	unsigned int getLastSteps() const;

	void setMaximumSteps(unsigned int _maxSteps);
	unsigned int getMaximumSteps() const;

private:
	reduced_color_vector_t p;
	reduced_color_vector_t r;
	reduced_color_vector_t tmp;

	real_t epsilon;
	real_t lastError;
	unsigned int lastSteps;
	unsigned int maxSteps;
};


}

#endif

