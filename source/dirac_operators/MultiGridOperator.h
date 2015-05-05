#ifndef MULTIGRIDOPERATOR_H
#define MULTIGRIDOPERATOR_H
#include <vector>
#include "DiracOperator.h"

namespace Update {

class MultiGridVectorLayout {
public:
	static void initialize();
	static void setBasisDimension(int _basisDimension);
	
	static int xBlockSize;
	static int yBlockSize;
	static int zBlockSize;
	static int tBlockSize;

	static int totalNumberOfBlocks;
	static int basisDimension;
	static int size;

	static int index(int site) {
		return (*blockIndex)[site];
	}
private:
	
	static reduced_index_lattice_t* blockIndex;
};

template<typename T,typename TLayout> class MultiGridVector {
public:
	MultiGridVector() : data(new T[TLayout::size]) { }
	MultiGridVector(const MultiGridVector& snd) : data(new T[TLayout::size]) {
		for (int i = 0; i < Layout::size; ++i) data[i] = snd.data[i];
	}
	~MultiGridVector() {
		delete[] data;
	}

	MultiGridVector& operator=(const MultiGridVector& snd) {
		data = new T[TLayout::size];
		for (int i = 0; i < Layout::size; ++i) data[i] = snd.data[i];
		return *this;
	}

	T& operator[](int index) {
		return data[index];
	}

	const T& operator[](int index) const {
		return data[index];
	}

	T& operator()(int vector, int site) {
		return data[Layout::totalNumberOfBlocks*vector + Layout::index(site)];
	}

	const T& operator()(int vector, int site) const {
		return data[Layout::totalNumberOfBlocks*vector + Layout::index(site)];
	}
	
	typedef TLayout Layout;

	vector_t asVector() {
		vector_t result(Layout::size);
		for (int i = 0; i < Layout::size; ++i) {
			result[i] = data[i];
		}
		return result;
	}

private:
	T* data;
};

typedef MultiGridVector< std::complex<real_t> , MultiGridVectorLayout > multigrid_vector_t;

class MultiGridOperator {
public:
	MultiGridOperator();

	void multiply(multigrid_vector_t& output, const multigrid_vector_t& input);

	//The add operation is applied only to the inner Dirac operator
	void multiplyAdd(multigrid_vector_t& output, const multigrid_vector_t& input, const complex& alpha);

	//To set the Dirac operator
	void setDiracOperator(DiracOperator* _dirac) {
		dirac = _dirac;
	}

	//Add a vector to the vector space used for block deflation
	void addVector(const reduced_dirac_vector_t& vector) {
		vectorspace.push_back(&vector);
	}

	//Remove all the vectors in the vector space
	void clearVectorSpace() {
		vectorspace.clear();
	}

	matrix_t asMatrix();

protected:
	std::vector<reduced_dirac_vector_t const*> vectorspace;
	DiracOperator* dirac;

	reduced_dirac_vector_t tmp_input;
	reduced_dirac_vector_t tmp_output;
};

class MultiGridProjector {
public:
	MultiGridProjector();

	void apply(multigrid_vector_t& output, const reduced_dirac_vector_t& input);

	void apply(reduced_dirac_vector_t& output, const multigrid_vector_t& input);

	//Add a vector to the vector space used for block deflation
	void addVector(const reduced_dirac_vector_t& vector) {
		vectorspace.push_back(&vector);
	}

	//Remove all the vectors in the vector space
	void clearVectorSpace() {
		vectorspace.clear();
	}

protected:
	std::vector<reduced_dirac_vector_t const*> vectorspace;

	reduced_dirac_vector_t tmp_input;
	reduced_dirac_vector_t tmp_output;
};

class BlockBasis {
public:
	BlockBasis(unsigned int _basisDimension);
	BlockBasis(const BlockBasis& toCopy);
	~BlockBasis();

	void orthogonalize();

	void orthogonalize(unsigned int index);

	reduced_dirac_vector_t& operator[](int index) {
		return vectorspace[index];
	}

	const reduced_dirac_vector_t& operator[](int index) const {
		return vectorspace[index];
	}

	int size() const {
		return basisDimension;
	}
	
private:
	unsigned int basisDimension;
	reduced_dirac_vector_t* vectorspace;
};

class MultiGridConjugateGradientSolver {
public:
	MultiGridConjugateGradientSolver() : precision(0.000001), maxSteps(300) { }

	bool solve(MultiGridOperator* multiGridOperator, const multigrid_vector_t& source_hat, multigrid_vector_t& solution_hat) {
		solution_hat = source_hat;
		multiGridOperator->multiply(tmp_hat,solution_hat);

#pragma omp parallel for
		for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
			r_hat[m] = source_hat[m] - tmp_hat[m];
			p_hat[m] = r_hat[m];
		}

		long_real_t norm = 0.;
#pragma omp parallel for reduction(+:norm)
		for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
			norm += real(conj(r_hat[m])*r_hat[m]);
		}
		reduceAllSum(norm);
		long_real_t norm_next = norm;

		for (unsigned int innerStep = 0; innerStep < maxSteps; ++innerStep) {
			multiGridOperator->multiply(tmp_hat,p_hat);
			norm = norm_next;
			long_real_t gammaRe = 0;
			long_real_t gammaIm = 0;
#pragma omp parallel for reduction(+:gammaRe,gammaIm)
			for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
				std::complex<real_t> result = conj(p_hat[m])*tmp_hat[m];
				gammaRe += real(result);
				gammaIm += imag(result);
			}
			reduceAllSum(gammaRe);
			reduceAllSum(gammaIm);
			std::complex<real_t> alpha = static_cast<real_t>(norm)/std::complex<real_t>(gammaRe,gammaIm);


#pragma omp parallel for
			for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
				solution_hat[m] = solution_hat[m] + alpha * p_hat[m];
				r_hat[m] = r_hat[m] - alpha * tmp_hat[m];
			}

			norm_next = 0.;
#pragma omp parallel for reduction(+:norm_next)
			for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
				norm_next += real(conj(r_hat[m])*r_hat[m]);
			}
			//Check for convergence
			if (norm_next < precision && innerStep > 5) {
				std::cout << "Inner convergence in " << innerStep << " steps."<< std::endl;
				break;
			} else if (innerStep == maxSteps - 1) {
				std::cout << "Failure in finding convergence in inner solver, last error " << norm_next << std::endl;
			}
			//std::cout << "Inner norm at step " << innerStep << ": " << norm_next << std::endl; 

			real_t beta = static_cast<real_t>(norm_next/norm);

#pragma omp parallel for
			for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
				p_hat[m] = r_hat[m] + beta * p_hat[m];
			}
		}

	}

	void setPrecision(const real_t& _precision) {
		precision = _precision;
	}

	void setMaximumSteps(int _maxSteps) {
		maxSteps = _maxSteps;
	}
private:
	real_t precision;
	int maxSteps;

	multigrid_vector_t r_hat;
	multigrid_vector_t p_hat;
	multigrid_vector_t tmp_hat;
};

class MultiGridBiConjugateGradientSolver {
public:
	MultiGridBiConjugateGradientSolver() : precision(0.000001), maxSteps(300) { }

	bool solve(MultiGridOperator* multiGridOperator, const multigrid_vector_t& source, multigrid_vector_t& solution) {
		solution = source;
		//Use p as temporary vector
		multiGridOperator->multiply(p,solution);
		
#pragma omp parallel for
		for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
			//Set the initial residual to source-A.solution and residual_hat accordingly
			residual[m] = source[m] - p[m];
			residual_hat[m] = source[m] + p[m];
			//Set then nu and p to zero
			p[m] = 0;
			nu[m] = 0;
		}

		//Set the initial parameter of the program
		std::complex<real_t> alpha = 1., omega = 1.;
		std::complex<long_real_t> rho = 1.;
		unsigned int step = 0;
		long_real_t norm_r = 0.;

		while (step < maxSteps) {
			//rho[k] = rhat.r[k-1]
			long_real_t rho_next_re = 0.;
			long_real_t rho_next_im = 0.;
#pragma omp parallel for reduction(+:rho_next_re, rho_next_im)
			for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
				complex partial = conj(residual_hat[m])*residual[m];
				rho_next_re += real(partial);
				rho_next_im += imag(partial);
			}
			reduceAllSum(rho_next_re);
			reduceAllSum(rho_next_im);

			std::complex<long_real_t> rho_next(rho_next_re,rho_next_im);

			if (norm(rho_next) == 0.) {
				if (isOutputProcess()) std::cout << "BiConjugateGradient::Fatal error in norm " << rho_next << " at step " << step << std::endl;
				return false;//TODO
			}

			std::complex<real_t> beta = static_cast< std::complex<real_t> >((rho_next/rho))*(alpha/omega);
			//p = r[[k - 1]] + beta*(p[[k - 1]] - omega[[k - 1]]*nu[[k - 1]])
#pragma omp parallel for
			for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
				p[m] = residual[m] + beta*(p[m] - omega*nu[m]);
			}

			//nu = A.p[[k]]
			multiGridOperator->multiply(nu,p);

			//alpha = rho[[k]]/(rhat[[1]].nu[[k]]);
			long_real_t alphatmp_re = 0.;
			long_real_t alphatmp_im = 0.;
#pragma omp parallel for reduction(+:alphatmp_re, alphatmp_im)
			for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
				complex partial = conj(residual_hat[m])*nu[m];
				alphatmp_re += real(partial);
				alphatmp_im += imag(partial);
			}
			reduceAllSum(alphatmp_re);
			reduceAllSum(alphatmp_im);

			std::complex<long_real_t> alphatmp(alphatmp_re,alphatmp_im);
			alpha = static_cast< std::complex<real_t> >(rho_next/alphatmp);

			//s = r[[k - 1]] - alpha*nu[[k]]
#pragma omp parallel for
			for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
				s[m] = residual[m] - alpha *(nu[m]);
			}
			//s.updateHalo();

			//t = A.s;
			multiGridOperator->multiply(t,s);

			//omega = (t.s)/(t.t)
			long_real_t tmp1_re = 0., tmp1_im = 0., tmp2_re = 0., tmp2_im = 0.;
#pragma omp parallel for reduction(+:tmp1_re, tmp1_im, tmp2_re, tmp2_im)
			for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
				complex partial1 = conj(t[m])*s[m];
				tmp1_re += real(partial1);
				tmp1_im += imag(partial1);
				complex partial2 = conj(t[m])*t[m];
				tmp2_re += real(partial2);
				tmp2_im += imag(partial2);
			}
			reduceAllSum(tmp1_re);
			reduceAllSum(tmp1_im);
			reduceAllSum(tmp2_re);
			reduceAllSum(tmp2_im);

			std::complex<long_real_t> tmp1(tmp1_re, tmp1_im), tmp2(tmp2_re, tmp2_im);
			omega = static_cast< std::complex<real_t> >(tmp1/tmp2);

			if (real(tmp2) == 0) {
#pragma omp parallel for
				for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
					solution[m] = source[m];
				}
				//solution.updateHalo();
				return true;//TODO, identity only?
			}

			//solution[[k]] = solution[[k - 1]] + alpha*p[[k]] + omega[[k]]*s
#pragma omp parallel for
			for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
				solution[m] += alpha*(p[m]) + omega*(s[m]);
			}
			//solution.updateHalo();

			//residual[[k]] = s - omega[[k]]*t
			//norm = residual[[k]].residual[[k]]
			norm_r = 0.;
#pragma omp parallel for reduction(+:norm_r)
			for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
				residual[m] = s[m] - omega*(t[m]);
				norm_r += real(conj(residual[m])*residual[m]);
			}
			reduceAllSum(norm_r);


			if (norm_r < precision && step > 4) {
				//lastSteps = step;
//#ifdef BICGLOG
				if (isOutputProcess()) std::cout << "MultiGridBiConjugateGradientSolver::Convergence in " << step << " - final error norm: " << norm_r << std::endl;
//#endif
				return true;
			}
			//#ifdef BICGLOG
			//		else if (isOutputProcess()) std::cout << "Error at step " << step << ": " << real(norm) << std::endl;
			//#endif

			rho = rho_next;

			++step;
		}

		if (isOutputProcess()) std::cout << "MultiGridBiConjugateGradientSolver::Failure in finding convergence in " << maxSteps << " - final error norm: " << norm_r << std::endl;

		return false;
	}

	void setPrecision(const real_t& _precision) {
		precision = _precision;
	}

	void setMaximumSteps(int _maxSteps) {
		maxSteps = _maxSteps;
	}
private:
	real_t precision;
	int maxSteps;

	multigrid_vector_t residual;
	multigrid_vector_t residual_hat;
	multigrid_vector_t p;
	multigrid_vector_t nu;
	multigrid_vector_t s;
	multigrid_vector_t t;
};


class MGLeftProjector {
public:
	MGLeftProjector() : mgSolver(new MultiGridBiConjugateGradientSolver()) {
		mgSolver->setMaximumSteps(50);
		mgSolver->setPrecision(1.e-11);
	}
	~MGLeftProjector() {
		delete mgSolver;
	}

	void apply(DiracOperator* dirac, MultiGridOperator* multiGridOperator, MultiGridProjector* multiGridProjector, reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
		multigrid_vector_t solution_hat, source_hat;
		multiGridProjector->apply(source_hat,input);

		mgSolver->solve(multiGridOperator, source_hat, solution_hat);

		multiGridProjector->apply(output,solution_hat);

		dirac->multiply(tmp,output);

#pragma omp parallel for
		for (int site = 0; site < output.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				output[site][mu] = input[site][mu] - tmp[site][mu];
			}
		}
	}
private:
	MGLeftProjector(const MGLeftProjector&);
	
	MultiGridBiConjugateGradientSolver* mgSolver;
	reduced_dirac_vector_t tmp;
};


class MGRightProjector {
public:
	MGRightProjector() : mgSolver(new MultiGridBiConjugateGradientSolver()) {
		mgSolver->setMaximumSteps(200);
		mgSolver->setPrecision(1.e-11);
	}
	~MGRightProjector() {
		delete mgSolver;
	}

	void apply(DiracOperator* dirac, MultiGridOperator* multiGridOperator, MultiGridProjector* multiGridProjector, reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
		dirac->multiply(output,input);
		multigrid_vector_t solution_hat, source_hat;
		multiGridProjector->apply(source_hat,output);

		mgSolver->solve(multiGridOperator, source_hat, solution_hat);

		multiGridProjector->apply(output,solution_hat);

#pragma omp parallel for
		for (int site = 0; site < output.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				output[site][mu] = input[site][mu] - output[site][mu];
			}
		}
	}
private:
	MGRightProjector(const MGLeftProjector&);
	
	MultiGridBiConjugateGradientSolver* mgSolver;
};

class MGDeflatedDirac : public DiracOperator {
public:
	MGDeflatedDirac(DiracOperator* _dirac, MGLeftProjector* _leftProjector, MultiGridOperator* _mgOperator, MultiGridProjector* _multiGridProjector) : DiracOperator(), dirac(_dirac), leftProjector(_leftProjector), mgOperator(_mgOperator), multiGridProjector(_multiGridProjector) { }

	void multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
		dirac->multiply(tmp,input);
		leftProjector->apply(dirac,mgOperator,multiGridProjector,output,tmp);
	}

	void multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const std::complex<real_t>& alpha) {
	}

	FermionForce* getForce() const {
		return NULL;
	}
private:
	DiracOperator* dirac;
	MGLeftProjector* leftProjector;
	MultiGridOperator* mgOperator;
	MultiGridProjector* multiGridProjector;

	reduced_dirac_vector_t tmp;
	
};

}


#endif
