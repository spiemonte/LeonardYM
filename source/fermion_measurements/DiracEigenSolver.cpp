/*
 * DiracEigenSolver.cpp
 *
 *  Created on: Jun 26, 2012
 *      Author: spiem_01
 */

#include "DiracEigenSolver.h"
#include "algebra_utils/AlgebraUtils.h"
#include "inverters/BiConjugateGradient.h"
#ifdef EIGEN
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#endif
#include <algorithm>
#include "utils/ToString.h"
#include "dirac_functions/ChebyshevRecursion.h"

namespace Update {

const std::complex<real_t> EigenvaluesMap[4] = {1., -1., std::complex<real_t>(0.,1.), std::complex<real_t>(0.,-1.)};

bool maxcomparison(const std::complex<real_t>& i, const std::complex<real_t>& j) { return (abs(i)>abs(j)); }
bool mincomparison(const std::complex<real_t>& i, const std::complex<real_t>& j) { return (abs(i)<abs(j)); }

inline void rotateVector(reduced_dirac_vector_t& vector, EigevaluesMode mode) {
	typedef reduced_dirac_vector_t::Layout Layout;
	if (mode != LargestReal) {
#pragma omp parallel for
		for (int site = 0; site < vector.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) vector[site][mu] = EigenvaluesMap[mode]*vector[site][mu];
		}
	}
}

DiracEigenSolver::DiracEigenSolver() : epsilon(0.00000001), inverterPrecision(0.0000000001), inverterMaximumSteps(10000), extra_steps(250), useChebyshev(false), maximalNumberOfRestarts(50), biConjugateGradient(0), chebyshevRecursion(new ChebyshevRecursion(0.2, 7., 15)) { }

DiracEigenSolver::DiracEigenSolver(const DiracEigenSolver& copy) : epsilon(copy.epsilon), inverterPrecision(copy.inverterPrecision), inverterMaximumSteps(copy.inverterMaximumSteps), extra_steps(copy.extra_steps), useChebyshev(copy.useChebyshev), maximalNumberOfRestarts(copy.maximalNumberOfRestarts), biConjugateGradient(0), chebyshevRecursion(new ChebyshevRecursion(*copy.chebyshevRecursion)) { }

DiracEigenSolver::~DiracEigenSolver() {
	if (biConjugateGradient) delete biConjugateGradient;
	if (chebyshevRecursion) delete chebyshevRecursion;
}

ChebyshevRecursion* DiracEigenSolver::getChebyshevRecursion() const {
	return chebyshevRecursion;
}

void DiracEigenSolver::setTolerance(const real_t& precision) {
	epsilon = precision;
}

real_t DiracEigenSolver::getTolerance() const {
	return epsilon;
}

void DiracEigenSolver::setInverterPrecision(const real_t& precision) {
	inverterPrecision = precision;
}

real_t DiracEigenSolver::getInverterPrecision() const {
	return inverterPrecision;
}

void DiracEigenSolver::setUseChebyshev(bool _useChebyshev) {
	useChebyshev = _useChebyshev;
}

bool DiracEigenSolver::getUseChebyshev() const {
	return useChebyshev;
}

void DiracEigenSolver::setMaximalNumberOfRestarts(unsigned int _maximalNumberOfRestarts) {
	maximalNumberOfRestarts = _maximalNumberOfRestarts;
}

unsigned int DiracEigenSolver::getMaximalNumberOfRestarts() const {
	return maximalNumberOfRestarts;
}

void DiracEigenSolver::extendArnoldi(DiracOperator* diracOperator, std::vector<reduced_dirac_vector_t>& V, reduced_dirac_vector_t& f,  matrix_t& H, unsigned int m, unsigned int k, EigevaluesMode mode) {
	reduced_dirac_vector_t w, tmp;

	for (unsigned int j = m; j < k-1; ++j) {
		//beta = norm(f)
		real_t beta = sqrt(AlgebraUtils::squaredNorm(f));
		//V[j] = f/beta
#pragma omp parallel for
		for (int site = 0; site < f.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				V[j+1][site][mu] = f[site][mu]/beta;
			}
		}
		V[j+1].updateHalo();
		//H(j+1,j) = beta
		H(j+1,j) = beta;
		//w = D.V[j+1]
		if (useChebyshev) chebyshevRecursion->evaluate(diracOperator, w, V[j+1]);
		else {
			tmp = V[j+1];
			rotateVector(tmp, mode);
			diracOperator->multiplyAdd(w,tmp,V[j+1],+7.);
		}
		
		//Gram schimdt
#pragma omp parallel for
		for (int site = 0; site < w.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				f[site][mu] = w[site][mu];
			}
		}
		f.updateHalo();
		for (unsigned int i = 0; i <= j+1; ++i) {
			std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(V[i],w));
			H(i,j+1) = proj;
#pragma omp parallel for
			for (int site = 0; site < f.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					f[site][mu] -= proj*V[i][site][mu];
				}
			}
			f.updateHalo();//TODO maybe not needed
		}
		//More stable Gram schimdt
		for (unsigned int i = 0; i <= j+1; ++i) {
			std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(V[i],f));
			H(i,j+1) += proj;
#pragma omp parallel for
			for (int site = 0; site < f.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					f[site][mu] -= proj*V[i][site][mu];
				}
			}
			f.updateHalo();//TODO maybe not needed
		}
	}
}

void DiracEigenSolver::startArnoldi(DiracOperator* diracOperator, std::vector<reduced_dirac_vector_t>& V, reduced_dirac_vector_t& f,  matrix_t& H, EigevaluesMode mode) {
	reduced_dirac_vector_t tmp, w;	

	AlgebraUtils::generateRandomVector(tmp);
	if (useChebyshev) chebyshevRecursion->evaluate(diracOperator,V[0],tmp);
	else V[0] = tmp;
	AlgebraUtils::normalize(V[0]);
	
	//w = D.V[0]	
	if (useChebyshev) chebyshevRecursion->evaluate(diracOperator, w, V[0]);
	else {
		tmp = V[0];
		rotateVector(tmp, mode);
		diracOperator->multiplyAdd(w,tmp,V[0],+7.);
	}
	
	std::complex<real_t> alpha = static_cast< std::complex<real_t> >(AlgebraUtils::dot(V[0], w));

	set_to_zero(H);
	//f = w - alpha*V[0]
#pragma omp parallel for
	for (int site = 0; site < w.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			f[site][mu] = w[site][mu] - alpha*V[0][site][mu];
		}
	}
	f.updateHalo();//TODO maybe not needed
	H(0,0) = alpha;
}

long_real_t DiracEigenSolver::finishArnoldi(DiracOperator* diracOperator, const std::vector<reduced_dirac_vector_t>& V, const matrix_t& H, std::vector< std::complex<real_t> >& eigenvalues, std::vector<reduced_dirac_vector_t>& eigenvectors, unsigned int steps, EigevaluesMode mode) {
	typedef reduced_dirac_vector_t::Layout Layout;

	eigenvalues.resize(steps);

	matrix_t eigvec(steps,steps);
#ifdef EIGEN
	Eigen::ComplexEigenSolver<matrix_t> ces(H, true);
	for (unsigned int i = 0; i < steps; ++i) {
		eigenvalues[i] = EigenvaluesMap[mode]*(ces.eigenvalues()[i] - static_cast<real_t>(7.));
	}
	eigvec = ces.eigenvectors();
#endif
#ifdef ARMADILLO
	vector_t eigval(steps);
	arma::eig_gen(eigval, eigvec, H);
	for (unsigned int i = 0; i < steps; ++i) {
		eigenvalues[i] = EigenvaluesMap[mode]*(eigval[i] - static_cast<real_t>(7.));
	}
#endif

	eigenvectors.resize(steps);
	for (unsigned int i = 0; i < steps; ++i) {
		AlgebraUtils::setToZero(eigenvectors[i]);
		for (unsigned int j = 0; j < steps; ++j) {
#pragma omp parallel for
			for (int site = 0; site < Layout::localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					eigenvectors[i][site][mu] += eigvec.at(j,i)*V[j][site][mu];
				}
			}
		}
		eigenvectors[i].updateHalo();
	
	}

	reduced_dirac_vector_t tmp, tmpe;

	long_real_t maximal_difference = 0.;

	for (unsigned int i = 0; i < steps; ++i) {
		diracOperator->multiply(tmp,eigenvectors[i]);
		eigenvalues[i] = AlgebraUtils::dot(eigenvectors[i],tmp)/AlgebraUtils::squaredNorm(eigenvectors[i]);

		if (i >=  extra_steps) {
			//Now we check the convergence
			long_real_t diffnorm = 0.;

#pragma omp parallel for reduction(+:diffnorm)
			for (int site = 0; site < Layout::localsize; ++site) {
				real_t sum = 0.;
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (unsigned int c = 0; c < diracVectorLength; ++c) {
						sum += std::fabs(tmp[site][mu][c] - eigenvalues[i]*eigenvectors[i][site][mu][c]);
					}
				}
				diffnorm += sum;
			}
			reduceAllSum(diffnorm);

			if (isOutputProcess()) {
				std::cout << "DiracEigenSolver::Convergence precision for eigenvalue " << i - extra_steps << " lambda=" << eigenvalues[i] << ": " << diffnorm;
				if (diffnorm > epsilon) std::cout << " > " << epsilon << std::endl;
				else std::cout << std::endl;
			}
		
			if (diffnorm > maximal_difference) maximal_difference = diffnorm;
		}
	}

	std::reverse(eigenvectors.begin(),eigenvectors.end());
	std::reverse(eigenvalues.begin(),eigenvalues.end());

	eigenvalues.erase(eigenvalues.end() - extra_steps, eigenvalues.end());
	eigenvectors.erase(eigenvectors.end() - extra_steps, eigenvectors.end());

	return maximal_difference;
	
/*
	//Now we check the convergence
	diracOperator->multiply(tmp,eigenvectors.back());
#pragma omp parallel for
	for (int site = 0; site < Layout::localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			tmpe[site][mu] = eigenvalues.back()*eigenvectors.back()[site][mu];
		}
	}

	long_real_t diffnorm = AlgebraUtils::differenceNorm(tmp,tmpe);

	if (isOutputProcess()) std::cout << "DiracEigenSolver::Convergence precision for the largest eigenvalue: " << abs(diffnorm) << std::endl;

	//Now we check the convergence
	diracOperator->multiply(tmp,eigenvectors[eigenvectors.size()-extra_steps]);
#pragma omp parallel for
	for (int site = 0; site < Layout::localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			tmpe[site][mu] = eigenvalues[eigenvalues.size()-extra_steps]*eigenvectors[eigenvectors.size()-extra_steps][site][mu];
		}
	}

	diffnorm = AlgebraUtils::differenceNorm(tmp,tmpe);

	if (isOutputProcess()) std::cout << "DiracEigenSolver::Convergence precision for the smallest eigenvalue: " << abs(diffnorm) << ", target: " << epsilon << std::endl;

	std::reverse(eigenvectors.begin(),eigenvectors.end());
	std::reverse(eigenvalues.begin(),eigenvalues.end());

	eigenvalues.erase(eigenvalues.end() - extra_steps, eigenvalues.end());
	eigenvectors.erase(eigenvectors.end() - extra_steps, eigenvectors.end());

	return diffnorm;*/
}

bool compare(const std::complex<real_t>& x, const std::complex<real_t>& y) { return (real(x)>real(y)); }

void checkArnoldi(DiracOperator* diracOperator, const std::vector<reduced_dirac_vector_t>& V, const reduced_dirac_vector_t& f,  const matrix_t& H, unsigned int index) {
	std::vector<reduced_dirac_vector_t> VD(V.size()), VH(V.size());
	for (unsigned int i = 0; i < V.size(); ++i) {
		reduced_dirac_vector_t tmp = V[i];
		rotateVector(tmp, SmallestReal);
		//diracOperator->multiplyAdd(VD[i],tmp,V[i],+7.);
		diracOperator->multiply(VD[i],tmp);
		//chebyshevRecursion->evaluate(diracOperator,VD[i],VD[i]);
	}

	for (unsigned int i = 0; i < V.size(); ++i) {
		AlgebraUtils::setToZero(VH[i]);// U(k,0)*V[0];
		for (unsigned int j = 0; j < V.size(); ++j) {
#pragma omp parallel for
			for (int site = 0; site < VH[j].completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					VH[i][site][mu] += H(j,i)*V[j][site][mu];
				}
			}
		}
	}
	
	for (unsigned int i = 0; i < index - 1; ++i) {
		long_real_t diffnorm = AlgebraUtils::differenceNorm(VD[i],VH[i]);
		if (isOutputProcess()) std::cout << "Check on vector " << i << " " << diffnorm << std::endl;
	}

#pragma omp parallel for
	for (int site = 0; site < VH[index-1].completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			VH[index-1][site][mu] += f[site][mu];
		}
	}

	long_real_t diffnorm = AlgebraUtils::differenceNorm(VD[index-1],VH[index-1]);
	if (isOutputProcess()) std::cout << "Check on vector " << index-1 << " " << diffnorm << std::endl;
}

void DiracEigenSolver::forceHermitianPairing(std::vector<reduced_dirac_vector_t>& V, std::vector< std::complex<real_t> >& eigenvalues) {

	std::cout << "Ma scusa: " << std::endl << V[0][0][0] << " - " << std::endl << V[0][0][1] << " - " << std::endl << V[0][0][2] << " - " << std::endl << V[0][0][3] << std::endl; 
	std::cout << "Ma scusa: " << std::endl << V[1][0][0] << " - " << std::endl << V[1][0][1] << " - " << std::endl << V[1][0][2] << " - " << std::endl << V[1][0][3] << std::endl;
	for (unsigned int i = 0; i < V.size(); i += 2) {
#pragma omp parallel for
		for (int site = 0; site < V[i].completesize; ++site) {
			for (unsigned i = 0; i < diracVectorLength; ++i) {
				V[i+1][site][0][i] = conj(-V[i][site][1][i]);
				V[i+1][site][1][i] = conj(V[i][site][0][i]);
				V[i+1][site][2][i] = conj(-V[i][site][3][i]);
				V[i+1][site][3][i] = conj(V[i][site][2][i]);
			}
		}

		eigenvalues[i+1] = eigenvalues[i];
	}
}

void DiracEigenSolver::restartArnoldi(DiracOperator* diracOperator, std::vector<reduced_dirac_vector_t>& V, reduced_dirac_vector_t& f,  matrix_t& H, unsigned int extra_steps) {
	unsigned int steps = H.rows();

	//std::cout << "H: " << toString< std::complex<real_t> >(H) << std::endl;

	//checkArnoldi(diracOperator, V, f,  H, V.size());

	std::vector< std::complex<real_t> > eigenvalues(steps);

	matrix_t eigvec(steps,steps);
#ifdef EIGEN
	Eigen::ComplexEigenSolver<matrix_t> ces(H, true);
	for (unsigned int i = 0; i < steps; ++i) {
		eigenvalues[i] = ces.eigenvalues()[i];
	}
	eigvec = ces.eigenvectors();
#endif
#ifdef ARMADILLO
	vector_t eigval(steps);
	arma::eig_gen(eigval, eigvec, H);
	for (unsigned int i = 0; i < steps; ++i) {
		eigenvalues[i] = eigval[i];
	}
#endif

	//Sort by largest real
	std::sort(eigenvalues.begin(), eigenvalues.end(), compare);
	matrix_t identity(steps,steps);
	for (unsigned int i = 0; i < steps; ++i) {
		for (unsigned int j = 0; j < steps; ++j) {
			identity(i,j) = 0.;
		}
		identity(i,i) = 1.;
	}

	std::vector<reduced_dirac_vector_t> Vp(V.size());

	
	for (unsigned int j = steps - 1; j >= steps - extra_steps; --j) {
		matrix_t Q, R;
#ifdef EIGEN
		matrix_t Hs = H - eigenvalues[j]*identity;
		Eigen::HouseholderQR<matrix_t> qr_decomposition(Hs);
		R = qr_decomposition.matrixQR().triangularView<Eigen::Upper>();
		Q = qr_decomposition.householderQ();
#endif

		H = R*Q + eigenvalues[j]*identity;

		for (unsigned int k = 0; k < V.size(); ++k) {
			AlgebraUtils::setToZero(Vp[k]);// U(k,0)*V[0];
			for (unsigned int l = 0; l < V.size(); ++l) {
				if (fabs(Q(l,k)) > 10e-17) {
#pragma omp parallel for
					for (int site = 0; site < Vp[k].completesize; ++site) {
						for (unsigned int mu = 0; mu < 4; ++mu) {
							Vp[k][site][mu] += Q(l,k)*V[l][site][mu];
						}
					}
				}
			}
		}

		V = Vp;
		
#pragma omp parallel for
		for (int site = 0; site <f.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				f[site][mu] = V[j][site][mu] * H(j, j - 1) + conj(Q(j, j - 1))*f[site][mu];
			}
		}

		for (unsigned int i = 0; i < steps; ++i) {
			H(i,j) = 0;
			H(j,i) = 0;
		}

		//checkArnoldi(diracOperator, V, f, H, j);		
	}
}



void DiracEigenSolver::maximumEigenvalues(DiracOperator* diracOperator, std::vector< std::complex<real_t> >& eigenvalues, std::vector<reduced_dirac_vector_t>& eigenvectors, unsigned int n, EigevaluesMode mode) {
	unsigned int steps = extra_steps + n;
	//The orthonormal vectors generated by the Arnoldi process
	std::vector<reduced_dirac_vector_t> V;
	reduced_dirac_vector_t w, f;
	//Reserve some memory
	V.resize(steps);

	//The upper Hessenberg matrix generated by the Arnoldi process
	matrix_t H(steps,steps);
	
	this->startArnoldi(diracOperator, V, f,  H, mode);
	this->extendArnoldi(diracOperator, V, f, H, 0, steps, mode);

	for (unsigned int m = 0; m < maximalNumberOfRestarts; ++m) {
		this->restartArnoldi(diracOperator, V, f, H, extra_steps);
		this->extendArnoldi(diracOperator, V, f, H, n-1, steps, mode);
		this->restartArnoldi(diracOperator, V, f, H, extra_steps);
		this->extendArnoldi(diracOperator, V, f, H, n-1, steps, mode);
		long_real_t convergence = this->finishArnoldi(diracOperator, V, H, eigenvalues, eigenvectors, steps, mode);
		if (convergence < epsilon) return;
	}

	this->finishArnoldi(diracOperator, V, H, eigenvalues, eigenvectors, steps, mode);
}



/*void DiracEigenSolver::maximumEigenvalues(DiracOperator* diracOperator, std::vector< std::complex<real_t> >& eigenvalues, std::vector<reduced_dirac_vector_t>& eigenvectors, unsigned int n, EigevaluesMode mode, unsigned int restarts) {
	typedef reduced_dirac_vector_t::Layout Layout;

	reduced_dirac_vector_t tmp, tmpe;
	
	unsigned int steps = extra_steps + n;
	//The orthonormal vectors generated by the arnoldi process
	std::vector<reduced_dirac_vector_t> V;
	reduced_dirac_vector_t w, f;
	//Reserve some memory
	V.resize(steps);
	
	AlgebraUtils::generateRandomVector(V[0]);

	for (unsigned int m = 0; m < restarts; ++m) {
		if (m != 0) {
			//Now we check the convergence
			

			//First we project away the eigenvectors that are already converged
			for (unsigned int i = 0; i < converged_eigenvectors.size(); ++i) {
				std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(converged_eigenvectors[i], V[0]));
				
#pragma omp parallel for
				for (int site = 0; site < w.localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						V[0][site][mu] -= proj*converged_eigenvectors[i][site][mu];
					}
				}
			}

			//Then we project away the unwanted eigenvectors
			for (unsigned int i = n; i < n + extra_steps; ++i) {
				std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(eigenvectors[i], V[0]));
				
#pragma omp parallel for
				for (int site = 0; site < w.localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						V[0][site][mu] -= proj*eigenvectors[i][site][mu];
					}
				}
			}

			//Finally we enhance the wanted eigenvectors
			for (unsigned int i = 0; i < n; ++i) {
				std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(eigenvectors[i], V[0]));

				//Now we check the convergence
				diracOperator->multiply(tmp,eigenvectors[i]);
#pragma omp parallel for
				for (int site = 0; site < Layout::localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						tmpe[site][mu] = eigenvalues[i]*eigenvectors[i][site][mu];
					}
				}

				long_real_t diffnorm = AlgebraUtils::differenceNorm(tmp,tmpe);
				if (diffnorm < 0.000000001) {
					if (isOutputProcess()) std::cout << "DiracEigenSolver::Eigenvalue " << i << " converged to " << eigenvalues[i] << " at restart iteration " << m << "." << std::endl;

#pragma omp parallel for
					for (int site = 0; site < w.localsize; ++site) {
						for (unsigned int mu = 0; mu < 4; ++mu) {
							V[0][site][mu] -= proj*eigenvectors[i][site][mu];
						}
					}

					converged_eigenvalues.push_back(eigenvalues[i]);
					converged_eigenvectors.push_back(eigenvectors[i]);

					if (converged_eigenvalues.size() == n) return;
				}
				else {
#pragma omp parallel for
					for (int site = 0; site < w.localsize; ++site) {
						for (unsigned int mu = 0; mu < 4; ++mu) {
							V[0][site][mu] += (proj/pow(eigenvalues[i],2))*eigenvectors[i][site][mu];
						}
					}
				}
			}
		}
	
		AlgebraUtils::normalize(V[0]);
	
		//w = D.V[0]
		reduced_dirac_vector_t tmpm = V[0];
		rotateVector(tmpm, mode);
		diracOperator->multiplyAdd(w,tmpm,V[0],+7.);
	
		std::complex<real_t> alpha = static_cast< std::complex<real_t> >(AlgebraUtils::dot(V[0], w));
		matrix_t H(steps,steps);
		set_to_zero(H);
		//f = w - alpha*V[0]
#pragma omp parallel for
		for (int site = 0; site < w.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				f[site][mu] = w[site][mu] - alpha*V[0][site][mu];
			}
		}
		f.updateHalo();//TODO maybe not needed
		H(0,0) = alpha;
		for (unsigned int j = 0; j < steps - 1; ++j) {
			//beta = norm(f)
			real_t beta = sqrt(AlgebraUtils::squaredNorm(f));
			//V[j] = f/beta
#pragma omp parallel for
			for (int site = 0; site < f.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					V[j+1][site][mu] = f[site][mu]/beta;
				}
			}
			V[j+1].updateHalo();
			//H(j+1,j) = beta
			H(j+1,j) = beta;
			//w = D.V[j+1]
			tmpm = V[j+1];
			rotateVector(tmpm, mode);
			diracOperator->multiplyAdd(w,tmpm,V[j+1],+7.);
			
		
			//Gram schimdt
#pragma omp parallel for
			for (int site = 0; site < w.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					f[site][mu] = w[site][mu];
				}
			}
			f.updateHalo();
			for (unsigned int i = 0; i <= j+1; ++i) {
				std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(V[i],w));
				H(i,j+1) = proj;
#pragma omp parallel for
				for (int site = 0; site < f.localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						f[site][mu] -= proj*V[i][site][mu];
					}
				}
				f.updateHalo();//TODO maybe not needed
			}
			//More stable Gram schimdt
			for (unsigned int i = 0; i <= j+1; ++i) {
				std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(V[i],f));
				H(i,j+1) += proj;
#pragma omp parallel for
				for (int site = 0; site < f.localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						f[site][mu] -= proj*V[i][site][mu];
					}
				}
				f.updateHalo();//TODO maybe not needed
			}
		}

		eigenvalues.resize(steps);

		matrix_t eigvec(steps,steps);
#ifdef EIGEN
		Eigen::ComplexEigenSolver<matrix_t> ces(H, true);
		for (unsigned int i = 0; i < steps; ++i) {
			eigenvalues[i] = EigenvaluesMap[mode]*(ces.eigenvalues()[i] - static_cast<real_t>(7.));
		}
		eigvec = ces.eigenvectors();
#endif
#ifdef ARMADILLO
		vector_t eigval(steps);
		arma::eig_gen(eigval, eigvec, H);
		for (unsigned int i = 0; i < steps; ++i) {
			eigenvalues[i] = EigenvaluesMap[mode]*(eigval[i] - static_cast<real_t>(7.));
		}
#endif

		eigenvectors.resize(steps);
		for (unsigned int i = 0; i < steps; ++i) {
			AlgebraUtils::setToZero(eigenvectors[i]);
			for (unsigned int j = 0; j < steps; ++j) {
#pragma omp parallel for
				for (int site = 0; site < Layout::localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						eigenvectors[i][site][mu] += eigvec.at(j,i)*V[j][site][mu];
					}
				}
			}
			eigenvectors[i].updateHalo();
		
		}
	
		//Now we check the convergence
		diracOperator->multiply(tmp,eigenvectors.back());
#pragma omp parallel for
		for (int site = 0; site < Layout::localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				tmpe[site][mu] = eigenvalues.back()*eigenvectors.back()[site][mu];
			}
		}

		std::complex<long_real_t> diffnorm = AlgebraUtils::differenceNorm(tmp,tmpe);

		if (isOutputProcess()) std::cout << "DiracEigenSolver::Convergence precision at restart " << m << ": " << abs(diffnorm) << std::endl;

		std::reverse(eigenvectors.begin(),eigenvectors.end());
		std::reverse(eigenvalues.begin(),eigenvalues.end());
	}

	for (unsigned int i = 0; i < n - converged_eigenvalues.size(); ++i) {
		converged_eigenvalues.push_back(eigenvalues[i]);
		converged_eigenvectors.push_back(eigenvectors[i]);
	}

	//eigenvalues.erase(eigenvalues.end() - extra_steps, eigenvalues.end());
	//eigenvectors.erase(eigenvectors.end() - extra_steps, eigenvectors.end());
	eigenvalues = converged_eigenvalues;
	eigenvectors = converged_eigenvectors;
}*/

void DiracEigenSolver::minimumEigenvalues(DiracOperator* diracOperator, std::vector< std::complex<real_t> >& eigenvalues, std::vector<reduced_dirac_vector_t>& eigenvectors/*, Polynomial& map*/, unsigned int n/*, int nmode*/) {
	typedef reduced_dirac_vector_t::Layout Layout;

	unsigned int steps = extra_steps + n;
	//The orthonormal vectors generated by the arnoldi process
	std::vector<reduced_dirac_vector_t> V;
	//Reserve some memory
	V.resize(steps);
	AlgebraUtils::generateRandomVector(V[0]);
	AlgebraUtils::normalize(V[0]);
	reduced_dirac_vector_t w, f;
	if (biConjugateGradient == 0) biConjugateGradient = new BiConjugateGradient();
	biConjugateGradient->setPrecision(inverterPrecision);
	biConjugateGradient->setMaximumSteps(inverterMaximumSteps);
	//w = D^(-1).V[0]
	biConjugateGradient->solve(diracOperator,V[0],w);
	std::complex<real_t> alpha = static_cast< std::complex<real_t> >(AlgebraUtils::dot(V[0], w));
	matrix_t H(steps,steps);
	H.zeros();
	//f = w - alpha*V[0]
#pragma omp parallel for
	for (int site = 0; site < w.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			f[site][mu] = w[site][mu] - alpha*V[0][site][mu];
		}
	}
	f.updateHalo();//TODO maybe not needed
	H(0,0) = alpha;
	for (unsigned int j = 0; j < steps - 1; ++j) {
		//beta = norm(f)
		real_t beta = sqrt(AlgebraUtils::squaredNorm(f));
		//V[j] = f/beta
#pragma omp parallel for
		for (int site = 0; site < f.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				V[j+1][site][mu] = f[site][mu]/beta;
			}
		}
		V[j+1].updateHalo();
		//H(j+1,j) = beta
		H(j+1,j) = beta;
		//w = D^(-1).V[j+1]
		biConjugateGradient->solve(diracOperator,V[j+1],w);
		//Gram schimdt
#pragma omp parallel for
		for (int site = 0; site < w.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				f[site][mu] = w[site][mu];
			}
		}
		f.updateHalo();//TODO maybe not needed
		for (unsigned int i = 0; i <= j+1; ++i) {
			std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(V[i],w));
			H(i,j+1) = proj;
#pragma omp parallel for
			for (int site = 0; site < f.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					f[site][mu] -= proj*V[i][site][mu];
				}
			}
			f.updateHalo();//TODO maybe not needed
		}
		//More stable Gram schimdt
		for (unsigned int i = 0; i <= j+1; ++i) {
			std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(V[i],f));
			H(i,j+1) += proj;
#pragma omp parallel for
			for (int site = 0; site < f.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					f[site][mu] -= proj*V[i][site][mu];
				}
			}
			f.updateHalo();//TODO maybe not needed
		}
	}
	//std::vector< std::complex<real_t> > result(steps);
	eigenvalues.resize(steps);
	
	matrix_t eigvec(steps,steps);
#ifdef EIGEN
	Eigen::ComplexEigenSolver<matrix_t> ces(H, true);
	for (unsigned int i = 0; i < steps; ++i) {
		eigenvalues[i] = static_cast<real_t>(1.)/ces.eigenvalues()[i];
	}
	eigvec = ces.eigenvectors();
#endif
#ifdef ARMADILLO
	vector_t eigval(steps);
	arma::eig_gen(eigval, eigvec, H);
	for (unsigned int i = 0; i < steps; ++i) {
		eigenvalues[i] = eigval[i];
	}
#endif

	eigenvectors.resize(steps);
	for (unsigned int i = 0; i < steps; ++i) {
		AlgebraUtils::setToZero(eigenvectors[i]);
		for (unsigned int j = 0; j < steps; ++j) {
#pragma omp parallel for
			for (int site = 0; site < Layout::localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					eigenvectors[i][site][mu] += eigvec.at(j,i)*V[j][site][mu];
				}
			}
		}
		eigenvectors[i].updateHalo();
		
	}

	reduced_dirac_vector_t tmp, tmpe;
	
	//Now we check the convergence
	diracOperator->multiply(tmp,eigenvectors.back());
#pragma omp parallel for
	for (int site = 0; site < Layout::localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			tmpe[site][mu] = eigenvalues.back()*eigenvectors.back()[site][mu];
		}
	}

	std::complex<long_real_t> diffnorm = AlgebraUtils::differenceNorm(tmp,tmpe);

	if (isOutputProcess()) std::cout << "DiracEigenSolver::Convergence precision: " << abs(diffnorm) << std::endl;

	std::reverse(eigenvectors.begin(),eigenvectors.end());
	std::reverse(eigenvalues.begin(),eigenvalues.end());
	
	eigenvalues.erase(eigenvalues.end() - extra_steps, eigenvalues.end());
	eigenvectors.erase(eigenvectors.end() - extra_steps, eigenvectors.end());

	//std::reverse(eigenvectors.begin(),eigenvectors.end());
	//std::reverse(eigenvalues.begin(),eigenvalues.end());
/*#ifdef EIGEN
	Eigen::ComplexEigenSolver<matrix_t> ces(H, true);
	for (unsigned int i = 0; i < steps; ++i) {
		eigenvalues[i] = ces.eigenvalues()[i];
	}
#endif
#ifdef ARMADILLO
	matrix_t eigvec(steps,steps);
	vector_t eigval(steps);
	arma::eig_gen(eigval, eigvec, H);
	for (unsigned int i = 0; i < steps; ++i) {
		eigenvalues[i] = eigval[i];
	}
#endif
	std::sort(eigenvalues.begin(),eigenvalues.end(),mincomparison);

#ifdef EIGEN
	//Now we test the results, looking at the last error
	unsigned int lastvectorindex = 0;
	for (unsigned int i = 0; i < steps; ++i) {
		if (ces.eigenvalues()[i] == eigenvalues[0]) {
			lastvectorindex = i;
			break;
		}
	}

	reduced_dirac_vector_t eigenvector;
#pragma omp parallel for
	for (int site = 0; site < w.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			set_to_zero(eigenvector[site][mu]);
		}
	}

	for (unsigned int i = 0; i < steps; ++i) {
#pragma omp parallel for
		for (int site = 0; site < w.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				eigenvector[site][mu] += V[i][site][mu]*ces.eigenvectors()(lastvectorindex,i);
			}
		}
	}

	eigenvector.updateHalo();

	reduced_dirac_vector_t tmp;

	diracOperator->multiply(tmp, eigenvector);

	std::cout << "Error in norm: " << (tmp[0][0][0] - eigenvalues[n]*eigenvector[0][0][0]) << std::endl;
#endif*/
	//result.erase(result.end() - extra_steps, result.end());
	//return result;
}

void DiracEigenSolver::minimumNonHermitianEigenvalues(DiracOperator* diracOperator, DiracOperator* squareHermitianDiracOperator, std::vector< std::complex<real_t> >& eigenvalues, std::vector<reduced_dirac_vector_t>& eigenvectors/*, Polynomial& map*/, unsigned int n/*, int nmode*/) {
	typedef reduced_dirac_vector_t::Layout Layout;

	unsigned int steps = extra_steps + n;
	//The orthonormal vectors generated by the arnoldi process
	std::vector<reduced_dirac_vector_t> V;
	//Reserve some memory
	V.resize(steps);
	AlgebraUtils::generateRandomVector(V[0]);
	AlgebraUtils::normalize(V[0]);
	reduced_dirac_vector_t w, f, tmp, tmpe;
	if (biConjugateGradient == 0) biConjugateGradient = new BiConjugateGradient();
	biConjugateGradient->setPrecision(epsilon);
	biConjugateGradient->setMaximumSteps(20000);
	//w = D^(-1).V[0]
	diracOperator->setGamma5(true);
	biConjugateGradient->solve(squareHermitianDiracOperator,V[0],tmp);
	diracOperator->multiply(w,tmp);
	AlgebraUtils::gamma5(w);
	std::complex<real_t> alpha = static_cast< std::complex<real_t> >(AlgebraUtils::dot(V[0], w));
	matrix_t H(steps,steps);
	H.zeros();
	//f = w - alpha*V[0]
#pragma omp parallel for
	for (int site = 0; site < w.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			f[site][mu] = w[site][mu] - alpha*V[0][site][mu];
		}
	}
	f.updateHalo();//TODO maybe not needed
	H(0,0) = alpha;
	for (unsigned int j = 0; j < steps - 1; ++j) {
		//beta = norm(f)
		real_t beta = sqrt(AlgebraUtils::squaredNorm(f));
		//V[j] = f/beta
#pragma omp parallel for
		for (int site = 0; site < f.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				V[j+1][site][mu] = f[site][mu]/beta;
			}
		}
		V[j+1].updateHalo();
		//H(j+1,j) = beta
		H(j+1,j) = beta;
		//w = D^(-1).V[j+1]
		biConjugateGradient->solve(squareHermitianDiracOperator,V[j+1],tmp);
		diracOperator->multiply(w,tmp);
		AlgebraUtils::gamma5(w);
		//Gram schimdt
#pragma omp parallel for
		for (int site = 0; site < w.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				f[site][mu] = w[site][mu];
			}
		}
		f.updateHalo();//TODO maybe not needed
		for (unsigned int i = 0; i <= j+1; ++i) {
			std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(V[i],w));
			H(i,j+1) = proj;
#pragma omp parallel for
			for (int site = 0; site < f.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					f[site][mu] -= proj*V[i][site][mu];
				}
			}
			f.updateHalo();//TODO maybe not needed
		}
		//More stable Gram schimdt
		for (unsigned int i = 0; i <= j+1; ++i) {
			std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(V[i],f));
			H(i,j+1) += proj;
#pragma omp parallel for
			for (int site = 0; site < f.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					f[site][mu] -= proj*V[i][site][mu];
				}
			}
			f.updateHalo();//TODO maybe not needed
		}
	}
	//std::vector< std::complex<real_t> > result(steps);
	eigenvalues.resize(steps);
	
	matrix_t eigvec(steps,steps);
#ifdef EIGEN
	Eigen::ComplexEigenSolver<matrix_t> ces(H, true);
	for (unsigned int i = 0; i < steps; ++i) {
		eigenvalues[i] = static_cast<real_t>(1.)/ces.eigenvalues()[i];
	}
	eigvec = ces.eigenvectors();
#endif
#ifdef ARMADILLO
	vector_t eigval(steps);
	arma::eig_gen(eigval, eigvec, H);
	for (unsigned int i = 0; i < steps; ++i) {
		eigenvalues[i] = eigval[i];
	}
#endif

	eigenvectors.resize(steps);
	for (unsigned int i = 0; i < steps; ++i) {
		AlgebraUtils::setToZero(eigenvectors[i]);
		for (unsigned int j = 0; j < steps; ++j) {
#pragma omp parallel for
			for (int site = 0; site < Layout::localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					eigenvectors[i][site][mu] += eigvec.at(j,i)*V[j][site][mu];
				}
			}
		}
		eigenvectors[i].updateHalo();
		
	}
	
	//Now we check the convergence
	diracOperator->setGamma5(false);
	diracOperator->multiply(tmp,eigenvectors.back());
#pragma omp parallel for
	for (int site = 0; site < Layout::localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			tmpe[site][mu] = eigenvalues.back()*eigenvectors.back()[site][mu];
		}
	}

	std::complex<long_real_t> diffnorm = AlgebraUtils::differenceNorm(tmp,tmpe);

	if (isOutputProcess()) std::cout << "DiracEigenSolver::Convergence precision: " << abs(diffnorm) << std::endl;

	std::reverse(eigenvectors.begin(),eigenvectors.end());
	std::reverse(eigenvalues.begin(),eigenvalues.end());
	
	eigenvalues.erase(eigenvalues.end() - extra_steps, eigenvalues.end());
	eigenvectors.erase(eigenvectors.end() - extra_steps, eigenvectors.end());

	//std::reverse(eigenvectors.begin(),eigenvectors.end());
	//std::reverse(eigenvalues.begin(),eigenvalues.end());
/*#ifdef EIGEN
	Eigen::ComplexEigenSolver<matrix_t> ces(H, true);
	for (unsigned int i = 0; i < steps; ++i) {
		eigenvalues[i] = ces.eigenvalues()[i];
	}
#endif
#ifdef ARMADILLO
	matrix_t eigvec(steps,steps);
	vector_t eigval(steps);
	arma::eig_gen(eigval, eigvec, H);
	for (unsigned int i = 0; i < steps; ++i) {
		eigenvalues[i] = eigval[i];
	}
#endif
	std::sort(eigenvalues.begin(),eigenvalues.end(),mincomparison);

#ifdef EIGEN
	//Now we test the results, looking at the last error
	unsigned int lastvectorindex = 0;
	for (unsigned int i = 0; i < steps; ++i) {
		if (ces.eigenvalues()[i] == eigenvalues[0]) {
			lastvectorindex = i;
			break;
		}
	}

	reduced_dirac_vector_t eigenvector;
#pragma omp parallel for
	for (int site = 0; site < w.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			set_to_zero(eigenvector[site][mu]);
		}
	}

	for (unsigned int i = 0; i < steps; ++i) {
#pragma omp parallel for
		for (int site = 0; site < w.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				eigenvector[site][mu] += V[i][site][mu]*ces.eigenvectors()(lastvectorindex,i);
			}
		}
	}

	eigenvector.updateHalo();

	reduced_dirac_vector_t tmp;

	diracOperator->multiply(tmp, eigenvector);

	std::cout << "Error in norm: " << (tmp[0][0][0] - eigenvalues[n]*eigenvector[0][0][0]) << std::endl;
#endif*/
	//result.erase(result.end() - extra_steps, result.end());
	//return result;
}

void DiracEigenSolver::setExtraSteps(unsigned int _extra_steps) {
	extra_steps = _extra_steps;
}

unsigned int DiracEigenSolver::getExtraSteps() const {
	return extra_steps;
}

} /* namespace Update */
