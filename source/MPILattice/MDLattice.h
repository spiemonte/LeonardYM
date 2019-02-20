#ifndef MDLATTICE_H
#define MDLATTICE_H
#include "Lattice.h"
#include "MatrixTypedef.h"


namespace Lattice {

template<typename T, typename TStored, typename TLayout, int M, int N> class ComplexVectorLattice {
public:
	ComplexVectorLattice() { }
	ComplexVectorLattice(const ComplexVectorLattice& copy) {
		for (int k = 0; k < M*N; ++k) {
			real_part[k] = copy.real_part[k];
			imag_part[k] = copy.imag_part[k];
		}
	}
	ComplexVectorLattice(const Lattice<TStored, TLayout>& copy) {
		for (int mu = 0; mu < M; ++mu) {
			for (int c = 0; c < N; ++c) {
#pragma omp parallel for
				for (int site = 0; site < TLayout::completesize; ++site) {
					real_part[c*M + mu][site] = copy[site][mu][c].real();
					imag_part[c*M + mu][site] = copy[site][mu][c].imag();
				}
			}
		}
	}

	ComplexVectorLattice& operator=(const Lattice<TStored, TLayout>& copy) {
		for (int mu = 0; mu < M; ++mu) {
			for (int c = 0; c < N; ++c) {
#pragma omp parallel for
				for (int site = 0; site < TLayout::completesize; ++site) {
					real_part[c*M + mu][site] = copy[site][mu][c].real();
					imag_part[c*M + mu][site] = copy[site][mu][c].imag();
				}
			}
		}

		return *this;
	}

	template<typename ULayout> void copy_to(Lattice<TStored,ULayout>& copy) const {
		for (int mu = 0; mu < M; ++mu) {
			for (int c = 0; c < N; ++c) {
				Lattice<T,ULayout> rtmp = real_part[c*M + mu];
				Lattice<T,ULayout> itmp = imag_part[c*M + mu];
#pragma omp parallel for
				for (int site = 0; site < TLayout::completesize; ++site) {
					copy[site][mu][c] = std::complex<Update::real_t>(rtmp[site], itmp[site]);
				}
			}
		}
	}

	void copy_to(Lattice<TStored, TLayout>& copy) const {
		for (int mu = 0; mu < M; ++mu) {
			for (int c = 0; c < N; ++c) {
#pragma omp parallel for
				for (int site = 0; site < TLayout::completesize; ++site) {
					copy[site][mu][c] = std::complex<Update::real_t>(real_part[c*M + mu][site], imag_part[c*M + mu][site]);
				}
			}
		}
	}

	void updateHalo() {
		for (int k = 0; k < M*N; ++k) {
			real_part[k].communicateHalo();
			imag_part[k].communicateHalo();
		}
		for (int k = 0; k < M*N; ++k) {
			real_part[k].waitHalo();
			imag_part[k].waitHalo();
		}
	}

	Lattice<T, TLayout> real_part[M*N];
	Lattice<T, TLayout> imag_part[M*N];

	bool operator==(const ComplexVectorLattice& snd) const {
		for (int k = 0; k < M*N; ++k) {
			if (real_part[k] != snd.real_part[k] || imag_part[k] != snd.imag_part[k]) {
				return false;
			}
		}
		return true;
	}
};

template<typename T, typename TStored, typename TLayout, int M, int N, int P> class ComplexMatrixLattice {
public:
	ComplexMatrixLattice() { }
	ComplexMatrixLattice(const ComplexMatrixLattice& copy) {
		for (int k = 0; k < M*N*P; ++k) {
			real_part[k] = copy.real_part[k];
			imag_part[k] = copy.imag_part[k];
		}
	}
	ComplexMatrixLattice(const Lattice<TStored, TLayout>& copy) {
		for (int mu = 0; mu < M; ++mu) {
			for (int c = 0; c < N; ++c) {
				for (int d = 0; d < P; ++d) {
#pragma omp parallel for
					for (int site = 0; site < TLayout::completesize; ++site) {
						real_part[d*N*M + c*M + mu][site] = copy[site][mu].at(c,d).real();
						imag_part[d*N*M + c*M + mu][site] = copy[site][mu].at(c,d).imag();
					}
				}
			}
		}
	}

	ComplexMatrixLattice& operator=(const Lattice<TStored, TLayout>& copy) {
		for (int mu = 0; mu < M; ++mu) {
			for (int c = 0; c < N; ++c) {
				for (int d = 0; d < P; ++d) {
#pragma omp parallel for
					for (int site = 0; site < TLayout::completesize; ++site) {
						real_part[d*N*M + c*M + mu][site] = copy[site][mu].at(c,d).real();
						imag_part[d*N*M + c*M + mu][site] = copy[site][mu].at(c,d).imag();
					}
				}
			}
		}

		return *this;
	}

	template<typename ULayout> void copy_to(Lattice<TStored,ULayout>& copy) const {
		for (int mu = 0; mu < M; ++mu) {
			for (int c = 0; c < N; ++c) {
				for (int d = 0; d < P; ++d) {
					Lattice<T,ULayout> rtmp = real_part[d*N*M + c*M + mu];
					Lattice<T,ULayout> itmp = imag_part[d*N*M + c*M + mu];
#pragma omp parallel for
					for (int site = 0; site < TLayout::completesize; ++site) {
						copy[site][mu].at(c,d) = std::complex<Update::real_t>(rtmp[site], itmp[site]);
					}
				}
			}
		}
	}

	void copy_to(Lattice<TStored, TLayout>& copy) const {
		for (int mu = 0; mu < M; ++mu) {
			for (int c = 0; c < N; ++c) {
				for (int d = 0; d < P; ++d) {
#pragma omp parallel for
					for (int site = 0; site < TLayout::completesize; ++site) {
						copy[site][mu].at(c,d) = std::complex<Update::real_t>(real_part[d*N*M + c*M + mu][site], imag_part[d*N*M + c*M + mu][site]);
					}
				}
			}
		}
	}

	void updateHalo() {
		for (int k = 0; k < M*N*P; ++k) {
			real_part[k].communicateHalo();
			imag_part[k].communicateHalo();
		}
		for (int k = 0; k < M*N*P; ++k) {
			real_part[k].waitHalo();
			imag_part[k].waitHalo();
		}
	}

	Lattice<T, TLayout> real_part[M*N*P];
	Lattice<T, TLayout> imag_part[M*N*P];
};

template<typename T, typename TStored, typename TLayout, int M, int N, int P> class RealMatrixLattice {
public:
	RealMatrixLattice() { }
	RealMatrixLattice(const RealMatrixLattice& copy) {
		for (int k = 0; k < M*N*P; ++k) {
			real_part[k] = copy.real_part[k];
		}
	}
	RealMatrixLattice(const Lattice<TStored, TLayout>& copy) {
		for (int mu = 0; mu < M; ++mu) {
			for (int c = 0; c < N; ++c) {
				for (int d = 0; d < P; ++d) {
#pragma omp parallel for
					for (int site = 0; site < TLayout::completesize; ++site) {
						real_part[d*N*M + c*M + mu][site] = copy[site][mu].at(c,d);
					}
				}
			}
		}
	}

	RealMatrixLattice& operator=(const Lattice<TStored, TLayout>& copy) {
		for (int mu = 0; mu < M; ++mu) {
			for (int c = 0; c < N; ++c) {
				for (int d = 0; d < P; ++d) {
#pragma omp parallel for
					for (int site = 0; site < TLayout::completesize; ++site) {
						real_part[d*N*M + c*M + mu][site] = copy[site][mu].at(c,d);
					}
				}
			}
		}

		return *this;
	}

	template<typename ULayout> void copy_to(Lattice<TStored,ULayout>& copy) const {
		for (int mu = 0; mu < M; ++mu) {
			for (int c = 0; c < N; ++c) {
				for (int d = 0; d < P; ++d) {
					Lattice<T,ULayout> rtmp = real_part[d*N*M + c*M + mu];
#pragma omp parallel for
					for (int site = 0; site < TLayout::completesize; ++site) {
						copy[site][mu].at(c,d) = rtmp[site];
					}
				}
			}
		}
	}

	void copy_to(Lattice<TStored, TLayout>& copy) const {
		for (int mu = 0; mu < M; ++mu) {
			for (int c = 0; c < N; ++c) {
				for (int d = 0; d < P; ++d) {
#pragma omp parallel for
					for (int site = 0; site < TLayout::completesize; ++site) {
						copy[site][mu].at(c,d) = real_part[d*N*M + c*M + mu][site];
					}
				}
			}
		}
	}

	void updateHalo() {
		for (int k = 0; k < M*N*P; ++k) {
			real_part[k].communicateHalo();
		}
		for (int k = 0; k < M*N*P; ++k) {
			real_part[k].waitHalo();
		}
	}

	Lattice<T, TLayout> real_part[M*N*P];
};

}

#endif
