#ifndef MULTIGRIDVECTOR_H
#define MULTIGRIDVECTOR_H
#include "MultiGridVectorLayout.h"

namespace Update {

template<typename T,typename TLayout> class MultiGridVector {
public:
	MultiGridVector() : data(new T[TLayout::size]) {
		++TLayout::instances;
	}
	MultiGridVector(const MultiGridVector& snd) : data(new T[TLayout::size]) {
		++TLayout::instances;
		for (int i = 0; i < Layout::size; ++i) data[i] = snd.data[i];
	}
	~MultiGridVector() {
		--TLayout::instances;
		delete[] data;
	}

	MultiGridVector& operator=(const MultiGridVector& snd) {
		delete[] data;
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

}

#endif

