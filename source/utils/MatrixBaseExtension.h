typedef unsigned int uint;
inline Scalar at(uint i, uint j) const { return this->operator()(i,j); }
inline Scalar& at(uint i, uint j) { return this->operator()(i,j); }
inline Scalar at(uint i) const { return this->operator[](i); }
inline Scalar& at(uint i) { return this->operator[](i); }
inline void zeros() {
	for (int i = 0; i < this->rows(); ++i) {
		for (int j = 0; j < this->cols(); ++j) {
			this->operator()(i,j) = 0;
		}
	}
}

