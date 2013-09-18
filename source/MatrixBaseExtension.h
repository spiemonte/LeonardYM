inline Scalar at(uint i, uint j) const { return this->operator()(i,j); }
inline Scalar& at(uint i, uint j) { return this->operator()(i,j); }
inline Scalar at(uint i) const { return this->operator[](i); }
inline Scalar& at(uint i) { return this->operator[](i); }
inline void zeros() {
	for (unsigned int i = 0; i < this->rows(); ++i) {
		for (unsigned int j = 0; j < this->cols(); ++j) {
			this->operator()(i,j) = 0;
		}
	}
}

