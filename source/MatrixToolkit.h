#ifndef MATRIXTOOLKIT_H
#define MATRIXTOOLKIT_H
#include <complex>
#include <iostream>

namespace matrix_toolkit {

template<typename T> class complex_container { };

template<typename T> class complex_container< std::complex<T> > {
public:
	static std::complex<T> conj(const std::complex<T>& c) {
		return std::complex<T>(real(c),-imag(c));
	}
};

template<> class complex_container< double > {
public:
	static double conj(const double& c) {
		return c;
	}
};

template<typename T, typename U> class tr { };

template<typename T> class tr<T,T> {
	public:
		typedef T type;
};

template<typename T> class tr< std::complex<T> , T> {
	public:
		typedef std::complex<T> type;
};

template<typename T> class tr< T, std::complex<T> > {
	public:
		typedef std::complex<T> type;
};

template<typename T, int NN> 
class Vector {
public:
	Vector() {}
	Vector(const Vector& snd) {
		for (int i = 0; i < NN; ++i) {
			data[i] = snd.data[i];
		}
	}
	Vector& operator=(const Vector& snd) {
		for (int i = 0; i < NN; ++i) {
			data[i] = snd.data[i];
		}
		return *this;
	}

	Vector& operator-() {
		for (int i = 0; i < NN; ++i) {
			data[i] = -data[i];
		}
		return *this;
	}

	Vector operator+(const Vector& snd) const {
		Vector result;
		for (int i = 0; i < NN; ++i) {
			result.data[i] = data[i] + snd.data[i];
		}
		return result;
	}

	Vector operator-(const Vector& snd) const {
		Vector result;
		for (int i = 0; i < NN; ++i) {
			result.data[i] = data[i] - snd.data[i];
		}
		return result;
	}

	Vector& operator+=(const Vector& snd) {
		for (int i = 0; i < NN; ++i) {
			data[i] += snd.data[i];
		}
		return *this;
	}

	Vector& operator-=(const Vector& snd) {
		for (int i = 0; i < NN; ++i) {
			data[i] -= snd.data[i];
		}
		return *this;
	}
	
	template<typename U> Vector operator/(const U& snd) const {
		Vector result;
		for (int i = 0; i < NN; ++i) {
			result.data[i] = data[i]/snd;
		}
		return result;
	}

	T& at(int i) {
		return data[i];
	}

	const T& at(int i) const {
		return data[i];
	}

	T& operator[](int i) {
		return data[i];
	}

	const T& operator[](int i) const {
		return data[i];
	}

	T& operator()(int i) {
		return data[i];
	}

	const T& operator()(int i) const {
		return data[i];
	}
	
private:
	T data[NN];
};

template<typename T, int NN> T vector_dot(const Vector<T, NN>& a, const Vector<T, NN>& b) {
	T result = 0;
	for (int i = 0; i < NN; ++i) {
		result += complex_container<T>::conj(a[i])*b[i];
	}
	return result;
}

template<typename U, typename T, int NN> Vector<T, NN> operator*(const U& a, const Vector<T, NN>& b) {
	Vector<T, NN> result;
	for (int i = 0; i < NN; ++i) {
		result[i] = static_cast<T>(a)*b[i];
	}
	return result;
}

template<typename U, typename T, int NN> Vector<T, NN> operator*(const Vector<T, NN>& a, const U& b) {
	Vector<T, NN> result;
	for (int i = 0; i < NN; ++i) {
		result[i] = a[i]*static_cast<T>(b);
	}
	return result;
}

template<typename T, int NN> void set_to_zero(Vector<T, NN>& a) {
	for (int i = 0; i < NN; ++i) {
		a[i] = 0.;
	}
}

template<typename T, int NN> 
class Matrix {
public:
	Matrix() {}
	Matrix(const Matrix& snd) {
		for (int i = 0; i < NN; ++i) {
			for (int j = 0; j< NN; ++j) {
				data[i][j] = snd.data[i][j];
			}
		}
	}
	Matrix& operator=(const Matrix& snd) {
		for (int i = 0; i < NN; ++i) {
			for (int j = 0; j< NN; ++j) {
				data[i][j] = snd.data[i][j];
			}
		}
		return *this;
	}

	Matrix operator+(const Matrix& snd) const {
		Matrix result;
		for (int i = 0; i < NN; ++i) {
			for (int j = 0; j< NN; ++j) {
				result.data[i][j] = data[i][j] + snd.data[i][j];
			}
		}
		return result;
	}

	Matrix operator-(const Matrix& snd) const {
		Matrix result;
		for (int i = 0; i < NN; ++i) {
			for (int j = 0; j< NN; ++j) {
				result.data[i][j] = data[i][j] - snd.data[i][j];
			}
		}
		return result;
	}

	template<typename U> Matrix< typename tr<U,T>::type , NN> operator*(const Matrix<U,NN>& snd) const {
		Matrix< typename tr<U,T>::type , NN> result;
		for (int i = 0; i < NN; ++i) {
			for (int j = 0; j< NN; ++j) {
				result.data[i][j] = 0.;
				for (int k = 0; k < NN; ++k) {
					result.data[i][j] += data[i][k]*snd[k][j];
				}
			}
		}
		return result;
	}

	Matrix& operator*=(const Matrix& snd) {
		Matrix result;
		for (int i = 0; i < NN; ++i) {
			for (int j = 0; j< NN; ++j) {
				result.data[i][j] = 0.;
				for (int k = 0; k < NN; ++k) {
					result.data[i][j] += data[i][k]*snd.data[k][j];
				}
			}
		}
		*this = result;
		return *this;
	}

	Matrix& operator+=(const Matrix& snd) {
		for (int i = 0; i < NN; ++i) {
			for (int j = 0; j < NN; ++j) {
				data[i][j] += snd.data[i][j];
			}
		}
		return *this;
	}

	Matrix& operator-=(const Matrix& snd) {
		for (int i = 0; i < NN; ++i) {
			for (int j = 0; j < NN; ++j) {
				data[i][j] -= snd.data[i][j];
			}
		}
		return *this;
	}

	Matrix operator-() const {
		Matrix result;
		for (int i = 0; i < NN; ++i) {
			for (int j = 0; j< NN; ++j) {
				result.data[i][j] = -data[i][j];
			}
		}
		return result;
	}

	T& at(int i, int j) {
		return data[i][j];
	}

	const T& at(int i, int j) const {
		return data[i][j];
	}

	T& operator()(int i, int j) {
		return data[i][j];
	}

	const T& operator()(int i, int j) const {
		return data[i][j];
	}

	T* operator[](int i) {
		return data[i];
	}

	const T* operator[](int i) const {
		return data[i];
	}
	
private:
	T data[NN][NN];
};

template<typename T, int NN> T det(const Matrix<T,NN>& ptr) {
    int i,j,k,h,c,s, acs = 1;
	int n = NN;
    Matrix<T, NN> a(NN,NN);
    for (i = 0; i < n;i++) for (j=0; j<n;j++) a[i][j] = ptr[i][j];
    T mi, tmp, det = 1;
    for (int i = 0; i < n; ++i) {
        int j = i;
        if (a[j][i] == 0.) {
            int s = j;
            while (a[s][i] == 0.) {
                s++;
                if (s == n) {
                    return 0;
                }    
            }
            for (int c = 0; c < n; c++) {//sostituzione righe caso 0
                tmp = a[s][c];
                a[s][c] = a[j][c];
                a[j][c] = tmp;
            }
            det = -det;
        }
        mi = a[j][i];
        det *= mi;
        for (int k = i; k < n; k++) a[j][k] /= mi;
        for (j=(i+1); j < n; j++) {
            //salta caso zero
            if(a[j][i] != 0.) {
                mi = a[j][i];
                det *= mi;
                for(int h = i; h < n; h++) a[j][h] = (a[j][h]/mi) - a[i][h];
            }
        }
    }
    if(a[n-1][n-1] == 0.) {
        return 0;
    }
    return det;
}

template<typename T> T det(const Matrix<T,2>& ptr) {
    return ptr.at(0,0)*ptr.at(1,1) - ptr.at(1,0)*ptr.at(0,1);
}

template<typename T, int NN> Matrix<T, NN> htrans(const Matrix<T, NN>& matrix) {
	Matrix<T, NN> result;
	for (int i = 0; i < NN; ++i) {
		for (int j = 0; j< NN; ++j) {
			result.at(i,j) = complex_container<T>::conj(matrix.at(j,i));
		}
	}
	return result;
}

/*template<typename T, int NN> Matrix<std::complex<T>, NN> htrans(const Matrix<std::complex<T>, NN>& matrix) {
	Matrix<std::complex<T>, NN> result;
	for (int i = 0; i < NN; ++i) {
		for (int j = 0; j< NN; ++j) {
			result.at(i,j) = conj(matrix.at(j,i));
		}
	}
	return result;
}*/

template<typename T, int NN> Matrix<T, NN> adj(const Matrix<T, NN>& matrix) {
	Matrix<T, NN> result;
	for (int i = 0; i < NN; ++i) {
		for (int j = 0; j< NN; ++j) {
			result.at(i,j) = complex_container<T>::conj(matrix.at(j,i));
		}
	}
	return result;
}

/*template<typename T, int NN> Matrix<std::complex<T>, NN> adj(const Matrix<std::complex<T>, NN>& matrix) {
	Matrix<std::complex<T>, NN> result;
	for (int i = 0; i < NN; ++i) {
		for (int j = 0; j< NN; ++j) {
			result.at(i,j) = conj(matrix.at(j,i));
		}
	}
	return result;
}*/

template<typename T, int NN> T trace(const Matrix<T, NN>& matrix) {
	T result = 0.;
	for (int i = 0; i < NN; ++i) {
		result += matrix.at(i,i);
	}
	return result;
}

template<typename U, typename T, int NN> Matrix< typename tr<U,T>::type, NN> operator*(const U& a, const Matrix<T, NN>& b) {
	Matrix<typename tr<U,T>::type, NN> result;
	for (int i = 0; i < NN; ++i) {
		for (int j = 0; j < NN; ++j) {
			result[i][j] = a*b[i][j];
		}
	}
	return result;
}

template<typename U, typename T, int NN> Matrix< typename tr<U,T>::type, NN> operator*(const Matrix<U, NN>& a, const T& b) {
	Matrix<typename tr<U,T>::type, NN> result;
	for (int i = 0; i < NN; ++i) {
		for (int j = 0; j < NN; ++j) {
			result[i][j] = a[i][j]*b;
		}
	}
	return result;
}

template<typename T, typename U, int NN> Vector< typename tr<U,T>::type , NN> operator*(const Matrix<T,NN>& matrix, const Vector<U, NN>& snd) const {
	Vector< typename tr<U,T>::type , NN> result;
	for (int i = 0; i < NN; ++i) {
		result[i] = 0.;
		for (int j = 0; j < NN; ++j) {
			result[i] += matrix[i][j]*snd[j];
		}
	}
	return result;
}

template<typename U, typename T, int NN> Matrix< typename tr<U,T>::type, NN> operator/(const Matrix<U, NN>& a, const T& b) {
	Matrix<typename tr<U,T>::type, NN> result;
	for (int i = 0; i < NN; ++i) {
		for (int j = 0; j < NN; ++j) {
			result[i][j] = a[i][j]/b;
		}
	}
	return result;
}

template<typename T, int NN> void set_to_zero(Matrix<T, NN>& a) {
	for (int i = 0; i < NN; ++i) {
		for (int j = 0; j < NN; ++j) {
			a[i][j] = 0.;
		}
	}
}

template<typename T, int NN> T trace(Matrix<T, NN>& a) {
	T result(0.);
	for (int i = 0; i < NN; ++i) {
		result += a[i][i];
	}
	return result;
}

template<typename T, int NN> void set_to_identity(Matrix<T, NN>& a) {
	for (int i = 0; i < NN; ++i) {
		for (int j = 0; j < NN; ++j) {
			if (i != j) a[i][j] = 0.;
			else a[i][j] = 1.;
		}
	}
}

template<typename T> void eigensystem(Vector<T, 2>& eigenvalues, Matrix<T, 2>& eigenvectors, const Matrix<T, 2>& matrix) {
	eigenvalues[0] = 0.5*(matrix.at(0,0) + matrix.at(1,1) - sqrt(matrix.at(0,0)*matrix.at(0,0) + 4.*matrix.at(0,1)*matrix.at(1,0) - 2.*matrix.at(0,0)*matrix.at(1,1) + matrix.at(1,1)*matrix.at(1,1)));
	eigenvalues[1] = 0.5*(matrix.at(0,0) + matrix.at(1,1) + sqrt(matrix.at(0,0)*matrix.at(0,0) + 4.*matrix.at(0,1)*matrix.at(1,0) - 2.*matrix.at(0,0)*matrix.at(1,1) + matrix.at(1,1)*matrix.at(1,1)));
	eigenvectors.at(0,0) = -((-matrix.at(0,0)+matrix.at(1,1)+sqrt(matrix.at(0,0)*matrix.at(0,0)+4.*matrix.at(0,1)*matrix.at(1,0) - 2.*matrix.at(0,0)*matrix.at(1,1)+matrix.at(1,1)*matrix.at(1,1)))/(2.*matrix.at(1,0)));
	eigenvectors.at(0,1) = 1.;
	eigenvectors.at(1,0) = -((-matrix.at(0,0)+matrix.at(1,1)-sqrt(matrix.at(0,0)*matrix.at(0,0)+4.*matrix.at(0,1)*matrix.at(1,0) - 2.*matrix.at(0,0)*matrix.at(1,1)+matrix.at(1,1)*matrix.at(1,1)))/(2.*matrix.at(1,0)));
	eigenvectors.at(1,1) = 1.;
	T norm0 = sqrt(conj(eigenvectors.at(0,0))*eigenvectors.at(0,0) + conj(eigenvectors.at(0,1))*eigenvectors.at(0,1));
	T norm1 = sqrt(conj(eigenvectors.at(1,0))*eigenvectors.at(1,0) + conj(eigenvectors.at(1,1))*eigenvectors.at(1,1));
	eigenvectors.at(0,0) /= norm0;
	eigenvectors.at(0,1) /= norm0;
	eigenvectors.at(1,0) /= norm1;
	eigenvectors.at(1,1) /= norm1;
}

template<typename T> Matrix< std::complex<T>, 2> matrix_exp(const Matrix< std::complex<T>, 2>& matrix) {
	Matrix< std::complex<T>, 2> result;
	T n1(imag(matrix(0,0)));
	T n2(imag(matrix(0,1)));
	T n3(real(matrix(0,1)));
	T norm(sqrt(n1*n1 + n2*n2 + n3*n3));
	if (norm < 0.0000001) {
		result(0,0) = std::complex<T>(1-norm*norm/2.,n1*(1-norm*norm/6.));
		result(0,1) = std::complex<T>(n3*(1-norm*norm/6.),n2*(1-norm*norm/6.));
		result(1,0) = std::complex<T>(-n3*(1-norm*norm/6.),n2*(1-norm*norm/6.));
		result(1,1) = std::complex<T>(1-norm*norm/2.,-n1*(1-norm*norm/6.));
	}
	else {
		result(0,0) = std::complex<T>(cos(norm),n1*sin(norm)/norm);
		result(0,1) = std::complex<T>(n3*sin(norm)/norm,n2*sin(norm)/norm);
		result(1,0) = std::complex<T>(-n3*sin(norm)/norm,n2*sin(norm)/norm);
		result(1,1) = std::complex<T>(cos(norm),-n1*sin(norm)/norm);
	}
	return result;
}

template<typename T> 
class Matrix<T, -1> {
public:
	Matrix() : data(0), m(0), n(0) {}
	Matrix(int _m, int _n) : m(_m) , n(_n) {
		data = new T*[m];
		for (int i = 0; i < m; ++i) data[i] = new T[n];
	}
	Matrix(const Matrix& snd) {
		if (m != snd.m || n != snd.n) {
			if (data) {
				for (int i = 0; i < m; ++i) delete[] data[i];
				delete[] data;
			}
			m = snd.m;
			n = snd.n;
			data = new T*[m];
			for (int i = 0; i < m; ++i) data[i] = new T[n];
		}
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				data[i][j] = snd.data[i][j];
			}
		}
	}
	Matrix& operator=(const Matrix& snd) {
		if (m != snd.m || n != snd.n) {
			if (data) {
				for (int i = 0; i < m; ++i) delete[] data[i];
				delete[] data;
			}
			m = snd.m;
			n = snd.n;
			data = new T*[m];
			for (int i = 0; i < m; ++i) data[i] = new T[n];
		}
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				data[i][j] = snd.data[i][j];
			}
		}
		return *this;
	}
	~Matrix() {
		if (data) {
			for (int i = 0; i < m; ++i) delete[] data[i];
			delete[] data;
		}
	}

	Matrix operator+(const Matrix& snd) const {
		Matrix result;
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j< n; ++j) {
				result.data[i][j] = data[i][j] + snd.data[i][j];
			}
		}
		return result;
	}

	Matrix operator-(const Matrix& snd) const {
		Matrix result;
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j< n; ++j) {
				result.data[i][j] = data[i][j] - snd.data[i][j];
			}
		}
		return result;
	}

	Matrix operator*(const Matrix& snd) const {
		Matrix result;
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j< n; ++j) {
				result.data[i][j] = 0.;
				for (int k = 0; k < n; ++k) {
					result.data[i][j] += data[i][k]*snd.data[k][j];
				}
			}
		}
		return result;
	}

	Matrix& operator+=(const Matrix& snd) {
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				data[i][j] += snd.data[i][j];
			}
		}
		return *this;
	}

	Matrix& operator-=(const Matrix& snd) {
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				data[i][j] -= snd.data[i][j];
			}
		}
		return *this;
	}

	T& at(int i, int j) {
		return data[i][j];
	}

	const T& at(int i, int j) const {
		return data[i][j];
	}

	T& operator()(int i, int j) {
		return data[i][j];
	}

	const T& operator()(int i, int j) const {
		return data[i][j];
	}

	T* operator[](int i) {
		return data[i];
	}

	const T* operator[](int i) const {
		return data[i];
	}

	int nRows() const {
		return m;
	}

	int nCols() const {
		return n;
	}
	
private:
	T** data;
	int m;
	int n;
};

template<typename T> void set_to_zero(Matrix<T, -1>& a) {
	for (int  i = 0; i < a.nRows(); ++i) {
		for (int  j = 0; j < a.nCols(); ++j) {
			a[i][j] = 0.;
		}
	}
}

template<typename T> Matrix<T, -1> inverse(const Matrix<T, -1>& ptr) {
	int n = ptr.nRows();
    Matrix<T, -1> temp(n,n);
    Matrix<T, -1> a(n,n);
    if (ptr.nRows() != ptr.nCols()) {
        return temp;
    }
    for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			a[i][j] = ptr[i][j];
		}
	}
    T mi, tmp;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i==j) temp[i][j] = 1;
            else temp[i][j] = 0;
        }
    }
    for (int i = 0; i < n; ++i) {
        int j=i;
        if (a[j][i] == 0.) {
            int s = j;
            while (a[s][i] != 0.) {
                s++;
                if (s == n) {
              		//Non invertible
                    return temp;
                }    
            }
            for (int c = 0; c < n; ++c) {//sostituzione righe caso 0
                tmp = a[s][c];
                a[s][c] = a[j][c];
                a[j][c] = tmp;
                tmp = temp[s][c];
                temp[s][c] = temp[j][c];
                temp[j][c] = tmp;
            }
        }
        mi = a[j][i];
        
        for (int k = i; k < n; ++k) a[j][k] /= mi;
        for (int c = 0; c < n; ++c) temp[i][c] /= mi;
        
        for (j=(i+1); j < n; ++j) {
            //salta caso zero
            if(a[j][i] != 0.) {
                mi = a[j][i];
                for (int h = i; h < n; ++h) a[j][h] = (a[j][h]/mi) - a[i][h];
                for (int c = 0; c < n; ++c) temp[j][c] = (temp[j][c]/mi)- temp[i][c];
            }
        }
    }
    if(a[n-1][n-1] == 0.) {
        return temp;
    }
    for (int i = n-1; i >= 0; --i) {
        int j = i - 1;
        mi = a[i][i];
        for (int k = 0; k < n; ++k) a[i][k] /= mi;
        for (int c = 0; c < n; ++c) temp[i][c] /= mi;
        
        for (; j >= 0; j--) {
            mi = a[j][i];
            for (int h = i; h >= 0; --h) a[j][h] = (a[j][h]/mi) - a[i][h];
            for (int c = 0; c < n; ++c) temp[j][c] = (temp[j][c]/mi)- temp[i][c];
        }
    }
    return temp;
}

template<typename T, int NN> std::ostream& operator<<(std::ostream& os, const Matrix<T, NN>& m) {
	os << std::endl;
	for (int i = 0; i < NN; ++i) {
		for (int j = 0; j < NN; ++j) os << m[i][j]<< " ";
		os << std::endl;
	}
	os << std::endl;
	return os;
}

template<typename T, int NN> std::ostream& operator<<(std::ostream& os, const Vector<T, NN>& v) {
	os << std::endl;
	for (int j = 0; j < NN; ++j) os << v[j]<< " ";
	os << std::endl;
	return os;
}

}

#endif
