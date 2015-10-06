#ifndef MPITYPE_H
#define MPITYPE_H

#ifdef ENABLE_MPI
#include <mpi.h>
#include "../MatrixTypedef.h"

template<typename T> class MpiType { };

template<> class MpiType<short int> {
	public:
		static const int size = sizeof(short int);
		static MPI_Datatype type;
};

template<> class MpiType<int> {
	public:
		static const int size = sizeof(int);
		static MPI_Datatype type;
};

template<> class MpiType<int[4]> {
	public:
		static const int size = sizeof(int);
		static MPI_Datatype type;
};


template<> class MpiType<double> {
	public:
		static const int size = 8;
		static MPI_Datatype type;
};

template<> class MpiType<double[4]> {
	public:
		static const int size = 8;
		static MPI_Datatype type;
};

template<> class MpiType<Update::FundamentalGroup[4]> {
	public:
		static const int size = 8;
		static MPI_Datatype type;
};

template<> class MpiType<Update::AdjointGroup[4]> {
	public:
		static const int size = 8;
		static MPI_Datatype type;
};

template<> class MpiType<Update::FundamentalGroup[6]> {
	public:
		static const int size = 8;
		static MPI_Datatype type;
};

template<> class MpiType<Update::AdjointGroup[6]> {
	public:
		static const int size = 8;
		static MPI_Datatype type;
};

template<> class MpiType<Update::FundamentalVector[4]> {
	public:
		static const int size = 8;
		static MPI_Datatype type;
};

template<> class MpiType<Update::AdjointVector[4]> {
	public:
		static const int size = 8;
		static MPI_Datatype type;
};

template<> class MpiType<Update::single_FundamentalVector[4]> {
	public:
		static const int size = sizeof(float);
		static MPI_Datatype type;
};

template<> class MpiType<Update::single_AdjointVector[4]> {
	public:
		static const int size = sizeof(float);
		static MPI_Datatype type;
};

template<> class MpiType<Update::FundamentalGroup> {
	public:
		static const int size = 8;
		static MPI_Datatype type;
};

template<> class MpiType<Update::AdjointGroup> {
	public:
		static const int size = 8;
		static MPI_Datatype type;
};

template<> class MpiType<Update::FermionicForceMatrix[4]> {
	public:
		static const int size = 8;
		static MPI_Datatype type;
};

#endif

#endif
