#ifndef LATTICE_H
#define LATTICE_H
#include "Site.h"
#include "MPIType.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#include <vector>
#include <fstream>
#include "utils/ToString.h"
#endif
#include <iostream>
#include <typeinfo>

//#define DEBUG_MEMORY_ALLOCATION

namespace Lattice {

template<typename T, typename TLayout> class Lattice {
	public:
		Lattice() : localsize(TLayout::localsize), completesize(TLayout::completesize), sharedsize(TLayout::sharedsize) {
#ifdef DEBUG_MEMORY_ALLOCATION
			try {
#endif
				localdata = new T[TLayout::completesize];
#ifdef DEBUG_MEMORY_ALLOCATION
				++allocationCounter;
				if (TLayout::this_processor == 0) std::cout << "Memory allocated in lattice constructor!" << std::endl;
				if (TLayout::this_processor == 0) std::cout << "Number of allocations: " << allocationCounter << std::endl;
				if (TLayout::this_processor == 0) std::cout << " for type: " << typeid(T).name() << std::endl;
			}
			catch (const std::bad_alloc&) {
				std::cout << "Memory error in lattice constructor!" << std::endl;
				std::cout << "Number of allocations: " << allocationCounter << std::endl;
				std::cout << " for type: " << typeid(T).name() << std::endl;
				exit(7);
			}
#endif
		}
		~Lattice() {
			delete[] localdata;
#ifdef DEBUG_MEMORY_ALLOCATION
			--allocationCounter;
#endif
		}
		Lattice(const Lattice& copy) : localsize(TLayout::localsize), completesize(TLayout::completesize), sharedsize(TLayout::sharedsize) {
#ifdef DEBUG_MEMORY_ALLOCATION
			try {
#endif
				localdata = new T[TLayout::completesize];
#ifdef DEBUG_MEMORY_ALLOCATION
				++allocationCounter;
				if (TLayout::this_processor == 0) std::cout << "Memory allocated in lattice copy constructor!" << std::endl;
				if (TLayout::this_processor == 0) std::cout << "Number of allocations: " << allocationCounter << std::endl;
				if (TLayout::this_processor == 0) std::cout << " for type: " << typeid(T).name() << std::endl;
			}
			catch (const std::bad_alloc&) {
				std::cout << "Memory error in lattice copy constructor!" << std::endl;
				std::cout << "Number of allocations: " << allocationCounter << std::endl;
				std::cout << " for type: " << typeid(T).name() << std::endl;
				exit(7);
			}
#endif
			memcpy(localdata, copy.localdata, TLayout::completesize*sizeof(T));
		}
		
		Lattice& operator=(const Lattice& copy) {
			memcpy(localdata, copy.localdata, TLayout::completesize*sizeof(T));
			return *this;
		}
		
		template<typename ULayout> Lattice& operator=(const Lattice<T, ULayout>& copy) {
			int* map = layout.getMapIndex(copy.getLayout());
#pragma omp parallel for
			for (int site = 0; site < layout.localsize; ++site) {
				memcpy(&localdata[site], &copy[map[site]], sizeof(T));
			}
			updateHalo();
			return *this;
		}
		
		template<typename ULayout> Lattice(const Lattice<T, ULayout>& copy) : localsize(TLayout::localsize), completesize(TLayout::completesize), sharedsize(TLayout::sharedsize) {
#ifdef DEBUG_MEMORY_ALLOCATION
			try {
#endif
				localdata = new T[TLayout::completesize];
#ifdef DEBUG_MEMORY_ALLOCATION
				++allocationCounter;
				if (TLayout::this_processor == 0) std::cout << "Memory allocated in lattice template copy constructor!" << std::endl;
				if (TLayout::this_processor == 0) std::cout << "Number of allocations: " << allocationCounter << std::endl;
				if (TLayout::this_processor == 0) std::cout << " for type: " << typeid(T).name() << std::endl;
			}
			catch (const std::bad_alloc&) {
				std::cout << "Memory error in lattice template copy constructor!" << std::endl;
				std::cout << "Number of allocations: " << allocationCounter << std::endl;
				std::cout << " for type: " << typeid(T).name() << std::endl;
				exit(7);
			}
#endif
			int* map = layout.getMapIndex(copy.getLayout());
#pragma omp parallel for
			for (int site = 0; site < layout.localsize; ++site) {
				memcpy(&localdata[site], &copy[map[site]], sizeof(T));
			}
			updateHalo();
		}
		
		typedef TLayout Layout;
		typedef T TData;
		
		void updateHalo() {
			communicateHalo();
			waitHalo();
		}

		void communicateHalo() {
#ifdef ENABLE_MPI
			int packsize = sizeof(T)/MpiType<T>::size;
			
			//Update halo!
			//First we initialize the send
			for (int i = 0; i < layout.numberChunks; ++i) {
				if (layout.latticeChunks[i].owner == layout.this_processor) {
					for (unsigned int j = 0; j < layout.latticeChunks[i].sharers.size(); ++j) {
						void* address = (void*)(&localdata[layout.latticeChunks[i].offset]);
						int dimension = layout.latticeChunks[i].size;
						int destination = layout.latticeChunks[i].sharers[j];
						int tag = layout.latticeChunks[i].tags[j];
						MPI_Request* send_request = new MPI_Request;
						MPI_Isend(address,packsize*dimension,MpiType<T>::type,destination,tag,MPI_COMM_WORLD,send_request);
						sendRequests.push_back(send_request);
						//outputfile << "Processor " << layout.this_processor << " is sending data with tag " << tag  << " to processor " << destination << " reading from offset " << layout.latticeChunks[i].offset << std::endl;
					}
				}
			}
			//Then we receive
			for (int i = 0; i < layout.numberChunks; ++i) {
				if (layout.latticeChunks[i].owner != layout.this_processor) {
					for (unsigned int j = 0; j < layout.latticeChunks[i].sharers.size(); ++j) {
						if (layout.latticeChunks[i].sharers[j] == layout.this_processor) {
							void* address = (void*)(&localdata[layout.latticeChunks[i].offset]);
							int dimension = layout.latticeChunks[i].size;
							int source = layout.latticeChunks[i].owner;
							int tag = layout.latticeChunks[i].tags[j];
							MPI_Request* recv_request = new MPI_Request;
							MPI_Irecv(address,packsize*dimension,MpiType<T>::type,source,tag,MPI_COMM_WORLD,recv_request);
							recvRequests.push_back(recv_request);
							//outputfile << "Processor " << layout.this_processor << " is receving data with tag " << tag << " from processor " << source << " writing to offset " << layout.latticeChunks[i].offset << std::endl;
						}
					}
				}
			}
#endif
		}

		void waitHalo() {
#ifdef ENABLE_MPI
			//Then we wait
			MPI_Status status;
			for (unsigned int i = 0; i < sendRequests.size(); ++i) MPI_Wait(sendRequests[i],&status);
			for (unsigned int i = 0; i < recvRequests.size(); ++i) MPI_Wait(recvRequests[i],&status);
			
			//delete everything
			for (unsigned int i = 0; i < sendRequests.size(); ++i) delete sendRequests[i];
			for (unsigned int i = 0; i < recvRequests.size(); ++i) delete recvRequests[i];

			sendRequests.clear();
			recvRequests.clear();
#endif
		}
		
		T& operator[](unsigned int index) {
			return localdata[index];
		}
		
		const T& operator[](unsigned int index) const {
			return localdata[index];
		}
		
		static int sup(unsigned int site, unsigned int mu) {
			return TLayout::sup_table[site][mu];
		}
		
		static int sdn(unsigned int site, unsigned int mu) {
			return TLayout::down_table[site][mu];
		}
		
		TLayout& getLayout() {
			return layout;
		}
		
		const TLayout& getLayout() const {
			return layout;
		}

		const T* getRawData() const {
			return localdata;
		}

		T* getRawData() {
			return localdata;
		}
		
	private:
		T* localdata;
		TLayout layout;
#ifdef ENABLE_MPI
		std::vector<MPI_Request*> sendRequests;
		std::vector<MPI_Request*> recvRequests;
#endif
				
	public:
		const int localsize;
		const int completesize;
		const int sharedsize;
#ifdef DEBUG_MEMORY_ALLOCATION
		static int allocationCounter;
#endif
};

#ifdef DEBUG_MEMORY_ALLOCATION
template<typename T, typename TLayout> int Lattice<T,TLayout>::allocationCounter = 0;
#endif

}

#endif
