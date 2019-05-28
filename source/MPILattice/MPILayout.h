#ifndef MPILAYOUT_H
#define MPILAYOUT_H
#include "Site.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include <iostream>
#include <algorithm>
#include <map>
#include <fstream>
#include <string>
#include <utility>
#include <rpc/xdr.h>
#include "utils/ToString.h"
#include "LatticeChunk.h"

namespace Lattice {

typedef int vect4[4];

template<typename Stencil> class MpiLayout {
	public:
		//The grid of processor
		static int pgrid_t;
		static int pgrid_x;
		static int pgrid_y;
		static int pgrid_z;
		//The global size
		static int glob_t;
		static int glob_x;
		static int glob_y;
		static int glob_z;
		static int glob[4];
		//The local size
		static int loc_t;
		static int loc_x;
		static int loc_y;
		static int loc_z;
		
		static int numberProcessors;
		static int globalVolume;
		static int glob_spatial_volume;
		
		static void initialize() {		
			numberProcessors = pgrid_x*pgrid_y*pgrid_z*pgrid_t;
			globalVolume = glob_x*glob_y*glob_z*glob_t;
			glob_spatial_volume = glob_x*glob_y*glob_z;

			glob[0] = glob_x;
			glob[1] = glob_y;
			glob[2] = glob_z;
			glob[3] = glob_t;
#ifdef ENABLE_MPI
			int sizeCommunicator = 0;
			MPI_Comm_size(MPI_COMM_WORLD, &sizeCommunicator);
			
			if (sizeCommunicator != numberProcessors) {
				std::cout << "Number of processors " << numberProcessors << " is different from the size of the communicator " << sizeCommunicator << std::endl;
				exit(1);
			}
			
			
			//The processor number that is given to every site
			int* processorNumberLattice = new int[globalVolume];
			//The Site to every integer index
			Site* globalMapSite = new Site[globalVolume];
			loc_x = glob_x / pgrid_x;
			loc_y = glob_y / pgrid_y;
			loc_z = glob_z / pgrid_z;
			loc_t = glob_t / pgrid_t;
			
			//Check that everything is ok
			if ((glob_x % pgrid_x) != 0) {
				std::cout << "The global size " << glob_x << " is not divisible in the grid size " << pgrid_x << std::endl;
				exit(3);
			}
			if ((glob_y % pgrid_y) != 0) {
				std::cout << "The global size " << glob_y << " is not divisible in the grid size " << pgrid_y << std::endl;
				exit(3);
			}
			if ((glob_z % pgrid_z) != 0) {
				std::cout << "The global size " << glob_z << " is not divisible in the grid size " << pgrid_z << std::endl;
				exit(3);
			}
			if ((glob_t % pgrid_t) != 0) {
				std::cout << "The global size " << glob_t << " is not divisible in the grid size " << pgrid_t << std::endl;
				exit(3);
			}
			
//			for (int site = 0; site < globalVolume; ++site) processorNumberLattice[site] = -1;
	
			//First we partition the lattice in processor and we set the coordinate
#pragma omp parallel for
			for (int x = 0; x < glob_x; ++x) {
				for (int y = 0; y < glob_y; ++y) {
					for (int z = 0; z < glob_z; ++z) {
						for (int t = 0; t < glob_t; ++t) {
							int index = getGlobalCoordinate(x, y, z, t);
							int processor = x/loc_x;
							processor = pgrid_y*processor + y/loc_y;
							processor = pgrid_z*processor + z/loc_z;
							processor = pgrid_t*processor + t/loc_t;
							processorNumberLattice[index] = processor;
							globalMapSite[index] = Site(x,y,z,t);
						}
					}
				}
			}

//#pragma omp parallel for
//			for (int site = 0; site < globalVolume; ++site) if (processorNumberLattice[site] == -1 || processorNumberLattice[site] >= numberProcessors) std::cout << "Fatal Error!" << std::endl;
			
			//For every site we store the processor that want to read it, if it is not the same
			std::vector<int>* exchangeTable = new std::vector<int>[globalVolume];
			for (int x = 0; x < glob_x; ++x) {
				for (int y = 0; y < glob_y; ++y) {
					for (int z = 0; z < glob_z; ++z) {
						for (int t = 0; t < glob_t; ++t) {
							for (unsigned int i = 0; i < Stencil::neighbourSites.size(); ++i) {
								Site site(x,y,z,t);
								int index = getGlobalCoordinate(site);
								int shiftedIndex = getGlobalCoordinate(site+Stencil::neighbourSites[i]);
								if (index == shiftedIndex) std::cout << "Fatal error in exchangeTable calculations!" << std::endl;
								if (processorNumberLattice[index] != processorNumberLattice[shiftedIndex]) exchangeTable[shiftedIndex].push_back(processorNumberLattice[index]);
							}
						}
					}
				}
			}

			//We sort and clean duplicates
#pragma omp parallel for
			for (int site = 0; site < globalVolume; ++site) {
				std::sort(exchangeTable[site].begin(), exchangeTable[site].end());
				exchangeTable[site].erase(std::unique(exchangeTable[site].begin(), exchangeTable[site].end()), exchangeTable[site].end());
			}

			/*int ll = 0;
			for (int site = 0; site < globalVolume; ++site) {
				if (processorNumberLattice[site] == this_processor) ++ll;
			}
			if (ll != loc_x*loc_y*loc_z*loc_t) std::cout << "Localsize mismatch" << std::endl;*/
			
			//Now we get the id; to every site it is associated an unique string depending on who would like to read it
			int* idLattice = new int[globalVolume];
			std::vector<std::string> ids;
			for (int site = 0; site < globalVolume; ++site) {
				std::string id = Update::toString(processorNumberLattice[site]);

				for (unsigned int i = 0; i < exchangeTable[site].size(); ++i) {
					id = id + "." + Update::toString(exchangeTable[site][i]);
				}
				
				//Differences also for the size of sharers
				id = id + "." + Update::toString(exchangeTable[site].size());
				
				bool found = false;
				for (unsigned int i = 0; i < ids.size(); ++i) {
					if (ids[i] == id) {
						idLattice[site] = i;
						found = true;
					}
				}
				if (!found) {
					idLattice[site] = ids.size();
					ids.push_back(id);
				}
			}
			
			//Now we divide the lattice in chunk, every chunk has the same owner and the same sharers
			numberChunks = ids.size();
			latticeChunks = std::vector<LatticeChunk>(numberChunks);
			for (int i = 0; i < numberChunks; ++i) {
				latticeChunks[i].size = 0;
				latticeChunks[i].id = i;
				latticeChunks[i].owner = -1;
			}

			//Ensure
			/*if (numberChunks > MPI_TAG_UB) {
				std::cout << "MPI tag maximum number " << MPI_TAG_UB << " is too small for " << numberChunks << " chunks! " << std::endl;
				exit(35);
			}*/
			
			//We take the number of sites in every chunk
			for (int site = 0; site < globalVolume; ++site) {
				++latticeChunks[idLattice[site]].size;
				if (latticeChunks[idLattice[site]].owner != -1 && latticeChunks[idLattice[site]].owner != processorNumberLattice[site]) std::cout << "Fatal error in reallocating numbers!" << latticeChunks[idLattice[site]].owner << " " << processorNumberLattice[site] << std::endl;
				latticeChunks[idLattice[site]].owner = processorNumberLattice[site];
				latticeChunks[idLattice[site]].sharers = exchangeTable[site];
				if (idLattice[site] >= numberChunks) std::cout << "Fatal error in id too big" << std::endl;
			}
			
			MPI_Comm_rank(MPI_COMM_WORLD, &this_processor);
			
			localsize = loc_x*loc_y*loc_z*loc_t;
			completesize = localsize;
			sharedsize = 0;
			//Now we get the completesize (localsize + sharedsize)
			for (int i = 0; i < numberChunks; ++i) {
				for (unsigned int j = 0; j < latticeChunks[i].sharers.size(); ++j) {
					if (latticeChunks[i].sharers[j] == this_processor) {
						completesize += latticeChunks[i].size;
					}
				}
			}

			//Now we create the tags for the mpi communications
			std::map< std::pair<int, int> , int > bridgeCommunications;
			for (int i = 0; i < numberChunks; ++i) {
				for (unsigned int j = 0; j < latticeChunks[i].sharers.size(); ++j) {
					std::pair<int, int> pr(latticeChunks[i].owner, latticeChunks[j].sharers[j]);
					if (bridgeCommunications.count(pr)) {
						latticeChunks[i].tags.push_back(bridgeCommunications[pr]+1);
						++bridgeCommunications[pr];
					}
					else {
						latticeChunks[i].tags.push_back(0);
						bridgeCommunications[pr] = 0;
					}
				}
			}
			
			
			//We get the localIndex, for converting global index to local index, -1 if it is not present locally
			localIndex = new int[globalVolume];
#pragma omp parallel for
			for (int site = 0; site < globalVolume; ++site) localIndex[site] = -1;
			
			int offset = 0;
			int index = 0;
			//We set now also where a chunk starts in the local array
			//We start first with the local part that is shared
			for (int i = 0; i < numberChunks; ++i) {
				latticeChunks[i].offset = -1;
				if (latticeChunks[i].owner == this_processor && latticeChunks[i].sharers.size() != 0) {
					sharedsize += latticeChunks[i].size;
					latticeChunks[i].offset = offset;
					offset += latticeChunks[i].size;
					for (int site = 0; site < globalVolume; ++site) {
						if (idLattice[site] == latticeChunks[i].id) {
							localIndex[site] = index;
							++index;
						}
					}
				}
			}
			
			//Now we put the local part not shared
			for (int i = 0; i < numberChunks; ++i) {
				if (latticeChunks[i].owner == this_processor && latticeChunks[i].sharers.size() == 0) {
					latticeChunks[i].offset = offset;
					offset += latticeChunks[i].size;
					for (int site = 0; site < globalVolume; ++site) {
						if (idLattice[site] == latticeChunks[i].id) {
							localIndex[site] = index;
							++index;
						}
					}
				}
			}
			
			//Now offset should be equal to localsize
			if (offset != localsize && this_processor == 0) {
				std::cout << "Estimated localsize is different: " << offset << " " << localsize << "!" << std::endl;
				exit(7);
			}
			
			//Now we put the part that belongs to other processes but that is accessed also by this process
			for (int i = 0; i < numberChunks; ++i) {
				for (unsigned int j = 0; j < latticeChunks[i].sharers.size(); ++j) {
					if (latticeChunks[i].sharers[j] == this_processor) {
						latticeChunks[i].offset = offset;
						offset += latticeChunks[i].size;
						for (int site = 0; site < globalVolume; ++site) {
							if (idLattice[site] == latticeChunks[i].id) {
								localIndex[site] = index;
								++index;
							}
						}
					}
				}
			}
			
			//Now offset should be equal to completesize
			if (offset != completesize && this_processor == 0) {
				std::cout << "Estimated completesize is different: " << offset << " " << completesize << "!" << std::endl;
				exit(7);
			}
			
			//Now we construct the global sup table
			std::vector<Site> deltaPlus;
			deltaPlus.push_back(Site(1,0,0,0));
			deltaPlus.push_back(Site(0,1,0,0));
			deltaPlus.push_back(Site(0,0,1,0));
			deltaPlus.push_back(Site(0,0,0,1));
	
			vect4* global_sup_table = new vect4[globalVolume];
#pragma omp parallel for
			for (int x = 0; x < glob_x; ++x) {
				for (int y = 0; y < glob_y; ++y) {
					for (int z = 0; z < glob_z; ++z) {
						for (int t = 0; t < glob_t; ++t) {
							Site site(x,y,z,t);
							int index = getGlobalCoordinate(site);
							for (unsigned int i = 0; i < 4; ++i) {
								global_sup_table[index][i] = getGlobalCoordinate(site + deltaPlus[i]);
							}
						}
					}
				}
			}
	
			//Now we construct the global down table
			std::vector<Site> deltaMinus;
			deltaMinus.push_back(Site(-1,0,0,0));
			deltaMinus.push_back(Site(0,-1,0,0));
			deltaMinus.push_back(Site(0,0,-1,0));
			deltaMinus.push_back(Site(0,0,0,-1));
	
			vect4* global_down_table = new vect4[globalVolume];
#pragma omp parallel for
			for (int x = 0; x < glob_x; ++x) {
				for (int y = 0; y < glob_y; ++y) {
					for (int z = 0; z < glob_z; ++z) {
						for (int t = 0; t < glob_t; ++t) {
							Site site(x,y,z,t);
							int index = getGlobalCoordinate(site);
							for (unsigned int i = 0; i < 4; ++i) {
								global_down_table[index][i] = getGlobalCoordinate(site + deltaMinus[i]);
							}
						}
					}
				}
			}
	
			//Now we construct the local sup/down table
			sup_table = new vect4[completesize];
			down_table = new vect4[completesize];
			globalCoordinate = new Site[completesize];
			for (int i = 0; i < numberChunks; ++i) {
				if (latticeChunks[i].offset != -1) {
					for (int site = 0; site < globalVolume; ++site) {
						if (idLattice[site] == latticeChunks[i].id) {
							for (int k = 0; k < 4; ++k) {
								if (localIndex[site] == -1 || localIndex[site] > completesize) std::cout << "Fatal error in localindex!" << localIndex[site] << " " << completesize << " " << localsize << std::endl;
								if (global_sup_table[site][k] < 0 || global_sup_table[site][k] > globalVolume) std::cout << "Fatal error in global sup/down table!" << std::endl;
								if (global_down_table[site][k] < 0 || global_down_table[site][k] > globalVolume) std::cout << "Fatal error in global sup/down table!" << std::endl;
								sup_table[localIndex[site]][k] = localIndex[global_sup_table[site][k]];
								down_table[localIndex[site]][k] = localIndex[global_down_table[site][k]];
							}
							globalCoordinate[localIndex[site]] = globalMapSite[site];
						}
					}
				}
			}
	
			delete[] processorNumberLattice;
			delete[] globalMapSite;
			delete[] exchangeTable;
			delete[] idLattice;
			delete[] global_sup_table;
			delete[] global_down_table;
#endif
#ifndef ENABLE_MPI
			std::cout << "MPI is not activated at compile time!" << std::endl;
			exit(4);
#endif
		}
		
		template<typename T> static int* getMapIndex(const MpiLayout<T>& layout) {
			if (coversion_map.count(T::id)) return coversion_map[T::id];
			else {
				int* result = new int[localsize];
#pragma omp parallel for
				for (int site = 0; site < localsize; ++site) result[site] = -1;
#pragma omp parallel for
				for (int site = 0; site < globalVolume; ++site)	{
					if (localIndex[site] != -1 && localIndex[site] < localsize) {
						result[localIndex[site]] = layout.localIndex[site];
					}
				}
#pragma omp parallel for
				for (int site = 0; site < localsize; ++site) if (result[site] == -1) std::cout << "Fatal error in exchangeTable conversions!" << std::endl;
				coversion_map[T::id] = result;
				return result;
			}
		}

		static std::map<int,int*> coversion_map;

		static void destroy() {
			for (std::map<int,int*>::iterator it = coversion_map.begin(); it != coversion_map.end(); ++it) {
				delete[] it->second;
			}
			delete[] localIndex;
			delete[] sup_table;
			delete[] down_table;
			delete[] globalCoordinate;
		}
	
		static std::vector<LatticeChunk> latticeChunks;
		static int numberChunks;
		static int localsize;
		static int completesize;
		static int sharedsize;
		static int this_processor;
		//To every global site it is associated a local site, -1 if it is not present locally
		static int* localIndex;
		//The sup table
		static vect4* sup_table;
		//The down table
		static vect4* down_table;
		//The globalCoordinate of a local site
		static Site* globalCoordinate;
		
		static void printReport() {
			std::fstream outputfile;
			std::string namefile = "logfile_mpitable_";
			namefile += Update::toString(this_processor);
	
			outputfile.open(namefile.c_str(), std::fstream::out | std::fstream::app);
			
			outputfile << "Global volume: " << globalVolume << std::endl;
			outputfile << "Number of processors: " << numberProcessors << std::endl;
	
	
			outputfile << "Number of chunks: " << numberChunks << std::endl;
			int totalsite = 0;
			for (int i = 0; i < numberChunks; ++i) {
				outputfile << "Number of sites in the chunk " << i << ": " << latticeChunks[i].size << std::endl;
				outputfile << " the owner is: " << latticeChunks[i].owner << std::endl;
				outputfile << " the sharers are: ";
				for (unsigned int j = 0; j < latticeChunks[i].sharers.size(); ++j) {
					outputfile << latticeChunks[i].sharers[j] << " ";
				}
				outputfile << std::endl;
				totalsite += latticeChunks[i].size;
			}
			outputfile << "The number of totalsite in the chunks is " << totalsite << std::endl;
			outputfile.close();
		}
		
		static int getGlobalCoordinate(int xp, int yp, int zp, int tp) {
			int x = modulus(xp,glob_x);
			int y = modulus(yp,glob_y);
			int z = modulus(zp,glob_z);
			int t = modulus(tp,glob_t);
			return glob_t*(glob_z*(glob_y*x + y) + z) + t;
		}

		static int getGlobalCoordinate(const Site& site) {
			int x = modulus(site.x,glob_x);
			int y = modulus(site.y,glob_y);
			int z = modulus(site.z,glob_z);
			int t = modulus(site.t,glob_t);
			return glob_t*(glob_z*(glob_y*x + y) + z) + t;
		}

		static int rankTable(int x, int y, int z, int t) {
			int processor = x/loc_x;
			processor = pgrid_y*processor + y/loc_y;
			processor = pgrid_z*processor + z/loc_z;
			processor = pgrid_t*processor + t/loc_t;
			return processor;
		}

		static int globalIndexX(int site) {
			return globalCoordinate[site].x;
		}

		static int globalIndexY(int site) {
			return globalCoordinate[site].y;
		}

		static int globalIndexZ(int site) {
			return globalCoordinate[site].z;
		}

		static int globalIndexT(int site) {
			return globalCoordinate[site].t;
		}

		static int globalIndex(int site, int mu) {
			return (reinterpret_cast<int*>(&(globalCoordinate[site])))[mu];
		}

		static void save(const std::string& basename) {
			FILE* fout(NULL);
			//First we save the sup table
			std::string output_file = basename + ".suptable_" + Update::toString(this_processor) + ".txt";
			
			fout = fopen(output_file.c_str(), "w");
			
			if (!fout) {
				std::cout << "File not writeble!" << std::endl;
				return;
			}

			XDR xout;
			xdrstdio_create(&xout, fout, XDR_ENCODE);

			for (int site = 0; site < completesize; ++site) {
				xdr_int(&xout, &sup_table[site][0]);
				xdr_int(&xout, &sup_table[site][1]);
				xdr_int(&xout, &sup_table[site][2]);
				xdr_int(&xout, &sup_table[site][3]);
			}

			xdr_destroy(&xout);
			fclose(fout);
			
			//then we store the down table
			output_file = basename + ".downtable_" + Update::toString(this_processor) + ".txt";

			fout = fopen(output_file.c_str(), "w");
			
			if (!fout) {
				std::cout << "File not writeble!" << std::endl;
				return;
			}

			xdrstdio_create(&xout, fout, XDR_ENCODE);

			for (int site = 0; site < completesize; ++site) {
				xdr_int(&xout, &down_table[site][0]);
				xdr_int(&xout, &down_table[site][1]);
				xdr_int(&xout, &down_table[site][2]);
				xdr_int(&xout, &down_table[site][3]);
			}

			xdr_destroy(&xout);
			fclose(fout);
			
			//Then we store the global sites
			output_file = basename + ".global_coordinate_" + Update::toString(this_processor) + ".txt";

			fout = fopen(output_file.c_str(), "w");
			
			if (!fout) {
				std::cout << "File not writeble!" << std::endl;
				return;
			}

			xdrstdio_create(&xout, fout, XDR_ENCODE);

			for (int site = 0; site < completesize; ++site) {
				xdr_int(&xout, &globalCoordinate[site].x);
				xdr_int(&xout, &globalCoordinate[site].y);
				xdr_int(&xout, &globalCoordinate[site].z);
				xdr_int(&xout, &globalCoordinate[site].t);
			}

			xdr_destroy(&xout);
			fclose(fout);

			//Then we store the global index
			output_file = basename + ".local_index_" + Update::toString(this_processor) + ".txt";

			fout = fopen(output_file.c_str(), "w");
			
			if (!fout) {
				std::cout << "File not writeble!" << std::endl;
				return;
			}

			xdrstdio_create(&xout, fout, XDR_ENCODE);

			for (int site = 0; site < globalVolume; ++site) {
				xdr_int(&xout, &localIndex[site]);
			}

			xdr_destroy(&xout);
			fclose(fout);

			//Now we store the datas of the chunks
			output_file = basename + ".chunk_" + Update::toString(this_processor) + ".txt";
			std::fstream chunksfile;
			chunksfile.open(output_file.c_str(), std::fstream::out);
			chunksfile << latticeChunks.size() << std::endl;
			for (unsigned int i = 0; i < latticeChunks.size(); ++i) {
				chunksfile << latticeChunks[i].id << " " << latticeChunks[i].owner << " " << latticeChunks[i].size << " " << latticeChunks[i].offset << " " << latticeChunks[i].sharers.size() << std::endl;
				for (unsigned int j = 0; j < latticeChunks[i].sharers.size(); ++j) {
					chunksfile << latticeChunks[i].sharers[j] << " " << latticeChunks[i].tags[j] << std::endl;
				}
			}
			chunksfile.close();
			
			if (this_processor == 0) {
				output_file = basename + ".descriptor.txt";
				std::fstream descriptor;
				descriptor.open(output_file.c_str(), std::fstream::out);
				descriptor << glob_x << " " << glob_y << " " << glob_z << " " << glob_t << std::endl;
				descriptor << pgrid_x << " " << pgrid_y << " " << pgrid_z << " " << pgrid_t << std::endl;
				descriptor << numberProcessors << std::endl;
				descriptor << localsize << " " << sharedsize << " " << completesize << std::endl;
				descriptor.close();
			}
		}

		static void load(const std::string& basename) {
			numberProcessors = pgrid_x*pgrid_y*pgrid_z*pgrid_t;
			globalVolume = glob_x*glob_y*glob_z*glob_t;
			glob_spatial_volume = glob_x*glob_y*glob_z;

			glob[0] = glob_x;
                        glob[1] = glob_y;
                        glob[2] = glob_z;
                        glob[3] = glob_t;
#ifdef ENABLE_MPI
			int sizeCommunicator = 0;
			MPI_Comm_size(MPI_COMM_WORLD, &sizeCommunicator);
			MPI_Comm_rank(MPI_COMM_WORLD, &this_processor);
			
			if (sizeCommunicator != numberProcessors) {
				std::cout << "Number of processors " << numberProcessors << " is different from the size of the communicator " << sizeCommunicator << std::endl;
				exit(1);
			}
			
			loc_x = glob_x / pgrid_x;
			loc_y = glob_y / pgrid_y;
			loc_z = glob_z / pgrid_z;
			loc_t = glob_t / pgrid_t;
			
			//Check that everything is ok
			if ((glob_x % pgrid_x) != 0) {
				std::cout << "The global size " << glob_x << " is not divisible in the grid size " << pgrid_x << std::endl;
				exit(3);
			}
			if ((glob_y % pgrid_y) != 0) {
				std::cout << "The global size " << glob_y << " is not divisible in the grid size " << pgrid_y << std::endl;
				exit(3);
			}
			if ((glob_z % pgrid_z) != 0) {
				std::cout << "The global size " << glob_z << " is not divisible in the grid size " << pgrid_z << std::endl;
				exit(3);
			}
			if ((glob_t % pgrid_t) != 0) {
				std::cout << "The global size " << glob_t << " is not divisible in the grid size " << pgrid_t << std::endl;
				exit(3);
			}

			std::string input_name = basename + ".descriptor.txt";
			int read_glob_x, read_glob_y, read_glob_z, read_glob_t;
			int read_pgrid_x, read_pgrid_y, read_pgrid_z, read_pgrid_t;
			int read_numberProcessors;

			std::fstream descriptor;
			descriptor.open(input_name.c_str(), std::fstream::in);

			descriptor >> read_glob_x >> read_glob_y >> read_glob_z >> read_glob_t;
			descriptor >> read_pgrid_x >> read_pgrid_y >> read_pgrid_z >> read_pgrid_t;
			descriptor >> read_numberProcessors;
			descriptor >> localsize >> sharedsize >> completesize;

			if (read_glob_x != glob_x || read_glob_y != glob_y || read_glob_z != glob_z || read_glob_t != glob_t) {
				std::cout << "Different lattice size in reading configuration!" << std::endl;
				std::cout << "Configured: " << glob_x << " " << glob_y << " " << glob_z << " " << glob_t << std::endl;
				std::cout << "Readed: " << read_glob_x << " " << read_glob_y << " " << read_glob_z << " " << read_glob_t << std::endl;
				exit(33);
			}

			if (read_pgrid_x != pgrid_x || read_pgrid_y != pgrid_y || read_pgrid_z != pgrid_z || read_pgrid_t != pgrid_t) {
				std::cout << "Different grid size in reading configuration!" << std::endl;
				std::cout << "Configured: " << pgrid_x << " " << pgrid_y << " " << pgrid_z << " " << pgrid_t << std::endl;
				std::cout << "Readed: " << read_pgrid_x << " " << read_pgrid_y << " " << read_pgrid_z << " " << read_pgrid_t << std::endl;
				exit(33);
			}

			descriptor.close();

			std::string input_file = basename + ".chunk_" + Update::toString(this_processor) + ".txt";
			std::fstream chunksfile;
			chunksfile.open(input_file.c_str(), std::fstream::in);
			chunksfile >> numberChunks;
			latticeChunks.resize(numberChunks);
			for (unsigned int i = 0; i < latticeChunks.size(); ++i) {
				int sharerSize; 
				chunksfile >> latticeChunks[i].id >> latticeChunks[i].owner >> latticeChunks[i].size >> latticeChunks[i].offset >> sharerSize;
				latticeChunks[i].sharers.resize(sharerSize);
				latticeChunks[i].tags.resize(sharerSize);
				for (unsigned int j = 0; j < latticeChunks[i].sharers.size(); ++j) {
					chunksfile >> latticeChunks[i].sharers[j] >> latticeChunks[i].tags[j];
				}
			}
			chunksfile.close();

			//Now we read the globalCoordinate
			input_file = basename + ".local_index_" + Update::toString(this_processor) + ".txt";
			
			FILE* fin(NULL);
			fin = fopen(input_file.c_str(), "r");

			if (!fin) {
				std::cout << "File not readble!" << std::endl;
				return;
			}

			XDR xin;
			xdrstdio_create(&xin, fin, XDR_DECODE);

			localIndex = new int[globalVolume];
			for (int site = 0; site < globalVolume; ++site) {
				xdr_int(&xin, &localIndex[site]);
			}

			xdr_destroy(&xin);
			fclose(fin);

			//Now we read the global coordinate
			input_file = basename + ".global_coordinate_" + Update::toString(this_processor) + ".txt";
			fin = fopen(input_file.c_str(), "r");

			if (!fin) {
				std::cout << "File not readble!" << std::endl;
				return;
			}

			xdrstdio_create(&xin, fin, XDR_DECODE);

			globalCoordinate = new Site[completesize];
			for (int site = 0; site < completesize; ++site) {
				xdr_int(&xin, &globalCoordinate[site].x);
				xdr_int(&xin, &globalCoordinate[site].y);
				xdr_int(&xin, &globalCoordinate[site].z);
				xdr_int(&xin, &globalCoordinate[site].t);
			}

			xdr_destroy(&xin);
			fclose(fin);
			
			//Now we read the downtable
			input_file = basename + ".downtable_" + Update::toString(this_processor) + ".txt";
			fin = fopen(input_file.c_str(), "r");

			if (!fin) {
				std::cout << "File not readble!" << std::endl;
				return;
			}

			xdrstdio_create(&xin, fin, XDR_DECODE);

			down_table = new vect4[completesize];
			for (int site = 0; site < completesize; ++site) {
				xdr_int(&xin, &down_table[site][0]);
				xdr_int(&xin, &down_table[site][1]);
				xdr_int(&xin, &down_table[site][2]);
				xdr_int(&xin, &down_table[site][3]);
			}

			xdr_destroy(&xin);
			fclose(fin);

			//Now we read the suptable
			input_file = basename + ".suptable_" + Update::toString(this_processor) + ".txt";
			fin = fopen(input_file.c_str(), "r");

			if (!fin) {
				std::cout << "File not readble!" << std::endl;
				return;
			}

			xdrstdio_create(&xin, fin, XDR_DECODE);

			sup_table = new vect4[completesize];
			for (int site = 0; site < completesize; ++site) {
				xdr_int(&xin, &sup_table[site][0]);
				xdr_int(&xin, &sup_table[site][1]);
				xdr_int(&xin, &sup_table[site][2]);
				xdr_int(&xin, &sup_table[site][3]);
			}

			xdr_destroy(&xin);
			fclose(fin);
#endif
#ifndef ENABLE_MPI
			std::cout << "MPI is not activated at compile time!" << std::endl;
			exit(4);
#endif
		}
		
	private:
		static int modulus(int value, int mod) {
			int ris = value;
			if (ris >= mod) return modulus(ris - mod, mod);
			else if (ris < 0) return modulus(ris + mod, mod);
			else return ris;
		}
		
};


template<typename Stencil> int MpiLayout<Stencil>::pgrid_t = 0;
template<typename Stencil> int MpiLayout<Stencil>::pgrid_x = 0;
template<typename Stencil> int MpiLayout<Stencil>::pgrid_y = 0;
template<typename Stencil> int MpiLayout<Stencil>::pgrid_z = 0;

template<typename Stencil> int MpiLayout<Stencil>::glob_t = 0;
template<typename Stencil> int MpiLayout<Stencil>::glob_x = 0;
template<typename Stencil> int MpiLayout<Stencil>::glob_y = 0;
template<typename Stencil> int MpiLayout<Stencil>::glob_z = 0;
template<typename Stencil> int MpiLayout<Stencil>::glob[4];
		
template<typename Stencil> int MpiLayout<Stencil>::loc_t = 0;
template<typename Stencil> int MpiLayout<Stencil>::loc_x = 0;
template<typename Stencil> int MpiLayout<Stencil>::loc_y = 0;
template<typename Stencil> int MpiLayout<Stencil>::loc_z = 0;
		
template<typename Stencil> int MpiLayout<Stencil>::numberProcessors = 0;
template<typename Stencil> int MpiLayout<Stencil>::globalVolume = 0;
template<typename Stencil> int MpiLayout<Stencil>::glob_spatial_volume = 0;

template<typename Stencil> std::vector<LatticeChunk> MpiLayout<Stencil>::latticeChunks;
template<typename Stencil> int MpiLayout<Stencil>::numberChunks = 0;
template<typename Stencil> int MpiLayout<Stencil>::localsize = 0;
template<typename Stencil> int MpiLayout<Stencil>::completesize = 0;
template<typename Stencil> int MpiLayout<Stencil>::sharedsize = 0;
template<typename Stencil> int MpiLayout<Stencil>::this_processor = 0;

template<typename Stencil> int* MpiLayout<Stencil>::localIndex = 0;
template<typename Stencil> vect4* MpiLayout<Stencil>::sup_table = 0;
template<typename Stencil> vect4* MpiLayout<Stencil>::down_table = 0;
template<typename Stencil> Site* MpiLayout<Stencil>::globalCoordinate = 0;

template<typename Stencil> std::map<int,int*> MpiLayout<Stencil>::coversion_map;

}

#endif
