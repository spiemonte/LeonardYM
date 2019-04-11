#ifndef LOCALLAYOUT_H
#define LOCALLAYOUT_H
#include "Site.h"
#include <vector>
#include <iostream>

namespace Lattice {

typedef int vect4[4];

class LocalLayout {
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
			numberProcessors = 1;
			this_processor = 0;
			globalVolume = glob_x*glob_y*glob_z*glob_t;
			glob_spatial_volume = glob_x*glob_y*glob_z;
		
			glob[0] = glob_x;
			glob[1] = glob_y;
			glob[2] = glob_z;
			glob[3] = glob_t;
			
			pgrid_x = 1;
			pgrid_y = 1;
			pgrid_z = 1;
			pgrid_t = 1;
			
			loc_x = glob_x;
			loc_y = glob_y;
			loc_z = glob_z;
			loc_t = glob_t;
			
			localsize = loc_x*loc_y*loc_z*loc_t;
			completesize = localsize;
			sharedsize = 0;
			
			//Now we construct the local sup/down table
			sup_table = new vect4[completesize];
			down_table = new vect4[completesize];
			globalCoordinate = new Site[completesize];
			
			//Now we construct the local sup table
			std::vector<Site> deltaPlus;
			deltaPlus.push_back(Site(1,0,0,0));
			deltaPlus.push_back(Site(0,1,0,0));
			deltaPlus.push_back(Site(0,0,1,0));
			deltaPlus.push_back(Site(0,0,0,1));
	
			for (int x = 0; x < glob_x; ++x) {
				for (int y = 0; y < glob_y; ++y) {
					for (int z = 0; z < glob_z; ++z) {
						for (int t = 0; t < glob_t; ++t) {
							Site site(x,y,z,t);
							int index = getGlobalCoordinate(site);
							for (unsigned int i = 0; i < 4; ++i) {
								sup_table[index][i] = getGlobalCoordinate(site + deltaPlus[i]);
							}
							globalCoordinate[index] = site;
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
	
			for (int x = 0; x < glob_x; ++x) {
				for (int y = 0; y < glob_y; ++y) {
					for (int z = 0; z < glob_z; ++z) {
						for (int t = 0; t < glob_t; ++t) {
							Site site(x,y,z,t);
							int index = getGlobalCoordinate(site);
							for (unsigned int i = 0; i < 4; ++i) {
								down_table[index][i] = getGlobalCoordinate(site + deltaMinus[i]);
							}
						}
					}
				}
			}

			localIndex = new int[completesize];
			for (int site = 0; site < completesize; ++site) localIndex[site] = site;
			
		}
	
		static int localsize;
		static int completesize;
		static int sharedsize;
		static int this_processor;
		//The sup table
		static vect4* sup_table;
		//The down table
		static vect4* down_table;
		//The globalCoordinate of a local site
		static Site* globalCoordinate;
		//The local index of a global one, not needed, only for compatibility with mpi
		static int* localIndex;

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
		
		static void printReport() {
			std::cout << "Global volume: " << globalVolume << std::endl;
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
	private:
		static int modulus(int value, int mod) {
			int ris = value;
			if (ris >= mod) return ris - mod;
			else if (ris < 0) return ris + mod;
			else return ris;
		}
};

}

#endif
