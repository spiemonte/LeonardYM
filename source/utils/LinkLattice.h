#ifndef LINKLATTICE_H_
#define LINKLATTICE_H_
#include <iostream>
#include <cstdlib>

namespace Update {

template<typename T> class PrivateLattice {
public:
	PrivateLattice() : sup_lookuptable(0), lengthX(4), lengthY(4), lengthZ(4), lengthT(4) {
		if (sup_lookuptable == 0) {
			initialize();
		}
		lattice = new T[lengthX*lengthY*lengthZ*lengthT];
	}

	PrivateLattice(const PrivateLattice& toCopy) : lengthX(toCopy.lengthX), lengthY(toCopy.lengthY), lengthZ(toCopy.lengthZ), lengthT(toCopy.lengthT), localsize(toCopy.localsize) {
		sup_lookuptable = new unsigned int[lengthX*lengthY*lengthZ*lengthT][4];
		sdn_lookuptable = new unsigned int[lengthX*lengthY*lengthZ*lengthT][4];
		for (int site = 0; site < lengthX*lengthY*lengthZ*lengthT; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				sup_lookuptable[site][mu] = toCopy.sup_lookuptable[site][mu];
				sdn_lookuptable[site][mu] = toCopy.sdn_lookuptable[site][mu];
			}
		}
		lattice = new T[lengthX*lengthY*lengthZ*lengthT];
		for (unsigned int site = 0; site < localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				lattice[site][mu] = toCopy.lattice[site][mu];
			}
		}
	}

	~PrivateLattice() {
		delete[] lattice;
		delete[] sup_lookuptable;
		delete[] sdn_lookuptable;
	}

	PrivateLattice& operator=(const PrivateLattice& toCopy) {
		if (lengthX != toCopy.lengthX || lengthY != toCopy.lengthY || lengthZ != toCopy.lengthZ || lengthT != toCopy.lengthT) {
			lengthX = toCopy.lengthX;
			lengthY = toCopy.lengthY;
			lengthZ = toCopy.lengthZ;
			lengthT = toCopy.lengthT;
			localsize = toCopy.localsize;

			delete[] lattice;
			delete[] sup_lookuptable;
			delete[] sdn_lookuptable;

			sup_lookuptable = new unsigned int[lengthX*lengthY*lengthZ*lengthT][4];
			sdn_lookuptable = new unsigned int[lengthX*lengthY*lengthZ*lengthT][4];
			for (int site = 0; site < lengthX*lengthY*lengthZ*lengthT; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					sup_lookuptable[site][mu] = toCopy.sup_lookuptable[site][mu];
					sdn_lookuptable[site][mu] = toCopy.sdn_lookuptable[site][mu];
				}
			}
			lattice = new T[lengthX*lengthY*lengthZ*lengthT];
		}
		for (unsigned int site = 0; site < localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				lattice[site][mu] = toCopy.lattice[site][mu];
			}
		}
		return *this;
	}

	unsigned int sup(unsigned int site, unsigned int mu) {
		return sup_lookuptable[site][mu];
	}

	unsigned int sdn(unsigned int site, unsigned int mu) {
		return sdn_lookuptable[site][mu];
	}

	T* operator[](unsigned int index) {
		return lattice[index];
	}

	const T* operator[](unsigned int index) const {
		return lattice[index];
	}

	T* operator()(int x, int y, int z, int t) {
		return lattice[(((t*lengthX)+x)*lengthY + y)*lengthZ + z];
	}

	const T* operator()(int x, int y, int z, int t) const {
		return lattice[(((t*lengthX)+x)*lengthY + y)*lengthZ + z];
	}

	unsigned int spatialVolume() const {
		return lengthX*lengthY*lengthZ;
	}

	unsigned int globalVolume() const {
		return lengthX*lengthY*lengthZ*lengthT;
	}

	void setLatticeSize(int _lengthX, int _lengthY, int _lengthZ, int _lengthT) {
		lengthX = _lengthX;
		lengthY = _lengthY;
		lengthZ = _lengthZ;
		lengthT = _lengthT;
		localsize = lengthX*lengthY*lengthZ*lengthT;
		this->initialize();
	}

	/*class Layout {
	public:
		static unsigned int pgrid_t;
		static unsigned int pgrid_x;
		static unsigned int pgrid_y;
		static unsigned int pgrid_z;
		static unsigned int glob_t;
		static unsigned int glob_x;
		static unsigned int glob_y;
		static unsigned int glob_z;
		static unsigned int loc_t;
		static unsigned int loc_x;
		static unsigned int loc_y;
		static unsigned int loc_z;

		static int rankT() { return 0; }
		static int index(int x, int y, int z, int t) {
			return (((t*loc_x)+x)*loc_y + y)*loc_z + z;
		}

		static unsigned int glob_spatial_volume;
	};*/
private:
	T *lattice;

	unsigned int (* sup_lookuptable)[4];
	unsigned int (* sdn_lookuptable)[4];

	struct Coordinate {
		int x;
		int y;
		int z;
		int t;

		Coordinate(int _x, int _y, int _z, int _t) : x(_x), y(_y), z(_z), t(_t) { }

		void normalize() {
			x = modulus(x, lengthX);
			y = modulus(y, lengthY);
			z = modulus(z, lengthZ);
			t = modulus(t, lengthT);
		}

		unsigned int toNumber() const {
			return (((t*lengthX)+x)*lengthY + y)*lengthZ + z;
		}

		int modulus(int value, int mod) const {
			int ris = value;
			if (ris >= mod) return ris - mod;
			else if (ris < 0) return ris + mod;
			else return ris;
		}
	};

	void initialize() {
		std::cout << "Lattice size: " << lengthX << "x" << lengthY << "x" <<  lengthZ << "x" <<  lengthT << std::endl;
		sup_lookuptable = new unsigned int[lengthX*lengthY*lengthZ*lengthT][4];
		sdn_lookuptable = new unsigned int[lengthX*lengthY*lengthZ*lengthT][4];
		for (int x = 0; x < lengthX; ++x) {
			for (int y = 0; y < lengthY; ++y) {
				for (int z = 0; z < lengthZ; ++z) {
					for (int t = 0; t < lengthT; ++t) {
						Coordinate coord(x,y,z,t);

						Coordinate coord_supx(x+1,y,z,t);
						coord_supx.normalize();
						sup_lookuptable[coord.toNumber()][0] = coord_supx.toNumber();
						Coordinate coord_supy(x,y+1,z,t);
						coord_supy.normalize();
						sup_lookuptable[coord.toNumber()][1] = coord_supy.toNumber();
						Coordinate coord_supz(x,y,z+1,t);
						coord_supz.normalize();
						sup_lookuptable[coord.toNumber()][2] = coord_supz.toNumber();
						Coordinate coord_supt(x,y,z,t+1);
						coord_supt.normalize();
						sup_lookuptable[coord.toNumber()][3] = coord_supt.toNumber();

						Coordinate coord_sdnx(x-1,y,z,t);
						coord_sdnx.normalize();
						sdn_lookuptable[coord.toNumber()][0] = coord_sdnx.toNumber();
						Coordinate coord_sdny(x,y-1,z,t);
						coord_sdny.normalize();
						sdn_lookuptable[coord.toNumber()][1] = coord_sdny.toNumber();
						Coordinate coord_sdnz(x,y,z-1,t);
						coord_sdnz.normalize();
						sdn_lookuptable[coord.toNumber()][2] = coord_sdnz.toNumber();
						Coordinate coord_sdnt(x,y,z,t-1);
						coord_sdnt.normalize();
						sdn_lookuptable[coord.toNumber()][3] = coord_sdnt.toNumber();
					}
				}
			}
		}

		for (unsigned int site = 0; site < lengthX*lengthY*lengthZ*lengthT; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				if (sup_lookuptable[sdn_lookuptable[site][mu]][mu] - site != 0) {
					std::cout << "Lookup table test failed!" << std::endl;
					exit(1);
				}
			}
		}
	}

	unsigned int lengthX;
	unsigned int lengthY;
	unsigned int lengthZ;
	unsigned int lengthT;
	unsigned int localsize;

};

} /* namespace Update */
#endif /* LINKLATTICE_H_ */
