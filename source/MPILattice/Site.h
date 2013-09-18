#ifndef SITE_H
#define SITE_H

namespace Lattice {

struct Site {
	int x;
	int y;
	int z;
	int t;
	Site(int _x, int _y, int _z, int _t) : x(_x), y(_y), z(_z), t(_t) { }
	Site operator+(const Site& snd) const {
		return Site(x+snd.x, y+snd.y, z+snd.z, t+snd.t);
	}
	Site operator-(const Site& snd) const {
		return Site(x-snd.x, y-snd.y, z-snd.z, t-snd.t);
	}
	Site operator-() const {
		return Site(-x, -y, -z, -t);
	}
	bool operator==(const Site& toCompare) const {
		return (x == toCompare.x) && (y == toCompare.y) && (z == toCompare.z) && (t == toCompare.t);
	}
	Site() : x(0), y(0), z(0), t(0) { }
};

}

#endif
