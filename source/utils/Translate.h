#ifndef TRANSLATE_H
#define TRANSLATE_H

namespace Update {

template<typename Lattice> void translate(const Lattice& toTranslate, Lattice& translated, int dx, int dy, int dz, int dt) {
	if (dx == 0 && dy == 0 && dz == 0 && dt == 0) {
		translated = toTranslate;
		return;
	}
	Lattice swap = toTranslate;
	for (int i = 0; i < dx; ++i) {
#pragma omp parallel for
		for (int site = 0; site < toTranslate.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				translated[site][mu] = swap[Lattice::sup(site,0)][mu];
			}
		}
		translated.updateHalo();
		swap = translated;
	}

	for (int i = 0; i < dy; ++i) {
#pragma omp parallel for
		for (int site = 0; site < toTranslate.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				translated[site][mu] = swap[Lattice::sup(site,1)][mu];
			}
		}
		translated.updateHalo();
		swap = translated;
	}

	for (int i = 0; i < dz; ++i) {
#pragma omp parallel for
		for (int site = 0; site < toTranslate.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				translated[site][mu] = swap[Lattice::sup(site,2)][mu];
			}
		}
		translated.updateHalo();
		swap = translated;
	}

	for (int i = 0; i < dt; ++i) {
#pragma omp parallel for
		for (int site = 0; site < toTranslate.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				translated[site][mu] = swap[Lattice::sup(site,3)][mu];
			}
		}
		translated.updateHalo();
		swap = translated;
	}
}

}

#endif

