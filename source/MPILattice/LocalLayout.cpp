#include "LocalLayout.h"

namespace Lattice {

int LocalLayout::pgrid_t = 0;
int LocalLayout::pgrid_x = 0;
int LocalLayout::pgrid_y = 0;
int LocalLayout::pgrid_z = 0;

int LocalLayout::glob_t = 0;
int LocalLayout::glob_x = 0;
int LocalLayout::glob_y = 0;
int LocalLayout::glob_z = 0;
		
int LocalLayout::loc_t = 0;
int LocalLayout::loc_x = 0;
int LocalLayout::loc_y = 0;
int LocalLayout::loc_z = 0;
		
int LocalLayout::numberProcessors = 0;
int LocalLayout::globalVolume = 0;
int LocalLayout::glob_spatial_volume = 0;

int LocalLayout::localsize = 0;
int LocalLayout::completesize = 0;
int LocalLayout::sharedsize = 0;
int LocalLayout::this_processor = 0;

vect4* LocalLayout::sup_table = 0;
vect4* LocalLayout::down_table = 0;
Site* LocalLayout::globalCoordinate = 0;

int* LocalLayout::localIndex = 0;

}