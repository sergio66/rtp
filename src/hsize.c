
/* return the size in bytes of some common HDF types
 * Created circa 18 March 2002 by Howard Motteler
 * Update: 20 Oct 2011, Scott Hannon - added DFNT_UINT8
 */

#include <stdio.h>
#include "hdf.h"

int hsize (int htype) {

  switch (htype) {
  case DFNT_INT32:
    return 4;
  case DFNT_FLOAT32:
    return 4;
  case DFNT_FLOAT64:
    return 8;
  case DFNT_UINT8:
    return 1;
  case DFNT_CHAR8:
    return 1;
  case DFNT_UCHAR8:
    return 1;
  default:
    fprintf(stderr, 
	    "hsize(): WARNING - unexpected data type code %d\n",
	    htype);
    return 0;
  }
}

