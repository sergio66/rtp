
/* rtpsize -- structure sizes and basic sanity checks
 *
 * USAGE
 *   rtpsize
 *
 * DESCRIPTION
 *   some sanity checks on the rtp structures and field lists
 *
 * BUGS
 *   rtp apps will often work even with size mismatches
 *
 *   the xspan calc should probably be commented out, it's just a
 *   sanity check for using the field offsets, and would need to be
 *   updated after any field mods
 *
 *   if (say after some MAX param mods) rtime does fall on an odd
 *   word boundary, this can be fixed by either adding or commenting
 *   out the fill field
 */

#include <stdio.h>
#include <stdlib.h>
#include "hdf.h"
#include "rtp.h"
#include "rtpfnx.h"

struct rtp_prof testprof;
struct rtp_head testhead;

main () {

  int i, j, k;
  int hfsize, pfsize, hspan, pspan, fspan, rspan, xspan;

  /* sanity check that flist and structure sizes match, note that
     we need the first and last fields here to calculate the actual
     size of the structure */

  /* header flist size */
  for (i = 0, j=0; i < NHFIELD; i++)
    j += hfield[i].order * hsize(hfield[i].htype);
  hfsize = j;

  /* header structure spanning size */
  hspan = (char *)&(testhead.itype) - (char *)&(testhead.ptype) + 4;

  if (hfsize != hspan) {
    fprintf(stdout,
	    "ERROR: header flist size does not match structure span!\n");
    fprintf(stdout,
	    "ERROR: header flist size = %d, header structure span = %d\n",
	    hfsize, hspan);
  }

  /* profile flist size */
  for (i = 0, j=0; i < NPFIELD; i++)
    j += pfield[i].order * hsize(pfield[i].htype);
  pfsize = j;

  /* profile structure spanning size */
  pspan = (char *)&(testprof.itype) - (char *)&(testprof.plat) + 4;

  /* profile structure span through rlat and rlon */
  fspan = (char *)&(testprof.rlon) - (char *)&(testprof.plat) + 4;

  /* profile structure span to start of rtime */
  rspan = (char *)&(testprof.rtime) - (char *)&(testprof.plat);

  if (pfsize != pspan) {
    fprintf(stdout,
	    "WARNING: profile flist size does not match structure span\n");
    fprintf(stdout,
	    "WARNING: profile flist size = %d, profile structure span = %d\n",
	    pfsize, pspan);
  }

  fprintf(stdout, "the RTP header structure size is %d bytes\n", 
	  i = sizeof(struct rtp_head));

  fprintf(stdout, "the RTP header structure span is %d bytes\n", 
	  hspan);

  fprintf(stdout, "the RTP profile structure size is %d bytes\n", 
	  i = sizeof(struct rtp_prof));

  fprintf(stdout, "the RTP profile structure span is %d bytes\n", 
	  pspan);

  fprintf(stdout, "the RTP profile struct span through rlon is %d words\n", 
	  fspan/4);

  fprintf(stdout, "the RTP profile span to start of rtime is %d words\n", 
	  rspan/4);

  /* the following is a "by hand" count of words in the rtp.h
     profile struct through the rlat and rlon fields; it should
     match fspan.  This is really just a one-time sanity check that
     we trust the offset calcs
  */
  xspan = 4 + 7 + 3*MAXEMIS + 1 + 3*MAXLEV + MAXGAS*MAXLEV + 2*MAXGAS
    + 4 + 3*MAXLEV + 2 + 2*MAXEMIS + 7 + 2*MAXEMIS + 16 + 2;

  if (fspan/4 != xspan) {
    fprintf(stdout,
	    "WARNING: xspan != fspan/4\n");
    fprintf(stdout,
	    "WARNING: xspan = %d, fspan/4 = %d\n",
	    xspan, fspan/4);
  }

  if ((rspan/4) % 2 == 1)
    fprintf(stdout, "WARNING: rtime falls on an odd word boundary\n");

  /*
  fprintf(stdout, "size of header flist is %d bytes, %d fields\n",
          i = sizeof((struct FLIST) hfield),
          i = sizeof((struct FLIST) hfield)/(sizeof(char *)+sizeof(int)*2));

  fprintf(stdout, "size of profile flist is %d bytes, %d fields\n",
          i = sizeof((struct FLIST) pfield),
          i = sizeof((struct FLIST) pfield)/(sizeof(char*)+sizeof(int)*2));

  fprintf(stdout, "static space for attributes is %d bytes\n", 
	  i = sizeof(struct rtpfatt) * MAXNATTR);
  */
}

