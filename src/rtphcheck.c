
/* rtphcheck -- sanity check for RTP header values
 * 
 * returns -1 on error, 0 on success
 * Version 2.01
 */

#include <stdlib.h>
#include "hdf.h"

#include "rtp.h"
#include "rtpfnx.h"

int rtphcheck(struct rtp_head *head) {

  int i;

  /* check ptype and pfields 
   */
  if (head->ptype != LEVPRO 
      && head->ptype != LAYPRO
      && head->ptype != AIRSLAY) {
    fprintf(stderr, "header field ptype = %d, bad value\n", 
	    head->ptype);
    return (-1);
  }

  if (head->pfields < 0 || PFIELDSMAX < head->pfields) {
    fprintf(stderr, "header field pfields = %d, out of range\n", 
	    head->pfields);
    return (-1);
  }


  /* check header size fields 
   */
  if (head->ngas  < 0 || MAXGAS  < head->ngas) {
    fprintf(stderr, "header field ngas = %d, out of range\n",
	    head->ngas);
    return (-1);
  }

  if (head->nchan < 0 || MAXCHAN < head->nchan ) {
    fprintf(stderr, "header field nchan = %d, out of range\n",
	    head->nchan);
    return (-1);
  }


  /* check the constituent list 
   */
  for (i=0; i < head->ngas; i++) {
    if (head->glist[i] < 0 || MAXGASID < head->glist[i]) {
      fprintf(stderr, "bad header gas ID %d\n ", head->glist[i]);
      return (-1);
    }
  }


  /* check for pfields and nchan consistency 
   */
  if ((head->pfields & IRCALCBIT || head->pfields & IROBSVBIT) && 
      head->nchan == 0) {
    fprintf(stderr, "header field pfields implies rad's but nchan = 0\n");
    fprintf(stderr, "pfields=%d, nchan=%d\n", head->pfields, head->nchan);
    return (-1);
  }


  /* check header max size fields -- these don't go into the vdata,
   * but they are used to set the size of the vdata profile fields
   */

  if (head->memis < 0 || MAXEMIS < head->memis) {
    fprintf(stderr, "header field memis = %d, out of range\n",
            head->memis);
    return (-1);
  }

  if (head->mlevs < 0 || MAXLEV  < head->mlevs) {
    fprintf(stderr, "header field mlevs = %d, out of range\n",
            head->mlevs);
    return (-1);
  }

  /* all tests passed OK, return success 
   */
  return (0);
}

