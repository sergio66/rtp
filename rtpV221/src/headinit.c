
/* initialize an RTP header structure
 */

#include "rtp.h"

void headinit(struct rtp_head *head) {

  int j;

  /* profile data */
  head->ptype           = LEVPRO;
  head->pfields         = PROFBIT;

  head->pmin            = BAD;
  head->pmax            = BAD;
  head->ngas            = 0;
  head->glist[0]        = BAD;
  head->gunit[0]        = BAD;

  /* radiance data */
  head->pltfid          = BAD;
  head->instid          = BAD;
  head->nchan           = 0;
  head->ichan[0]        = BAD;
  head->vchan[0]        = BAD;
  head->vcmin           = BAD;
  head->vcmax           = BAD;

  /* maxes for profile fields */
  head->memis           = 0;
  head->mlevs           = 0;

  /* user-defined fields */
  for (j=0; j < MAXIUDEF; j++)  head->iudef[j]  = BAD;
  head->itype           = BAD;

}

