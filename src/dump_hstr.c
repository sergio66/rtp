
/* dump the RTP header structure Version 2.01 */

#include <stdio.h>
#include "rtp.h"

void dump_hstr(struct rtp_head *head) {
  int i, j;

  fprintf(stdout, "------ Header Fields ---------\n");

  fprintf(stdout, "ptype   = %d\n",  head->ptype);
  fprintf(stdout, "pfields = %d\n",  head->pfields);
  fprintf(stdout, "pmin    = %g\n",  head->pmin);
  fprintf(stdout, "pmax    = %g\n",  head->pmax);

  fprintf(stdout, "ngas    = %d\n",  head->ngas);

  fprintf(stdout, "glist  = ");
  for (j=0; j < head->ngas; j++) 
    fprintf(stdout, " %d", head->glist[j]);
  fprintf(stdout, "\n");

  fprintf(stdout, "gunit = ");
  for (j=0; j < head->ngas; j++) 
    fprintf(stdout, " %d", head->gunit[j]);
  fprintf(stdout, "\n");

  /* radiance fields */
  fprintf(stdout, "pltfid = %d\n", head->pltfid);
  fprintf(stdout, "instid = %d\n", head->instid);
  fprintf(stdout, "nchan  = %d\n", head->nchan);

  fprintf(stdout, "vchan  =");
  for (j=0; j < head->nchan; j++) 
    fprintf(stdout, " %g", head->vchan[j]);
  fprintf(stdout, "\n");

  fprintf(stdout, "ichan  =");
  for (j=0; j < head->nchan; j++) 
    fprintf(stdout, " %d", head->ichan[j]);
  fprintf(stdout, "\n");

  /* user defined fields */
  fprintf(stdout, "iudef =");
  for (j=0; j < MAXIUDEF; j++) 
    fprintf(stdout, " %d", head->iudef[j]);
  fprintf(stdout, "\n");
  fprintf(stdout, "itype = %d\n", head->itype);

  fprintf(stdout, "-----------------------------\n");
}

