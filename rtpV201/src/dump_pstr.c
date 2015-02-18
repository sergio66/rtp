
/* dump_pstr -- dump an RTP profile structure
 * Version 2.01
 * note: this routine is intended mainly for test purposes, only a
 * subset of the profile fields are printed out.
 *
 * size bounds from the header struct are used to restrict printing
 * to relevant values
 */

#include <stdio.h>
#include "rtp.h"

void dump_pstr(
	       struct rtp_head *head,
	       struct rtp_prof *prof
	       ) {
  int i, j, k;

  fprintf(stdout, "------ Profile Surface & Atmosphere Fields -----------\n");

  fprintf(stdout, "plat   = %g\n",  prof->plat);
  fprintf(stdout, "plon   = %g\n",  prof->plon);
  fprintf(stdout, "ptime  = %g\n",  prof->ptime);

  fprintf(stdout, "stemp  = %g\n",  prof->stemp);
  fprintf(stdout, "salti  = %g\n",  prof->salti);
  fprintf(stdout, "spres  = %g\n",  prof->spres);

  fprintf(stdout, "nemis  = %d\n",  prof->nemis);

  fprintf(stdout, "efreq  =");
  for (k=0; k < prof->nemis; k++)
    fprintf(stdout, " %g", prof->efreq[k]);
  fprintf(stdout, "\n");

  fprintf(stdout, "emis   =");
  for (k=0; k < prof->nemis; k++)
    fprintf(stdout, " %g", prof->emis[k]);
  fprintf(stdout, "\n");

  fprintf(stdout, "rho    =");
  for (k=0; k < prof->nemis; k++)
    fprintf(stdout, " %g", prof->rho[k]);
  fprintf(stdout, "\n");

  fprintf(stdout, "nlevs  = %d\n", prof->nlevs);

  fprintf(stdout, "plevs  =");
  for (k=0; k < prof->nlevs; k++)
    fprintf(stdout, " %g", prof->plevs[k]);
  fprintf(stdout, "\n");

  fprintf(stdout, "ptemp  =");
  for (k=0; k < prof->nlevs; k++)
    fprintf(stdout, " %g", prof->ptemp[k]);
  fprintf(stdout, "\n");

  for (j=0; j < head->ngas; j++) {
    fprintf(stdout, "gas %d  =", head->glist[j]);
    for (k=0; k < prof->nlevs; k++)
      fprintf(stdout, " %g", prof->gamnt[j][k]);
    fprintf(stdout, "\n");
  }

  fprintf(stdout, "gxover =");
  for (j=0; j < head->ngas; j++)
    fprintf(stdout, " %g", prof->gxover[j]);
  fprintf(stdout, "\n");

  fprintf(stdout, "scanang = %g\n", prof->scanang);
  fprintf(stdout, "satzen = %g\n", prof->satzen);

  fprintf(stdout, "------ Profile Radiance Fields --------------------\n");

  if (head->pfields & IRCALCBIT) {
    fprintf(stdout, "rcalc = ");
    for (k=0; k < head->nchan; k++)
      fprintf(stdout, " %g", prof->rcalc[k]);
    fprintf(stdout, "\n");
  }

  if (head->pfields & IROBSVBIT) {
    fprintf(stdout, "rlat = %g\n", prof->rlat);
    fprintf(stdout, "rlon = %g\n", prof->rlon);
    fprintf(stdout, "rtime = %g\n", prof->rtime);
  }

  if (head->pfields & IROBSVBIT) {
    fprintf(stdout, "robs1 = ");
    for (k=0; k < head->nchan; k++)
      fprintf(stdout, " %g", prof->robs1[k]);
    fprintf(stdout, "\n");
    fprintf(stdout, "calflag  = ");
    for (k=0; k < head->nchan; k++)
      fprintf(stdout, " %d", prof->calflag[k]);
    fprintf(stdout, "\n");
  }


  fprintf(stdout, "------ Profile User Defined Fields -------------------\n");

  fprintf(stdout, "pnote = %s\n", prof->pnote);
  fprintf(stdout, "udef =");
  for (j=0; j < MAXUDEF; j++) 
    fprintf(stdout, " %g", prof->udef[j]);
  fprintf(stdout, "\n");

  fprintf(stdout, "------------------------------------------------------\n");
}

