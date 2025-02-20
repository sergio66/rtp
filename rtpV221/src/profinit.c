
/* initialize an RTP profile structure
 */

#include "rtp.h"

void profinit(struct rtp_prof *prof) {

  int j;

   /* surface data */
  prof->plat            = BAD;
  prof->plon            = BAD;
  prof->ptime           = BAD;

  prof->stemp           = BAD;
  prof->salti           = BAD;
  prof->spres           = BAD;
  prof->landfrac        = BAD;
  prof->landtype        = BAD;
  prof->wspeed          = BAD;
  prof->nemis           = 0;
  prof->efreq[0]        = BAD;
  prof->emis[0]         = BAD;
  prof->rho[0]          = BAD;

  /* atmospheric data */
  prof->nlevs           = 0;
  prof->plevs[0]        = BAD;
  prof->palts[0]        = BAD;
  prof->ptemp[0]        = BAD;
  prof->gamnt[0][0]     = BAD;
  prof->gtotal[0]	= BAD;
  prof->gxover[0]	= BAD;
  prof->txover		= BAD;
  prof->co2ppm          = BAD;

  /* clear flag data */
  prof->clrflag         = BAD;

  /* cloud data */
  prof->tcc             = BAD;
  prof->cc[0]           = BAD;
  prof->ciwc[0]         = BAD;
  prof->clwc[0]         = BAD;

  /* cloud1 data */
  prof->ctype           = BAD;
  prof->cfrac           = BAD;
  prof->cemis[0]        = BAD;
  prof->crho[0]         = BAD;
  prof->cprtop          = BAD;
  prof->cprbot          = BAD;
  prof->cngwat          = BAD;
  prof->cpsize          = BAD;
  prof->cstemp          = BAD;

  /* cloud2 data */
  prof->ctype2          = BAD;
  prof->cfrac2          = BAD;
  prof->cemis2[0]       = BAD;
  prof->crho2[0]        = BAD;
  prof->cprtop2         = BAD;
  prof->cprbot2         = BAD;
  prof->cngwat2         = BAD;
  prof->cpsize2         = BAD;
  prof->cstemp2         = BAD;
  prof->cfrac12         = BAD;

  /* common radiance data */
  prof->pobs            = BAD;
  prof->zobs            = BAD;
  prof->upwell          = BAD;
  prof->scanang         = BAD;
  prof->satzen          = BAD;
  prof->satazi          = BAD;

  prof->solzen          = BAD;
  prof->solazi          = BAD;
  prof->sundist         = BAD;
  prof->glint           = BAD;

  /* observed radiance data */
  prof->rlat            = BAD;
  prof->rlon            = BAD;
  prof->rtime           = BAD;

  prof->findex          = BAD;
  prof->atrack          = BAD;
  prof->xtrack          = BAD;
  prof->ifov            = BAD;

  /* observed radiance data */
  prof->robs1[0]        = BAD;
  for (j=0; j < MAXCHAN; j++)  prof->calflag[j] = (char) 0;
  prof->robsqual        = BAD;
  prof->freqcal         = BAD;

  /* calculated radiance data */
  prof->rcalc[0]        = BAD;

  /* user-defined fields */
  for (j=0; j < MAXPNOTE; j++) prof->pnote[j] = ' ';
  for (j=0; j < MAXUDEF; j++)  prof->udef[j]  = BAD;
  for (j=0; j < MAXIUDEF; j++) prof->iudef[j] = BAD;
  prof->itype                                 = BAD;

}

