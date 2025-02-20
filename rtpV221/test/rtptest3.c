/*
 * rtptest3 -- big file and field test
 *
 * writes and then reads a large rtp file, and checks read values
 * for the fields plevs, ptemp, gamnt, robs1, and rcalc.
 */

#include <stdio.h>
#include <stdlib.h>
#include "hdf.h"
#include "rtp.h"
#include "pvfnx.h"
#include "rtpfnx.h"

main (int argc, char *argv[]) {

  double x, y, z;
  int i, j, k, s1, ci;
  int32 file_id;
  int32 vhead_id, vprof_id;
  int hnrec, hnfield, hnattr;
  int pnrec, pnfield, pnattr;
  struct rtp_head head1, head2;
  struct rtp_prof prof1, prof2;
  struct FLIST (*hflist)[], (*pflist)[];
  struct ALIST (*halist)[], (*palist)[];
  struct ALIST halist1[8], palist1[8];
  int npro, dpro;

  fprintf(stdout, "============ write test ===========\n");

  npro = 36000;   /* number of profiles */
  dpro = 5;       /* dump or check this profile */

  /* header test values */
  headinit(&head1);

  head1.memis = 10;
  head1.mlevs  = 102;

 head1.ptype = LEVPRO;
  head1.pfields = PROFBIT + IROBSVBIT + IRCALCBIT;
  /*  head1.pfields = PROFBIT + IROBSVBIT; */

  head1.pmin  = 0.01;
  head1.pmax  = 1100;

  head1.ngas  = 16;
  for (i=0; i < head1.ngas; i++) {
    head1.glist[i] = 2*i+1;
    head1.gunit[i] = 1;
  }

  head1.nchan = 4300;
  for (i=0; i < head1.nchan; i++) {
    head1.vchan[i] = 500 + i/4.0;
    head1.ichan[i] = i;
  }

  /* header test attributes  */
  halist1[0].fname = "header";
  halist1[0].aname = "title";
  halist1[0].atext = "sample header attribute";
  halist1[1].fname = "ngas";
  halist1[1].aname = "units";
  halist1[1].atext = "(count)";

  /* profile test attributes */
  palist1[0].fname = "plevs";
  palist1[0].aname = "units";
  palist1[0].atext = "millibars";
  palist1[1].fname = "gas_1";
  palist1[1].aname = "units";
  palist1[1].atext = "PPMV";

  rtpwrite1("rtptest3.hdf", 
	    &head1,
	    &halist1, 2, 	    
	    &palist1, 2, 	    
	    &ci);

  /* loop on profiles
   */
 for (k=0; k < npro; k++) {

    /* generate profile test values */
    profinit(&prof1);

    prof1.plat = 32;
    prof1.plon = 55;
    prof1.ptime = 1000 + k;

    prof1.stemp = 280;
    prof1.salti = 10;
    prof1.spres = 1000;
    prof1.efreq[0] = 700; prof1.efreq[1] = 1400;
    prof1.emis[0] = .95;  prof1.emis[1] = .96;
    prof1.rho[0] = .05;   prof1.rho[1] = .04;
    prof1.nlevs = head1.mlevs;

    for (j=0; j < head1.mlevs; j++) {
      prof1.plevs[j] = j * 10 + 1;
      prof1.ptemp[j] = 200 + j;
    }

    for (i=0; i < head1.ngas; i++) 
          for (j=0; j < head1.mlevs; j++)
	    prof1.gamnt[i][j] = (i+1) * 100 + j + 1;

    prof1.ctype=0;
    prof1.cfrac=0.25;
    prof1.cprtop=800;

    prof1.ctype2=1;
    prof1.cfrac2=0.35;
    prof1.cprtop2=300;

    prof1.scanang = 42;
    prof1.satzen = 45;

    prof1.rtime = 1000 + k;

    for (i=0; i < head1.nchan; i++) {
      prof1.robs1[i] = i % 97;
      prof1.calflag[i] = 254;
    }

    prof1.robsqual = 1;
    prof1.freqcal = -13.5;

    for (i=0; i < head1.nchan; i++) {
      prof1.rcalc[i] = i % 101;
    }

    strcpy((char *) prof1.pnote, "rtptest3 comment string");

    prof1.udef[0] = 42; prof1.udef[1] = 43; 

    /*
    if (k == dpro)
      dump_pstr(&head1, &prof1);
    */
    /* add a deliberate error, see if we catch it */
    /* if (k == dpro)
      /* prof1.robs1[head1.nchan/2] = 42; */
      /* prof1.ptemp[20] = 42; */
      /* prof1.gamnt[2][20] = 42; */

    rtpwrite2(ci, (char *) &prof1);

  } /* end profile write loop */

  /* dump_chan(ci); */
  rtpclose1(ci);

  fprintf(stdout, "============ read test ===========\n");

  headinit(&head2);

  rtpread1("rtptest3.hdf", 
	   &head2,
	   &halist, &hnattr,
	   &palist, &pnattr, 	    
	   &ci);

  /* loop on profiles 
   */
  for (k=0; k < npro; k++) {

    profinit(&prof2);
    rtpread2(ci, (char *) &prof2);

    /* check a few key values */
    for (j=0; j < head2.mlevs; j++)
      if ((prof2.plevs[j] != j * 10 + 1 ) ||
          (prof2.ptemp[j] != 200 + j)) {
        fprintf(stdout, "plevs or ptemp error prof %d lev %d\n", k, j);
          exit(-1);
      }
    for (i=0; i < head2.ngas; i++) 
      for (j=0; j < head2.mlevs; j++)
        if (prof2.gamnt[i][j] != (i+1) * 100 + j + 1) {
          fprintf(stdout, "gamnt error prof %d ind %d lev %d\n", k, i, j);
          exit(-1);
        }
    for (i=0; i < head2.nchan; i++)
      if (prof2.robs1[i] != i % 97) {
        fprintf(stdout, "robs1 error prof %d chan %d\n", k, i);
        exit(-1);
      }
    for (i=0; i < head2.nchan; i++)
      if (prof2.rcalc[i] != i % 101) {
        fprintf(stdout, "rcalc error prof %d chan %d\n", k, i);
        exit(-1);
      }
    /*
    if (k == dpro)
      dump_pstr(&head2, &prof2);
    */
  }

  /* dump_chan(ci); */
  rtpclose1(ci);

  /*
  dump_attrs(halist, hnattr, "rtptest() header alist dump");
  dump_attrs(palist, pnattr, "rtptest() profile alist dump");
  */

  printf("head.mlevs = %d\n", head2.mlevs);
  printf("MAXCALF = %d\n", MAXCALF);
  printf("MAXPN4 = %d\n", MAXPN4);

  return(0);
}
