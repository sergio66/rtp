
/* rtptest1 -- basic RTP write and read test
 *
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
  struct rtp_prof prof1[4], prof2[4];
  struct FLIST (*hflist)[], (*pflist)[];
  struct ALIST (*halist)[], (*palist)[];
  struct ALIST halist1[8], palist1[8];
  int npro;

  /* fill in some header values
   */
  headinit(&head1);
  npro = 2;

  head1.memis = 10;
  head1.mlevs  = 80;

  head1.ptype = LEVPRO;			/* level profile */
  head1.pfields = PROFBIT + IROBSVBIT;  /* profile + obs radiances */
  head1.pfields = PROFBIT + IROBSVBIT + IRCALCBIT;

  head1.pmin  = 0.01;
  head1.pmax  = 1100;

  head1.ngas  = 6;
  for (i=0; i < head1.ngas; i++) {
    head1.glist[i] = 2*i+1;
    head1.gunit[i] = 1;
  }

  head1.nchan = 200;
  for (i=0; i < head1.nchan; i++) {
    head1.vchan[i] = 500 + i/4.0;
    head1.ichan[i] = i;
  }

  /* fill in some profile values
   */
  for (k=0; k < npro; k++) {

    profinit(&prof1[k]);

    prof1[k].plat = 32;
    prof1[k].plon = 55;
    prof1[k].ptime = 31536000;

    prof1[k].stemp = 300+k;
    prof1[k].salti = 10;
    prof1[k].spres = 1000;
    prof1[k].efreq[0] = 700; prof1[k].efreq[1] = 1400;
    prof1[k].emis[0] = .95;  prof1[k].emis[1] = .96;
    prof1[k].rho[0] = .05;    prof1[k].rho[1] = .04;
    prof1[k].nlevs = head1.mlevs;

    for (j=0; j < head1.mlevs; j++) {
      prof1[k].plevs[j] = j * 10 + 1;
      prof1[k].ptemp[j] = 200 + j;
    }

    for (i=0; i < head1.ngas; i++) 
          for (j=0; j < head1.mlevs; j++)
	    prof1[k].gamnt[i][j] = (i+1) * 1000 + j + 1;

    prof1[k].ctype=0;
    prof1[k].cfrac=0.25;
    prof1[k].cprtop=800;

    prof1[k].ctype2=1;
    prof1[k].cfrac2=0.35;
    prof1[k].cprtop2=300;

    prof1[k].scanang = 42;
    prof1[k].satzen = 45;

    prof1[k].rtime = 31536000;

    for (i=0; i < head1.nchan; i++) {
      prof1[k].robs1[i] = .03;
      prof1[k].calflag[i] = 254;
    }
    prof1[k].robs1[head1.nchan-1] = .99;

    prof1[k].robsqual = 1;
    prof1[k].freqcal = -13.5;

    for (i=0; i < head1.nchan; i++) {
      prof1[k].rcalc[i] = .05;
    }
    prof1[k].rcalc[head1.nchan-1] = .99;

    prof1[k].udef[0] = 42; prof1[k].udef[1] = 43; 

    strcpy((char *) prof1[k].pnote, "rtptest1 comment string");
  }

  /* fill in some attribute data */
  halist1[0].fname = "header";
  halist1[0].aname = "title";
  halist1[0].atext = "sample headder attribute";
  halist1[1].fname = "ngas";
  halist1[1].aname = "units";
  halist1[1].atext = "(count)";

  palist1[0].fname = "plevs";
  palist1[0].aname = "units";
  palist1[0].atext = "millibars";
  palist1[1].fname = "gas_1";
  palist1[1].aname = "units";
  palist1[1].atext = "PPMV";

  dump_pstr(&head1, &prof1[0]);

  fprintf(stdout, "============ write test ===========\n");

  rtpwrite1("rtptest1.hdf", 
	    &head1,
	    &halist1, 2, 	    
	    &palist1, 2, 	    
	    &ci);

  rtpwrite2(ci, (char *) &prof1[0]);
  rtpwrite2(ci, (char *) &prof1[1]);

  rtpclose1(ci);

  fprintf(stdout, "============ read test ===========\n");

  rtpread1("rtptest1.hdf", 
	   &head2,
	   &halist, &hnattr,
	   &palist, &pnattr, 	    
	   &ci);

  rtpread2(ci, (char *) &prof2[0]);
  rtpread2(ci, (char *) &prof2[1]);

  dump_chan(ci);

  rtpclose1(ci);

  dump_attrs(halist, hnattr, "rtptest() header alist dump");
  dump_attrs(palist, pnattr, "rtptest() profile alist dump");

  dump_pstr(&head2, &prof2[0]);

  printf("head.mlevs = %d\n", head2.mlevs);
  printf("MAXCALF = %d\n", MAXCALF);
  printf("MAXPN4 = %d\n", MAXPN4);

  return(0);
}
