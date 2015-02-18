
/* RTP Fortran compatibility header and profile structures
 * Version 2.01
 * The C and Fortran structure definitions and the C field "flists"
 * giving string names to fields must all match exactly.  Also, the
 * #define constants given below should match those in the Fortran
 * parameter declarations.
 *
 * Note: total record size for header or profile may not exceed 50 kB
 *
 * H. Motteler
 * 5 Mar 02
 * Update: 20 October 2011, Scott Hannon - change prof.calflag and
 *    prof.pnote from DFNT_UCHAR8 and DFNT_CHAR8 respectively to
 *    DFNT_UINT8 for both.  The C declaration of pnote was changed
 *    from char8 ro to uchar8; calflag was already uchar8.  These
 *    changes are intended to avoid problems with MATLABs conversion
 *    of character data to 16 bits, and thus make these two fields
 *    easier to work with in MALTAB as generic byte data.  In
 *    theory, C uchar8 is consistent with FORTRAN character*1.
 */

#include "hdf.h"
#include "pvdefs.h"

/* set the library version 
 */
#define VERSION "RTP version 2.01, 20 Oct 2011"


/* --------------
 * RTP parameters 
 * --------------
 */

/* the JPL "BAD" value 
 */
#define BAD -9999

/* header field "ptype" defined values 
 */
#define LEVPRO   0
#define LAYPRO   1
#define AIRSLAY  2

/* heder field "pfields" bit flags for profile field groups 
 */
#define PROFBIT	    1   /* profiles                             */
#define IRCALCBIT   2   /* calculated radiances                 */
#define IROBSVBIT   4   /* observed radiances                   */
#define PFIELDSMAX  7   /* max value of pfields (1+2+4)		*/

/* RTP structure array limits
 *
 * These are max size limits for the static header and profile
 * structures defined below; the actual sizes of the corresponding
 * vdata fields are set by the relevant header structure fields 
 */
#define MAXEMIS   100	/* max number of emis/rho points 	*/
#define MAXGAS     80	/* max number of gases  		*/
#define MAXGASID  303	/* max valid HITRAN gas ID 		*/
#define MAXLEV    120	/* max number of levels or layers	*/
#define MAXCHAN   800	/* max number of chanels		*/
#define MAXPNOTE   80   /* max profile comment string		*/ 
#define MAXUDEF    20   /* max profile  udef values	        */
#define MAXIUDEF   10   /* max profile and header iudef values	*/

/* RTP API array bounds 
 */
#define MAXOPEN     8	/* max number of open RTP channels 	*/
#define MAXNATTR   32   /* max number of attributes per file	*/

/* derived parameters
 */
#define MAXCALF    (((MAXCHAN-1)/4+1)*4)
#define MAXPN4     (((MAXPNOTE-1)/4+1)*4)


/* ----------------
 * header structure
 * ----------------
 */

struct rtp_head {

  /* profile data (7 fields)
   */
  int32   ptype;   	  	/* 0 = level prof., 1 = layer prof.	*/
  int32   pfields;	  	/* field-set code 			*/

  float32 pmin;		  	/* lowest profile press level 		*/
  float32 pmax;		  	/* hightst profile press level		*/
  int32   ngas;   	  	/* number of explicit constituents 	*/
  int32   glist[MAXGAS];  	/* constituent gas list, ngas-vector	*/
  int32   gunit[MAXGAS];  	/* constituent gas units, ngas-vector 	*/

  /* radiance data (7 fields)
   */
  int32   pltfid;		/* platform (satellite) ID number  	*/
  int32   instid;		/* instrument ID number                	*/
  int32   nchan;		/* number of IR radiance channels	*/
  int32   ichan[MAXCHAN]; 	/* IR channel numbers, nchan-vector 	*/
  float32 vchan[MAXCHAN]; 	/* IR chan center freq's, nchan-vector	*/ 
  float32 vcmin;	  	/* chan set min freq (including wings)	*/
  float32 vcmax;	  	/* chan set max freq (including wings)	*/

  /* maxes for profile fields (2 fields)
   * these fields are not saved explicitly in the HDF file
   */
  int32   memis;          	/* max number of emissivity points	*/
  int32   mlevs;          	/* max number of pressure levels	*/

  /* user-defined fields (2 fields)
   */
  int32   iudef[MAXIUDEF];	/* user-defined integer array		*/ 
  int32   itype;        	/* user-defined integer                	*/ 
};

/* -----------------
 * header field list
 * ----------------- 
 */

#ifdef RTPDEF

/* define the header field list ("flist")
 *
 * sizes and HDF types must match structure fields exactly;
 */
struct FLIST hfield[] = {

  /* profile data (7 fields) */
  "ptype",   DFNT_INT32,    1,
  "pfields", DFNT_INT32,    1,
  "pmin",    DFNT_FLOAT32,  1,
  "pmax",    DFNT_FLOAT32,  1,
  "ngas",    DFNT_INT32,    1,
  "glist",   DFNT_INT32,    MAXGAS,
  "gunit",   DFNT_INT32,    MAXGAS,

  /* radiance data (7 fields) */
  "pltfid",  DFNT_INT32,    1,
  "instid",  DFNT_INT32,    1,
  "nchan",   DFNT_INT32,    1,
  "ichan",   DFNT_INT32,    MAXCHAN,
  "vchan",   DFNT_FLOAT32,  MAXCHAN,
  "vcmin",   DFNT_FLOAT32,  1,
  "vcmax",   DFNT_FLOAT32,  1,

  /* max profile size fields (2 fields) 
   * these fields are not saved explicitly in the HDF file */
  "memis",   DFNT_INT32,    1,
  "mlevs",   DFNT_INT32,    1,

  /* user-defined fields (2 fields) */
  "iudef",   DFNT_INT32,    MAXIUDEF,
  "itype",   DFNT_INT32,    1
};

#else

/* declare the header field list as an extern 
*/
extern struct FLIST hfield[];

#endif

/* specify the number of header fields 
 */
#define NHFIELD 18


/* -----------------
 * profile structure
 * -----------------
 */

struct rtp_prof {

  /* profile location data (3 fields) */
  float32 plat;   		/* profile latitude     	*/
  float32 plon;   		/* profile longitude    	*/
  float64 ptime;  		/* profile time         	*/

  /* surface data (10 fields) */
  float32 stemp;  		/* surface temperature  	*/
  float32 salti;  		/* surface altitude     	*/
  float32 spres;  		/* surface pressure     	*/
  float32 landfrac;		/* land fraction		*/
  int32   landtype;		/* land type code		*/
  float32 wspeed;      		/* wind speed			*/
  int32   nemis; 		/* number of emis. pts		*/
  float32 efreq[MAXEMIS];	/* emissivity freq's    	*/
  float32 emis[MAXEMIS];	/* surface emissivities 	*/
  float32 rho[MAXEMIS];  	/* surface reflectance  	*/

  /* atmospheric data (9 fields) */
  int32   nlevs; 		/* number of press levels	*/
  float32 plevs[MAXLEV];	/* pressure levels      	*/
  float32 palts[MAXLEV];  	/* level altitudes      	*/
  float32 ptemp[MAXLEV];  	/* temperature profile  	*/
  float32 gamnt[MAXGAS][MAXLEV]; /* gas amounts		        */
  float32 gtotal[MAXGAS];	/* total column gas amount	*/
  float32 gxover[MAXGAS];	/* gas crossover press	        */
  float32 txover;		/* temperature crossover press	*/	
  float32 co2ppm;		/* CO2 PPMV			*/

  /* clear flag (1 field) */
  int32   clrflag;		/* clear flag/code		*/

  /* cloud1 data (9 fields) */
  int32   ctype;                /* cloud type code		*/
  float32 cfrac; 	        /* cloud fraction       	*/
  float32 cemis[MAXEMIS]; 	/* cloud top emissivity 	*/
  float32 crho[MAXEMIS]; 	/* cloud top reflectivity 	*/
  float32 cprtop; 	        /* cloud top pressure   	*/
  float32 cprbot; 	        /* cloud bottom pressure   	*/
  float32 cngwat; 	        /* cloud non-gas water  	*/
  float32 cpsize; 	        /* cloud particle size   	*/
  float32 cstemp;               /* cloud surface temperature    */

  /* cloud2 data (10 fields) */
  int32   ctype2;               /* cloud2 type code		*/
  float32 cfrac2; 	        /* cloud2 fraction       	*/
  float32 cemis2[MAXEMIS]; 	/* cloud2 top emissivity 	*/
  float32 crho2[MAXEMIS]; 	/* cloud2 top reflectivity 	*/
  float32 cprtop2; 	        /* cloud2 top pressure   	*/
  float32 cprbot2; 	        /* cloud2 bottom pressure   	*/
  float32 cngwat2; 	        /* cloud2 non-gas water  	*/
  float32 cpsize2; 	        /* cloud2 particle size   	*/
  float32 cstemp2;              /* cloud2 surface temperature   */
  float32 cfrac12;              /* cloud1+2 fraction		*/

  /* radiance orientation data (6 fields) */
  float32 pobs;			/* observation pressure		*/
  float32 zobs;			/* observation height		*/
  int32   upwell;		/* radiation direction		*/
  float32 scanang;		/* IR scan angle		*/
  float32 satzen;		/* IR zenith angle 		*/
  float32 satazi;		/* sat azimuth angle		*/

  /* sun data (4 fields) */
  float32 solzen;		/* sun zenith angle 		*/
  float32 solazi;		/* sun azimuth angle		*/
  float32 sundist;		/* sun-Earth distance		*/
  float32 glint;		/* sun glint distance 		*/

  /* observed radiance data (7 fields) */
  float32 rlat;           	/* obs rad lat.         	*/
  float32 rlon;           	/* obs rad lon.         	*/
  /* int32   rfill; */		/* align rtime on 8 byte bndry	*/
  float64 rtime;          	/* radiance obs time    	*/

  int32   findex;		/* file (granule) index		*/
  int32   atrack;		/* along-track index		*/
  int32   xtrack;     		/* cross-track index		*/
  int32   ifov;     		/* field of view index		*/

  /* observed radiance data (4 fields) */
  float32 robs1[MAXCHAN]; 	/* observed radiance           	*/
  uchar8  calflag[MAXCALF];   	/* calibration flags		*/
  int32   robsqual;		/* quality code/flag		*/
  float32 freqcal;              /* frequency calibration        */

  /* calculated radiance data (1 field) */
  float32 rcalc[MAXCHAN]; 	/* calculated radiance         	*/

  /* user-defined fields (4 fields) */
  uchar8  pnote[MAXPN4];   	/* profile annotation		*/
  float32 udef[MAXUDEF];	/* user-defined real array      */ 
  int32   iudef[MAXIUDEF];	/* user-defined integer array 	*/
  int32   itype;                /* user-defined integer         */
};

/* ------------------
 * profile field list
 * ------------------ 
 */

#ifdef RTPDEF

/* define the profile field list ("flist")
 *
 * sizes and HDF types must match structure fields exactly;
 */
struct FLIST pfield[] = {

  /* profile location data (3 fields) */
  "plat",     DFNT_FLOAT32,  1,
  "plon",     DFNT_FLOAT32,  1,
  "ptime",    DFNT_FLOAT64,  1,

  /* surface data (10 fields) */
  "stemp",    DFNT_FLOAT32,  1,
  "salti",    DFNT_FLOAT32,  1,
  "spres",    DFNT_FLOAT32,  1,
  "landfrac", DFNT_FLOAT32,  1,
  "landtype", DFNT_INT32,    1,
  "wspeed",   DFNT_FLOAT32,  1,
  "nemis",    DFNT_INT32,    1,
  "efreq",    DFNT_FLOAT32,  MAXEMIS,
  "emis",     DFNT_FLOAT32,  MAXEMIS,
  "rho",      DFNT_FLOAT32,  MAXEMIS,

  /* atmospheric data (9 fields) */
  "nlevs",    DFNT_INT32,    1,
  "plevs",    DFNT_FLOAT32,  MAXLEV,
  "palts",    DFNT_FLOAT32,  MAXLEV,
  "ptemp",    DFNT_FLOAT32,  MAXLEV,
  "gamnt",    DFNT_FLOAT32,  MAXGAS*MAXLEV,
  "gtotal",   DFNT_FLOAT32,  MAXGAS,
  "gxover",   DFNT_FLOAT32,  MAXGAS,
  "txover",   DFNT_FLOAT32,  1,
  "co2ppm",   DFNT_FLOAT32,  1,

  /* clear flag (1 field) */
  "clrflag",  DFNT_INT32,    1,

  /* cloud1 data (9 fields) */
  "ctype",    DFNT_INT32,    1,
  "cfrac",    DFNT_FLOAT32,  1,
  "cemis",    DFNT_FLOAT32,  MAXEMIS,
  "crho",     DFNT_FLOAT32,  MAXEMIS,
  "cprtop",   DFNT_FLOAT32,  1,
  "cprbot",   DFNT_FLOAT32,  1,
  "cngwat",   DFNT_FLOAT32,  1,
  "cpsize",   DFNT_FLOAT32,  1,
  "cstemp",   DFNT_FLOAT32,  1,

  /* cloud2 data (10 fields) */
  "ctype2",   DFNT_INT32,    1,
  "cfrac2",   DFNT_FLOAT32,  1,
  "cemis2",   DFNT_FLOAT32,  MAXEMIS,
  "crho2",    DFNT_FLOAT32,  MAXEMIS,
  "cprtop2",  DFNT_FLOAT32,  1,
  "cprbot2",  DFNT_FLOAT32,  1,
  "cngwat2",  DFNT_FLOAT32,  1,
  "cpsize2",  DFNT_FLOAT32,  1,
  "cstemp2",  DFNT_FLOAT32,  1,
  "cfrac12",  DFNT_FLOAT32,  1,

  /* radiance orientation data (6 fields) */
  "pobs",     DFNT_FLOAT32,  1,
  "zobs",     DFNT_FLOAT32,  1,
  "upwell",   DFNT_INT32,    1,
  "scanang",  DFNT_FLOAT32,  1,
  "satzen",   DFNT_FLOAT32,  1,
  "satazi",   DFNT_FLOAT32,  1,

  /* sun data (4 fields) */
  "solzen",   DFNT_FLOAT32,  1,
  "solazi",   DFNT_FLOAT32,  1,
  "sundist",  DFNT_FLOAT32,  1,
  "glint",    DFNT_FLOAT32,  1,

  /* observed radiance data (7 fields) */
  "rlat",     DFNT_FLOAT32,  1,
  "rlon",     DFNT_FLOAT32,  1,
  /*"rfill",  DFNT_INT32,    1, */
  "rtime",    DFNT_FLOAT64,  1,
  "findex",   DFNT_INT32,    1,
  "atrack",   DFNT_INT32,    1,
  "xtrack",   DFNT_INT32,    1,
  "ifov",     DFNT_INT32,    1,

  /* observed radiance data (4 fields) */
  "robs1",    DFNT_FLOAT32,  MAXCHAN,
  "calflag",  DFNT_UINT8,    ((MAXCHAN-1)/4+1)*4,
  "robsqual", DFNT_INT32,    1,
  "freqcal",  DFNT_FLOAT32,  1,

  /* calculated radiance data (1 field) */
  "rcalc",    DFNT_FLOAT32,  MAXCHAN,

  /* user-defined fields (4 fields) */
  "pnote",    DFNT_UINT8,    ((MAXPNOTE-1)/4+1)*4,
  "udef",     DFNT_FLOAT32,  MAXUDEF,
  "iudef",    DFNT_INT32,    MAXIUDEF,
  "itype",    DFNT_INT32,    1
};

#else

/* declare the profile field list as an extern
 */
extern struct FLIST pfield[];

#endif

/* specify the number of profile fields 
 */
/* #define NPFIELD (sizeof(pfield) / (sizeof(char*) + sizeof(int)*2)) */
#define NPFIELD 68


/* ---------------------------
 * Fortran attribute structure
 * ---------------------------
 */

/* 
 * attribute structure (for fortran attributes)
 * this must match the fortran structure RTPATTR 
 */
struct rtpfatt {
  char 	fname[MAXVNAME];
  char	aname[MAXANAME];
  char	atext[MAXATEXT];
};


/* ---------------------
 * RTP channel structure
 * ---------------------
 */
struct rtp_chan {
  int mode;		    /* 1=hread, 2=hwrite, 3=pread, 4=pwrite */
  int32 file_id;	    /* HDF file ID */
  int32 vdata_id;	    /* HDF vdata header ID */
  int nvf;		    /* number of vdata fields */
  int *p1, *p2, *p3; 	    /* vdata pointer arrays */
  char *vbuf;		    /* vdata buffer */
  struct FLIST (*flist)[];  /* pointer field list (for testing) */
};

#ifdef RTPDEF

/* define the open channel structure
 */
struct rtp_chan chan[MAXOPEN];

/* oflag should be initialized with MAXOPEN zeros 
 */
int oflag[] = {0,0,0,0,0,0,0,0};

#else

/* declare the channel structures as externs
 */
extern struct rtp_chan chan[];
extern int oflag[];

#endif

