
/* RTP definitions that reserve space
 *
 * included only by rtpopen.  The header and profile field lists,
 * sizes, and HDF types must match the structure fields in rtp.h
 * exactly.
 */

/* -----------------
 * header field list
 * ----------------- 
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

/* ------------------
 * profile field list
 * ------------------ 
 */

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

  /* cloud data (4 fields) */
  "tcc",      DFNT_FLOAT32,  1,
  "cc",       DFNT_FLOAT32,  MAXLEV, 
  "ciwc",     DFNT_FLOAT32,  MAXLEV, 
  "clwc",     DFNT_FLOAT32,  MAXLEV, 

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

  /* observed radiance data (8 fields) */
  "rlat",     DFNT_FLOAT32,  1,
  "rlon",     DFNT_FLOAT32,  1,
  "rfill",    DFNT_INT32,    1,
  "rtime",    DFNT_FLOAT64,  1,
  "findex",   DFNT_INT32,    1,
  "atrack",   DFNT_INT32,    1,
  "xtrack",   DFNT_INT32,    1,
  "ifov",     DFNT_INT32,    1,

  /* observed radiance data (4 fields) */
  "robs1",    DFNT_FLOAT32,  MAXCHAN,
  "calflag",  DFNT_UINT8,    MAXCALF,
  "robsqual", DFNT_INT32,    1,
  "freqcal",  DFNT_FLOAT32,  1,

  /* calculated radiance data (1 field) */
  "rcalc",    DFNT_FLOAT32,  MAXCHAN,

  /* user-defined fields (4 fields) */
  "pnote",    DFNT_UINT8,    MAXPN4,
  "udef",     DFNT_FLOAT32,  MAXUDEF,
  "iudef",    DFNT_INT32,    MAXIUDEF,
  "itype",    DFNT_INT32,    1
};

/* define the open channel structure
 */
struct rtp_chan chan[MAXOPEN];

/* oflag should be initialized with MAXOPEN zeros 
 */
int oflag[] = {0,0,0,0,0,0,0,0};

