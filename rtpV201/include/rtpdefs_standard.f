
c
c     RTP Fortran API structures and parameters
c     Version 2.01
c     The Fortran structures defined here, RTPHEAD, RTPPROF, and
c     RTPATTR, must match the corresponding C structures rtp_head,
c     rtp_prof, and rtpfatt, and the parameters set below must have
c     the same values as the corresponding C #define parameters.
c
c     Note: total record size for header or profile records
c     may not exceed 50 kB
c
c     See rtp.h and rtpspec.pdf for more information on the fields
c     defined below.
c

c --------------
c RTP parameters
c --------------

        integer BAD         ! value to be used if no data
        integer LEVPRO      ! levels-type profile flag
        integer LAYPRO      ! layers-type profile flag
        integer AIRSLAY     ! AIRS layers-type profile flag

        integer PROFBIT     ! Summed radiance bit flags
        integer IRCALCBIT   ! bit flag for calculated radiance
        integer IROBSVBIT   ! bit flag for observed radiance
        integer PFIELDSMAX  ! max allowed value of PROFBIT

        integer MAXEMIS     ! max num of emissivity/rho points
        integer MAXGAS      ! max number of gases
        integer MAXGASID    ! max gas ID value
        integer MAXLEV      ! max number of levels
        integer MAXCHAN     ! max number of channels
        integer MAXPNOTE    ! max profile comment string
        integer MAXUDEF     ! max profile udef values
        integer MAXIUDEF    ! max profile and header iudef values

        integer MAXOPEN     ! max num of open RTP files
        integer MAXNATTR    ! max number of attributes
        integer MAXVNAME    ! max num of chars in field or vdata name
        integer MAXANAME    ! max num of chars in attribute name
        integer MAXATEXT    ! max num of chars in attribute text

        integer MAXCALF     ! max cal flag bytes, ceil(MAXCHAN/4)*4
        integer MAXPN4      ! true max for pnote, ceil(MAXPNOTE/4)*4

        ! the following parameters must match the values set in rtp.h
        !
        parameter ( BAD = -9999 )
        parameter ( LEVPRO  = 0 )
        parameter ( LAYPRO  = 1 )
        parameter ( AIRSLAY = 2 )

        parameter ( PROFBIT    = 1 )
        parameter ( IRCALCBIT  = 2 )
        parameter ( IROBSVBIT  = 4 )
        parameter ( PFIELDSMAX = 7 )

        parameter ( MAXEMIS   =  100 )
        parameter ( MAXGAS    =   16 )
        parameter ( MAXGASID  =  303 )
        parameter ( MAXLEV    =  120 )
        parameter ( MAXCHAN   = 4231 )
        parameter ( MAXPNOTE  =   80 )
        parameter ( MAXUDEF   =   20 )
        parameter ( MAXIUDEF   =  10 )

        parameter ( MAXOPEN   =    8 )
        parameter ( MAXNATTR  =   32 )

        parameter ( MAXCALF   = ((MAXCHAN-1)/4+1)*4 )
        parameter ( MAXPN4    = ((MAXPNOTE-1)/4+1)*4 )

        ! the following parameters must match the values set in pvdefs.h
        ! and also the field sizes in the RTPATTR structure, defined below
        !
        parameter ( MAXVNAME  =   64 )
        parameter ( MAXANAME  =   64 )
        parameter ( MAXATEXT  =  1024 )


c --------------------
c RTP header structure
c --------------------
c
        STRUCTURE /RTPHEAD/

	  ! profile data
          integer*4  ptype		 ! profile type
          integer*4  pfields		 ! profile field set
          real*4     pmin 		 ! min plevs value
          real*4     pmax 		 ! max plevs value
          integer*4  ngas 		 ! number of gases
          integer*4  glist(MAXGAS)	 ! constituent gas list
          integer*4  gunit(MAXGAS)	 ! constituent gas units
                
	  ! radiance data
          integer*4  pltfid              ! platform ID/code number
          integer*4  instid              ! instrument ID/code number
          integer*4  nchan		 ! number of channels
          integer*4  ichan(MAXCHAN)	 ! channel ID numbers
          real*4     vchan(MAXCHAN)	 ! channel center freq.
          real*4     vcmin 		 ! chan set min freq, including wings
          real*4     vcmax 		 ! chan set max freq, including wings

	  ! maxes for profile fields
          ! these fields are not saved explicitly in the HDF file
          integer*4  memis		 ! max number of emis/rho points
          integer*4  mlevs		 ! max number of pressure level

	  ! user defined fields
          integer*4  iudef(MAXIUDEF)	 ! user-defined integer array
          integer*4  itype       	 ! user-defined integer

        END STRUCTURE


c ---------------------
c RTP profile structure
c ---------------------
c  
        STRUCTURE /RTPPROF/

	  ! profile location/time
          real*4     plat                 ! profile latitude
          real*4     plon                 ! profile longitude 
          real*8     ptime                ! profile time

          ! surface data
          real*4     stemp                ! surface temperature
          real*4     salti                ! surface altitude
          real*4     spres                ! surface pressure
          real*4     landfrac             ! land fraction
          integer*4  landtype             ! land type code
          real*4     wspeed               ! wind speed
          integer*4  nemis                ! number of emis. pts
          real*4     efreq(MAXEMIS)       ! emissivity freq's
          real*4     emis(MAXEMIS)        ! surface emissivities
          real*4     rho(MAXEMIS)         ! surface reflectance

          ! atmospheric data
          integer*4  nlevs                ! number of press levels
          real*4     plevs(MAXLEV)        ! pressure levels
          real*4     palts(MAXLEV)        ! level altitudes
          real*4     ptemp(MAXLEV)        ! temperature profile
          real*4     gamnt(MAXLEV,MAXGAS) ! gas amounts
          real*4     gtotal(MAXGAS)       ! total column gas amount
          real*4     gxover(MAXGAS)       ! gas crossover press
          real*4     txover               ! temperature crossover press
          real*4     co2ppm               ! CO2 mixing ratio

          ! clear flag/code
          integer*4  clrflag              ! clear flag/code
                
          ! cloud1 data
          integer*4  ctype                ! cloud type code
          real*4     cfrac                ! cloud fraction 
          real*4     cemis(MAXEMIS)       ! cloud top emissivity
          real*4     crho(MAXEMIS)        ! cloud top reflectivity
          real*4     cprtop               ! cloud top pressure
          real*4     cprbot               ! cloud bottom pressure
          real*4     cngwat               ! cloud non-gas water
          real*4     cpsize               ! cloud particle size
          real*4     cstemp               ! cloud surface temperature

          ! cloud2 data
          integer*4  ctype2               ! cloud2 type code
          real*4     cfrac2               ! cloud2 fraction 
          real*4     cemis2(MAXEMIS)      ! cloud2 top emissivity
          real*4     crho2(MAXEMIS)       ! cloud2 top reflectivity
          real*4     cprtop2              ! cloud2 top pressure
          real*4     cprbot2              ! cloud2 bottom pressure
          real*4     cngwat2              ! cloud2 non-gas water
          real*4     cpsize2              ! cloud2 particle size
          real*4     cstemp2              ! cloud2 surface temperature
          real*4     cfrac12              ! cloud1+2 fraction

	  ! radiance orientation data
          real*4     pobs                 ! observation pressure
          real*4     zobs                 ! observation height
          integer*4  upwell               ! radiation direction
          real*4     scanang              ! scan angle
          real*4     satzen               ! satellite zenith angle
          real*4     satazi               ! satellite azimuth angle

          ! sun info
          real*4     solzen               ! sun zenith angle
          real*4     solazi               ! sun azimuth angle
          real*4     sundist              ! Earth-Sun distance
          real*4     glint                ! glint distance or flag

          ! observation location/time
          real*4     rlat                 ! radiance obs lat.
          real*4     rlon                 ! radiance obs lon.
C         integer*4  rfill                ! align rtime on 8 byte bndry
          real*8     rtime                ! radiance obs time

          ! observation indices
          integer*4  findex		  ! file (granule) index 
          integer*4  atrack		  ! along-track index
          integer*4  xtrack		  ! cross-track index
          integer*4  ifov                 ! field of view index

	  ! observed radiance data
          real*4     robs1(MAXCHAN)       ! obs radiance
          character*1 calflag(MAXCALF)    ! obs rad per chan calib/qual flags
          integer*4  robsqual             ! obs rad overall quality flag/code
          real*4     freqcal              ! frequency calibration

	  ! calculated radiance data
          real*4     rcalc(MAXCHAN)       ! calc radiance

	  ! user defined fields
          character*80  pnote             ! profile annotation, size MAXPN4
          real*4     udef(MAXUDEF)	  ! user-defined real array
          integer*4  iudef(MAXIUDEF)	  ! user-defined integer array
          integer*4  itype                ! user0defined integer

        END STRUCTURE


c -----------------------
c RTP attribute structure
c -----------------------
c
c fname is the name of the field the attribute is to be associated
c       with, 'header' for a general header attribute, or 'profiles'
c       for a general profile attribute.  Its size declaration should
c       be the same as the parameter MAXVNAME.
c
c aname is the attribute name, e.g., 'units' for a field attribute, 
c       or 'TITLE', for a general header attribute.  Its size should
c       also be the same as the parameter MAXANAME.
c
c atext is the attribute text, 'e.g., '48 Fitting Profiles' might be
c       the atext of the header 'TITLE' attribute.  Its size should be 
c       the same as the parameter MAXATEXT.
c
        STRUCTURE /RTPATTR/

          character*64 fname	! associated field name
          character*64 aname	! attribute name
          character*1024 atext	! attribute text

        END STRUCTURE
