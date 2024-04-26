!
!     RTP Fortran API structures and parameters
!     Version 2.01
!     The Fortran structures defined here, RTPHEAD, RTPPROF, and
!     RTPATTR, must match the corresponding C structures rtp_head,
!     rtp_prof, and rtpfatt, and the parameters set below must have
!     the same values as the corresponding C #define parameters.
!
!     Note: total record size for header or profile records
!     may not exceed 50 kB
!
!     See rtp.h and rtpspec.pdf for more information on the fields
!     defined below.
!

! --------------
! RTP parameters
! --------------

integer :: BAD              ! value to be used if no data                            
integer :: LEVPRO           ! levels-type profile flag
integer :: LAYPRO           ! layers-type profile flag                            
integer :: AIRSLAY          ! AIRS layers-type profile flag

integer :: PROFBIT          ! Summed radiance bit flags
integer :: IRCALCBIT        ! bit flag for calculated radiance
integer :: IROBSVBIT        ! bit flag for observed radiance
integer :: PFIELDSMAX       ! max allowed value of PROFBIT
                            
integer :: MAXEMIS          ! max num of emissivity/rho points
integer :: MAXGAS           ! max number of gases
integer :: MAXGASID         ! max gas ID value
integer :: MAXLEV           ! max number of levels
integer :: MAXCHAN          ! max number of channels
integer :: MAXPNOTE         ! max profile comment string
integer :: MAXUDEF          ! max profile udef values
integer :: MAXIUDEF         ! max profile and header iudef values
                            
integer :: MAXOPEN          ! max num of open RTP files
integer :: MAXNATTR         ! max number of attributes
integer :: MAXVNAME         ! max num of chars in field or vdata name
integer :: MAXANAME         ! max num of chars in attribute name
integer :: MAXATEXT         ! max num of chars in attribute text
                            
integer :: MAXCALF          ! max cal flag bytes, ceil(MAXCHAN/4)*4
integer :: MAXPN4           ! true max for pnote, ceil(MAXPNOTE/4)*4

! the following parameters must match the values set in rtp.h
!
parameter (BAD = - 9999)  
parameter (LEVPRO = 0)  
parameter (LAYPRO = 1)  

parameter (AIRSLAY = 2)  
parameter (PROFBIT = 1)  
parameter (IRCALCBIT = 2)  
parameter (IROBSVBIT = 4)  

parameter (PFIELDSMAX = 7)  
parameter (MAXEMIS = 100)  
parameter (MAXGAS = 16)  
parameter (MAXGASID = 303)  
parameter (MAXLEV = 120)  
parameter (MAXCHAN = 4231)  
parameter (MAXPNOTE = 80)  
parameter (MAXUDEF = 20)  
parameter (MAXIUDEF = 10)  

parameter (MAXOPEN = 8)  
parameter (MAXNATTR = 32)  

parameter (MAXCALF = ( (MAXCHAN - 1) / 4 + 1) * 4)  
parameter (MAXPN4 = ( (MAXPNOTE-1) / 4 + 1) * 4)  

! the following parameters must match the values set in pvdefs.h
! and also the field sizes in the RTPATTR structure, defined bel
!
parameter (MAXVNAME = 64)  
parameter (MAXANAME = 64)  
parameter (MAXATEXT = 1024)  

! --------------------
! RTP header structure
! --------------------
!

STRUCTURE / RTPHEAD /  

! profile data
integer :: ptype	     ! profile type	  
integer :: pfields	     ! profile field set	  
real :: pmin		     ! min plevs value
real :: pmax		     ! max plevs value
integer :: ngas		     ! number of gases
integer :: glist (MAXGAS)    ! constituent gas list
integer :: gunit (MAXGAS)    ! constituent gas units

! radiance data
                             
integer :: pltfid            ! platform ID/code number
integer :: instid            ! instrument ID/code number
integer :: nchan	     ! number of channels	  
integer :: ichan (MAXCHAN)   ! channel ID numbers
real :: vchan (MAXCHAN)      ! channel center freq.
real :: vcmin		     ! chan set min freq, including wings
real :: vcmax		     ! chan set max freq, including wings

! maxes for profile fields
! these fields are not saved explicitly in the HDF file
integer :: memis	     ! max number of emis/rho points	  
integer :: mlevs	     ! max number of pressure level

! user defined fields      
integer :: iudef (MAXIUDEF)  ! user-defined integer array 
integer :: itype	     ! user-defined integer

END STRUCTURE  

! ---------------------
! RTP profile structure
! ---------------------

STRUCTURE / RTPPROF /  
! profile location/time
                                          
real :: plat                              ! profile latitude
real :: plon                              ! profile longitude
real (8) :: ptime                         ! profile time

! surface data
real :: stemp                             ! surface temperature
real :: salti                             ! surface altitude
real :: spres                             ! surface pressure
real :: landfrac                          ! land fraction
integer :: landtype                       ! land type code
real :: wspeed                            ! wind speed
integer :: nemis                          ! number of emis. pts
real :: efreq (MAXEMIS)                   ! emissivity freq's
real :: emis (MAXEMIS)                    ! surface emissivities
real :: rho (MAXEMIS)                     ! surface reflectance

! atmospheric data
integer :: nlevs                          ! number of press levels
real :: plevs (MAXLEV)                    ! pressure levels
real :: palts (MAXLEV)                    ! level altitudes
real :: ptemp (MAXLEV)                    ! temperature profile
real :: cc(MAXLEV)                        ! cloud cover
real :: ciwc(MAXLEV)                      ! cloud ice water content
real :: clwc(MAXLEV)                      ! cloud liq water content	
real :: gamnt (MAXLEV, MAXGAS)            ! gas amounts
real :: gtotal (MAXGAS)                   ! total column gas amount
real :: gxover (MAXGAS)                   ! gas crossover press
real :: txover                            ! temperature crossover press
real :: co2ppm                            ! CO2 mixing ratio

! clear flag/code
integer :: clrflag                        ! clear flag/code

! cloud1 data                                          
integer :: ctype                          ! cloud type code
real :: cfrac                             ! cloud fraction
real :: cemis (MAXEMIS)                   ! cloud top emissivity
real :: crho (MAXEMIS)                    ! cloud top reflectivity
real :: cprtop                            ! cloud top pressure
real :: cprbot                            ! cloud bottom pressure
real :: cngwat                            ! cloud non-gas water
real :: cpsize                            ! cloud particle size
real :: cstemp                            ! cloud surface temperature

! cloud2 data
integer :: ctype2                          ! cloud type code
real :: cfrac2                             ! cloud fraction
real :: cemis2 (MAXEMIS)                   ! cloud top emissivity
real :: crho2 (MAXEMIS)                    ! cloud top reflectivity
real :: cprtop2                            ! cloud top pressure
real :: cprbot2                            ! cloud bottom pressure
real :: cngwat2                            ! cloud non-gas water
real :: cpsize2                            ! cloud particle size
real :: cstemp2                            ! cloud surface temperature
real :: cfrac12                            ! cloud overlap

! radiance orientation data
real :: pobs                              ! observation pressure
real :: zobs                              ! observation height
integer :: upwell                         ! radiation direction
real :: scanang                           ! scan angle
real :: satzen                            ! satellite zenith angle
real :: satazi                            ! satellite azimuth angle

! sun info
real :: solzen                            ! sun zenith angle
real :: solazi                            ! sun azimuth angle
real :: sundist                           ! Earth-Sun distance
real :: glint                             ! glint distance or flag

! observation location/time
real :: rlat                              ! radiance obs lat.
real :: rlon                              ! radiance obs lon.
!integer*4  rfill                         ! align rtime on 8 byte bndry
real (8) :: rtime                         ! radiance obs time

! observation indices                               
integer :: findex	               ! file (granule) index	                                
integer :: atrack	               ! along-track index	                                 
integer :: xtrack	               ! cross-track index	  
integer :: ifov                        ! field of view index

! observed radiance data                                          
real :: robs1 (MAXCHAN)                ! obs radiance
                                       
character (len=1) :: calflag (MAXCALF) ! obs rad per chan calib/qual                                           
integer :: robsqual                    ! obs rad overall quality flag
real :: freqcal                        ! frequency calibration

! calculated radiance data                                          
real :: rcalc (MAXCHAN)                ! calc radiance
 
! user defined fields                                          
character (len=80) :: pnote            ! profile annotation, size MAX
real :: udef (MAXUDEF) 	               ! user-defined real array                            
integer :: iudef (MAXIUDEF) 	       ! user-defined integer array
integer :: itype                       ! user-defined integer


END STRUCTURE  
! -----------------------
! RTP attribute structure
! -----------------------
!
! fname is the name of the field the attribute is to be associated
!       with, 'header' for a general header attribute, or 'profiles'
!       for a general profile attribute.  Its size declaration should
!       be the same as the parameter MAXVNAME.
!
! aname is the attribute name, e.g., 'units' for a field attribute,
!       or 'TITLE', for a general header attribute.  Its size should
!       also be the same as the parameter MAXANAME.
!
! atext is the attribute text, 'e.g., '48 Fitting Profiles' might be
!       the atext of the header 'TITLE' attribute.  Its size should be
!       the same as the parameter MAXATEXT.
!

STRUCTURE / RTPATTR /  

character (len=64) :: fname    ! associated field name
character (len=64) :: aname    ! attribute name	  	  
character (len=1024) :: atext  ! attribute text	  
end STRUCTURE
