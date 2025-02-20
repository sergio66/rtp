! rtptest1 -- basic RTP write and read test
!
!

!      this takes in a TestPointProf text file and changes it to a RTP file
!      Sergio De Souza-Machado
!      very closely linked to the ftest1.f file that Howard Motteler wrote
!
! so eg      makeRTPfile.x fin=day_py4cats.atmX                      will produce  makeRTPfile.ip.rtp
! so eg      makeRTPfile.x fin=day_py4cats.atmX fout=junk.ip.rtp     will produce  junk.ip.rtp

      program makeRTPfile
      implicit none

! RTP declarations
      include 'rtpdefs.f90'
      integer rtpopen, rtpread, rtpwrite, rtpclose
      record /RTPHEAD/ head
      record /RTPPROF/ prof,profR,profR1
      record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)

! other variables
      integer status, IOPCI, IOPCO, iPreSet, iI, iP

! for vararg
      CHARACTER*80 FIN, FOUT, mode
      INTEGER IARGC, IOERR
      INTEGER J
      INTEGER NARGS
 
      CHARACTER*80 BUF
      CHARACTER*80 VAL
      CHARACTER*80 VAR

!************************************************************************
       IOPCI = -1
       IOPCO = -1
       print *,'IOPCI, IOPCO = ',IOPCI,IOPCO

! see /asl/packages/rtpV201/utils/rdinfo_subset.f
!      ------------
!      Set defaults
!      ------------
       IOERR = 6
       FIN =  'profile.txt'           ! input filename
       FOUT = 'makeRTPfile.ip.rtp'    ! output filename
       FOUT = 'profile_ip.rtp'        ! output filename

!      -----------------------------------------------------------------
!      Loop on program parameters
!      --------------------------
!      Determine the number of command-line arguments
       NARGS=IARGC()

       IF (NARGS .EQ. 0) THEN
          WRITE(IOERR,1010)
 1010     FORMAT('makeRTPfile.f90 must be run with at least 1 argument')
          WRITE(IOERR,1011)
 1011     FORMAT('   fin  = <filename> txt input file  {mandatory}','//'    &
                 '   fout = <filename> rtp output file {optinal, set to profile_ip.rtp}')
          STOP
       ENDIF

!      Loop over the command-line arguments
       DO iI = 1, NARGS

!         Pull out the ith argument
          CALL GETARG(iI, BUF)

!         Find the "=" character in the command-line argument string
          J=INDEX(BUF, '=')

          IF (J .NE. 0) THEN

!            Name of variable
             VAR = BUF(1:J-1)
             CALL UPCASE(VAR)

!            Specified value
             VAL = BUF(J+1:LEN(BUF))

!            Big "IF" to set parameters
!            ----------------------------
             IF (VAR(1:3) .EQ. 'FIN') THEN
                FIN=VAL
             ELSEIF (VAR(1:4) .EQ. 'FOUT') THEN
                FOUT=VAL
             ELSE
                WRITE(IOERR,1020) VAR
 1020           FORMAT('Unknown command-line argument: ',A6)
                STOP

             ENDIF

          ENDIF
       ENDDO  ! end of loop over command-line arguments
       write(*,'(A,A)') 'input  file name = ',fin
       write(*,'(A,A)') 'output file name = ',fout

!************************************************************************
      iPreset = -1
      IF (iPreset .gt. 0) THEN
        ! set only one header attribute
        hatt(1).fname = 'glist'//char(0)
        hatt(1).aname = 'note'//char(0)
        hatt(1).atext = 'translated from text sonde profiles'//char(0)
        ! set some sample profile attributes
        patt(1).fname = 'ptemp'//char(0)
        patt(1).aname = 'units'//char(0)
        patt(1).atext = 'degrees K'//char(0)
      ELSE
        CALL RTPINIT(head,prof)
      END IF

!      print *,FIN
!      print *,FOUT

! read in the point profile
      CALL ReadProfile(prof,head,FIN)
      print *,'head ngas  = ',head.ngas
      print *,'     glist = ',head.glist(1:head.ngas)
      print *,'     gunit = ',head.gunit(1:head.ngas)
      write(*,'(A,I3,2(F12.4),ES12.4)') 'prof  nlevs,stemp,spres,sum(WV) = ', &
        prof.nlevs,prof.stemp,prof.spres,sum(prof.gamnt(1:prof.nlevs,1))
      print *,'     rlat,rlon                 = ',prof.rlat,prof.rlon
      print *,'        iI      P(z)           T(z)           WV(z)'
      print *,'-----------------------------------------------------------'
      do iI = 1,prof.nlevs
        print *,'write',iI,prof.plevs(iI),prof.ptemp(iI),prof.gamnt(iI,1)
      end do
      print *,'-----------------------------------------------------------'

        ! set only one header attribute
        hatt(1).fname = 'glist'//char(0)
        hatt(1).aname = 'note'//char(0)
        hatt(1).atext = 'translated from text sonde profiles'//char(0)
        ! set some sample profile attributes
!        patt(1).fname = 'ptemp'//char(0)
!        patt(1).aname = 'units'//char(0)
!        patt(1).atext = 'degrees K'//char(0)

! create the file, write the attributes, write the (mostly unfilled) 
! profile structure, and close the file
! 'c'=create, 'r'=read : see eg /home/sergio/OTHERSTUFF/rtpV201/src/rtpopen.c
      mode = 'c'
      write(*,'(A,A)') 'Writing output in running dir .. output name = ',FOUT

!      IOPCO = 1
      status = rtpopen(FOUT, mode, head, hatt, patt, IOPCO)
      print *, 'write open status, IOPCO  = ', status,IOPCO
      status = rtpwrite2(IOPCO, prof)
      print *, 'write status = ', status,IOPCO
      status = rtpclose(IOPCO)
      print *, 'write close status = ', status,IOPCO

      print *,' '
      print *,' '
      print *,' '

!************************************************************************

      IP = +1
      if (iP > 0) THEN
        mode = 'r'
!       FOUT = 'pin_feb2002_sea_airs45deg_op.so2.latlon.const_emiss_0.9.rtp'
        write(*,'(A,A)') 'now reading in file ',FOUT
        status = rtpopen(FOUT, mode, head, hatt, patt, IOPCI)
        print *, 'read open status, IOPCI  = ', status,IOPCI

        iP = 0
        do while (.true.)
          status = rtpread(IOPCI, profR1)
          if (status .eq. -1) goto 22

          profR = profR1
          iP = iP + 1
          write(*,'(A,I3,2(F12.4),ES12.4)') 'xprof nlevs,xstemp,xspres,sum(WV) = ', &
           profR.nlevs,profR.stemp,profR.spres,sum(profR.gamnt(1:profR.nlevs,1))
          print *,'     rlat,rlon                 = ',profR.rlat,profR.rlon

!           do iI = 1,profR.nlevs
!             print *,'read',iI,profR.plevs(iI),profR.ptemp(iI),profR.gamnt(iI,1)
!           end do

        end do
  22    continue
        write(*,'(A,I3,A)') 'successfully read ',iP,' profiles .... here is P(z),T(z),WV(z) for last one'
        write(*,'(A,I3,2(F12.4),ES12.4)') 'xprof nlevs,xstemp,xspres,sum(WV) = ', & 
          profR.nlevs,profR.stemp,profR.spres,sum(profR.gamnt(1:profR.nlevs,1))
        print *,'     rlat,rlon                 = ',profR.rlat,profR.rlon
        print *,'        iI      P(z)           T(z)           WV(z)'
        print *,'-----------------------------------------------------------'
        do iI = 1,profR.nlevs
          print *,'read',iP,iI,profR.plevs(iI),profR.ptemp(iI),profR.gamnt(iI,1)
        end do
        print *,'-----------------------------------------------------------'

        status = rtpclose(IOPCI)
        print *, 'read close status = ', status
      end if

      stop
      end

!************************************************************************
! this subroutine prompts the user for the name of the point profile,
! opens it, reads in the info and saves it in the "prof" record

! the info in the levels profile is as follows 
! !!! first few lines are comments, starting with ! or % or #
! cType,nlevs,ngas    where cType = 'p' or 'P'
! gasID list
! lat lon
! spres stemp         -9999 if just set it to max pres, temp at max pres
! loop   ii = 1 : nlevs
!   rPress,rTemp,rDensity,(raGasAmts(iI),iI = 1,ngas)

      SUBROUTINE ReadProfile(prof,head,caFname)

      implicit none

      include 'rtpdefs.f90'

      record /RTPHEAD/ head
      record /RTPPROF/ prof

      INTEGER ngas,nlevs,glist(MaxGas),thetype,ispres,istemp
      REAL pmin,pmax
      CHARACTER*80 caFName
      CHARACTER*280 caStr
      INTEGER iIOUN2,iErrIO
      INTEGER iI,iZ
      REAL rLat,rLon,spres,stemp
      REAL rPress,rTemp,rDensity,raGasAmts(MAXGAS),rZ
      CHARACTER*1 cType

      pmin = 1.0e30
      pmax = -1.0e30

      iIOUN2 = 15

! now loop iNpath/iNumGases  times for each gas in the user supplied profile 
!      print *,caFName

      OPEN(UNIT=iIOun2,FILE=caFName,STATUS='OLD',FORM='FORMATTED',IOSTAT=iErrIO) 
      IF (iErrIO .NE. 0) THEN 
        WRITE(*,1070) iErrIO
      END IF 
 1070 FORMAT('Error number ',I5,' opening your profile')

      !!!!!!!! go thru the comments at the top of the file
 30   CONTINUE 
      READ (iIOUN2,5030,ERR=13,END=13) caStr 
      IF ((caStr(1:1) .EQ. '!') .OR. (caStr(1:1) .EQ. '%') .OR. (caStr(1:1) .EQ. '#')) THEN 
        GO TO 30 
      END IF

      !!read type of sonde profile (pressure or height), num levels, num gases
      READ (caStr,*,ERR=13,END=13) cType,nlevs,ngas 
      IF ((cType .NE. 'p') .AND. (cType .NE. 'P') .AND. (cType .NE. 'h') .AND. (cType .NE. 'H'))  THEN
        print *,'Error : makeRTPfile assumes pressure info in the file'
        CLOSE(iIOUN2)
        STOP
      END IF
      IF ((cType .EQ. 'p') .OR. (cType .EQ. 'P')) THEN
        thetype = 1
      ELSEIF ((cType .EQ. 'h') .OR. (cType .EQ. 'H')) THEN
        thetype = 2
      END IF
      
      IF (ngas .GT. MAXGAS) THEN
        print *,'Error : makeRTPfile can only handle ',MAXGAS,' gases'
        CLOSE(iIOUN2)
        STOP
      END IF

      IF (nlevs .GT. MAXLEV) THEN
        print *,'Error : makeRTPfile can only handle ',MAXLEV,' levels'
        CLOSE(iIOUN2)
        STOP
      END IF

      !! read the gasID list
      READ (iIOUN2,5030,ERR=13,END=13) caStr 
      READ (caStr,*,ERR=13,END=13) (glist(iI),iI = 1,ngas)
 
      !! read the latitude and longitude
      READ (iIOUN2,5030,ERR=13,END=13) caStr 
      READ (caStr,*,ERR=13,END=13) rLat,rLon

      ispres = -1
      istemp = -1
      !! read the spres and stemp
      READ (iIOUN2,5030,ERR=13,END=13) caStr 
      READ (caStr,*,ERR=13,END=13) spres,stemp
      if (spres > 0) ispres = +1
      if (stemp > 0) istemp = +1
        
      prof.rlat = rLat
      prof.rlon = rLon
      prof.nlevs = nlevs

      !!read the profile info
      DO iZ = 1,nlevs
        READ (iIOUN2,5030,ERR=13,END=13) caStr
        IF (thetype .EQ. 1) THEN
          READ (caStr,*,ERR=13,END=13) rPress,rTemp,rDensity,(raGasAmts(iI),iI = 1,ngas)
!          print *,rPress,rTemp,rDensity,(raGasAmts(iI),iI = 1,ngas)        
        ELSEIF (thetype .EQ. 2) THEN
          READ (caStr,*,ERR=13,END=13) rZ,rPress,rTemp,(raGasAmts(iI),iI = 1,ngas)
!          print *,rZ,rPress,rTemp,(raGasAmts(iI),iI = 1,ngas)        
        END IF
        
        prof.plevs(iZ) = rPress
        prof.ptemp(iZ) = rTemp
        IF (thetype .EQ. 2) THEN
          prof.palts(iZ) = rZ
        END IF

        DO iI = 1,ngas
          prof.gamnt(iZ,iI) = raGasAmts(iI)
        END DO
        IF (rPress .LT. pmin) THEN
          pmin = rPress
        END IF
        IF (rPress .GT. pmax) THEN
          pmax = rPress
          IF (ispres < 0) spres = rPress
          IF (istemp < 0) stemp = rTemp
        END IF

      ENDDO

 13   CONTINUE
      CLOSE(iIOUN2)
 5030 FORMAT(A280)

      prof.spres = spres
      prof.stemp = stemp

      prof.scanang = 0.0
      prof.satzen  = 0.0
      prof.solzen  = 0.0


! set HEAD values
      head.ptype = LEVPRO
      head.pfields = PROFBIT
!      head.mrho = 0
      head.memis = 0
      head.ngas = ngas
! see /asl/packages/klayersV205/Doc/gas_units_code.txt     10 = ppmv     
      DO iI = 1,ngas
        head.glist(iI) = glist(iI)
        head.gunit(iI) = 10
      END DO

      prof.nlevs = nlevs

      head.nchan = 1
      head.ichan = 1520
      head.vchan = 1231.22
      head.pmin = pmin
      head.pmax = pmax
!      head.iudef = -9999
!      head.itype = -9999

      print*,'<----------- SUMMARY of INPUT TXT file -----------> '
      print *,' the gasIDs are : '
      print *, (glist(iI),iI = 1,ngas)
      print *,'  rLat,rLon = ', prof.rLat,prof.rLon
      print *,'  number of levels, number of gases =', prof.nlevs,head.ngas
      print *,'  surf pres, surf temp = ', prof.spres,prof.stemp
      print*,'<----------- SUMMARY of INPUT TXT file -----------> '
      print *,'  '

      RETURN
      END
!************************************************************************
       SUBROUTINE UPCASE(BUF)

!      Convert a string to upper case

       CHARACTER*(*) BUF

       INTEGER I
       INTEGER IC
       INTEGER J

       DO I=1,LEN(BUF)
          IC=ICHAR( BUF(I:I) )
          IF (IC .GE. ICHAR('a') .AND. IC .LE. ICHAR('z')) THEN
             J=IC + (ICHAR('A') - ICHAR('a'))
             BUF(I:I)=CHAR(J)
          ENDIF
       ENDDO

       RETURN
       END

!----------------------------------------------------------------------

  program rtptest1
  implicit none

  define RTPDEF
  include "rtp.h"

  double x, y, z
  int i, j, k, s1, ci
  int32 file_id
  int32 vhead_id, vprof_id
  int hnrec, hnfield, hnattr
  int pnrec, pnfield, pnattr
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

  head1.memis = 2;
  head1.mlevs  = 4;

  head1.ptype = LEVPRO;			/* level profile */
  head1.pfields = PROFBIT + IROBSVBIT;  /* profile + obs radiances */

  head1.ngas  = 2;
  head1.pmin  = 10;
  head1.pmax  = 1000;
  head1.glist[0] = 1; head1.glist[1] = 3;

  head1.nchan = 3;
  head1.vchan[0] = 700; head1.vchan[1] = 701; head1.vchan[2] = 702;
  head1.ichan[0] = 10; head1.ichan[1] = 11; head1.ichan[2] = 12;

  /* fill in some profile values
   */
  for (k=0; k < npro; k++) {

    profinit(&prof1[k]);

    prof1[k].plat = 32;
    prof1[k].plon = 55;
    prof1[k].ptime = 30001;

    prof1[k].stemp = 300+k;
    prof1[k].salti = 10;
    prof1[k].spres = 1000;
    prof1[k].efreq[0] = 700; prof1[k].efreq[1] = 1400;
    prof1[k].emis[0] = .95;  prof1[k].emis[1] = .96;
    prof1[k].rho[0] = .5;    prof1[k].rho[1] = .6;
    prof1[k].nlevs = 3;

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

    prof1[k].rtime = 30099;

    for (i=0; i < head1.nchan; i++) {
      prof1[k].robs1[i] = .02;
      prof1[k].calflag[i] = 254;
    }
    prof1[k].robsqual = 1;
    prof1[k].freqcal = -13.5;

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

  rtpwrite2(ci, &prof1[0]);
  rtpwrite2(ci, &prof1[1]);

  rtpclose1(ci);

  fprintf(stdout, "============ read test ===========\n");

  rtpread1("rtptest1.hdf", 
	   &head2,
	   &halist, &hnattr,
	   &palist, &pnattr, 	    
	   &ci);

  rtpread2(ci, &prof2[0]);
  rtpread2(ci, &prof2[1]);

  dump_chan(ci);

  rtpclose1(ci);

  dump_attrs(halist, hnattr, "rtptest() header alist dump");
  dump_attrs(palist, pnattr, "rtptest() profile alist dump");

  dump_pstr(&head2, &prof2[0]);

  printf("head.mlevs = %d\n", head2.mlevs);

  return(0);
}
