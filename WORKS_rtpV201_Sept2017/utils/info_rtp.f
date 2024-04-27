C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              INFO
C
!F77====================================================================


!ROUTINE NAME: INFO


!ABSTRACT: Return some basic info about a RTP file


!CALL PROTOCOL: main program


!INPUT PARAMETERS: none


!OUTPUT PARAMETERS: none


!INPUT/OUTPUT PARAMETERS: none


!RETURN VALUES: main program


!PARENT(S): main program


!ROUTINES CALLED: see description


!FILES ACCESSED:
C unit5 : input filename
C unit6 : output info summary


!COMMON BLOCKS:
C    none


!DESCRIPTION:
C    INFO program is basically a driver for the routines.
C    rtpopen: open input RTP and read "head"
C    rtpread: read the first "prof" (only) from input RTP
C    rtpclose: close RTP file
C    Write an info summary to unit6 (ie the screen).
C

!ALGORITHM REFERENCES: see DESCRIPTION


!KNOWN BUGS AND LIMITATIONS:
C    None


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 06 Oct 2005 Scott Hannon      Created
C 22 Oct 2008 Scott Hannon      Minor update for RTP v2.0


!END====================================================================

C      =================================================================
       PROGRAM INFO
C      =================================================================


C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE


C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
      include 'rtpdefs.f'


C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
C      none


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      none (main program)


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
C
       CHARACTER*1 MODE        ! RTP file mode (r=read or w=write)
       CHARACTER*70 FIN        ! input RTP filename
       INTEGER I               ! generic loop variable
       INTEGER STATUS          ! status of RTP file open
       INTEGER IOPCI           ! input RTP file I/O channel
C
C      for N2BITS
       INTEGER*4 NUMBER
       LOGICAL LFLAGS(32)
       LOGICAL LPROF           ! file includes profile data?
       LOGICAL LRCAL          ! file includes calc'ed radiance?
       LOGICAL LROBS          ! file includes observed radiance?
C
C      RTP structures (see "rtpdefs.f")
       INTEGER rtpopen                  ! function rtpopen
       INTEGER rtpclose                 ! function rtpclose
       INTEGER rtpread                  ! function rtpread
       RECORD /RTPHEAD/ HEAD            ! header data
       RECORD /RTPATTR/ HATT(MAXNATTR)  ! header attributes
       RECORD /RTPATTR/ PATT(MAXNATTR)  ! profile attributes
       RECORD /RTPPROF/ PROF

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE
C***********************************************************************
C***********************************************************************
C
C      -----------------------------
C      Get filename
C      -----------------------------
       READ(5,1001) FIN
 1001  FORMAT(A70)

C      -------------------
C      Open RTP input file
C      -------------------
       MODE='r'
       STATUS=rtpopen(FIN, MODE, HEAD, HATT, PATT, IOPCI)
ccc
c       print *, 'read open status = ', STATUS
ccc
C
C      -------------
C      Check pfields
C      -------------
       NUMBER=HEAD.pfields
       CALL N2BITS(NUMBER, LFLAGS)
       LPROF =LFLAGS(1)  ! PROFBIT is bit1
       LRCAL=LFLAGS(2)   ! RCALCBIT is bit2
       LROBS=LFLAGS(3)   ! ROBSVBIT is bit3

C      ----------------------
C      Read the first profile
C      ----------------------
       STATUS=rtpread(IOPCI, PROF)
C
C      ----------
C      Close file
C      ----------
       STATUS=rtpclose(IOPCI)

C      ---------------------
C      Write an info summary
C      ---------------------
       WRITE(6,1001) FIN
       WRITE(6,2010) 'ptype    =', HEAD.ptype
       WRITE(6,2010) 'pfields  =', HEAD.pfields
       WRITE(6,2020) 'PROFBIT  =', LPROF
       WRITE(6,2020) 'RCALCBIT=',  LRCAL
       WRITE(6,2020) 'ROBSVBIT=',  LROBS
       WRITE(6,2010) 'nchan    =', HEAD.nchan
       WRITE(6,2010) 'ngas     =', HEAD.ngas
       WRITE(6,2030) 'glist    =', (HEAD.glist(I),I=1,HEAD.ngas)
C
       WRITE(6,2010) 'findex(1)=', PROF.findex
       WRITE(6,2040) 'rtime(1) =', PROF.rtime
C
 2010  FORMAT(A10,1X,I4)
 2020  FORMAT(A10,1X,L1)
 2030  FORMAT(A10,44(1X,I2))
 2040  FORMAT(A10,1X,1PD15.8)
C
       STOP
       END
