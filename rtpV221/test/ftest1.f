c
c simple RTP attribute demo
c
        program ftest
        implicit none
c
c RTP declarations
c
        include 'rtpdefs.f'
        integer rtpopen, rtpread, rtpwrite, rtpclose
        record /RTPHEAD/ head
        record /RTPPROF/ prof
        record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
c 
c other variables
c
        integer status, rchan, i
        character*18 fname, mode
        character*64 version
c
c return our version string
c
        version = ''
        call rtpvers(version)
        print *, version
c
c initialize header and profile fields with sensible values
c
       call rtpinit(head, prof)
c
c set some sample header attributes
c
c we set three general 'header' attributes and an attribute
c for the header field 'glist'
c
        hatt(1).fname = 'header'//char(0)
        hatt(1).aname = 'title'//char(0)
        hatt(1).atext = 'Simple Attribute Demo'//char(0)

        hatt(2).fname = 'header'//char(0)
        hatt(2).aname = 'author'//char(0)
        hatt(2).atext = 'Hobart Muffler'//char(0)

        hatt(3).fname = 'header'//char(0)
        hatt(3).aname = 'date'//char(0)
        hatt(3).atext = '23 Oct. 2008'//char(0)

        hatt(4).fname = 'glist'//char(0)
        hatt(4).aname = 'note'//char(0)
        hatt(4).atext = 'sample fast model constituent set'//char(0)

        hatt(5).fname = char(0)
c
c set some sample profile attributes
c
        patt(1).fname = 'ptemp'//char(0)
        patt(1).aname = 'units'//char(0)
        patt(1).atext = 'degrees K'//char(0)
        
        patt(2).fname = 'stemp'//char(0)
        patt(2).aname = 'comment'//char(0)
        patt(2).atext = 'simulated surface temp'//char(0)

        patt(3).fname = 'satzen'//char(0)
        patt(3).aname = 'comment'//char(0)
        patt(3).atext = 'satellite zenith angle'//char(0)

        patt(4).fname = char(0)
c
c just set a few demo values, for this test
c
        head.ptype = LEVPRO
        head.pfields = PROFBIT + IRCALCBIT
        head.ngas = 2
        head.glist(1) = 1
        head.glist(2) = 3
        head.nchan = 10

        head.memis = 0
        head.mlevs = 8

        prof.nlevs = 8
        prof.stemp = 290
        prof.scanang = 42

        prof.plevs(1) = 10
        prof.plevs(2) = 11
        prof.plevs(3) = 12
        prof.plevs(8) = 20

        prof.ptemp(1) = 280
        prof.ptemp(2) = 281
        prof.ptemp(3) = 282
        prof.ptemp(8) = 300

        prof.txover = 280
        prof.co2ppm = 410
        prof.ctype = 6

        prof.tcc = 0.87
        prof.cc(1) = 0.6
        prof.cc(2) = 0.61
        prof.cc(8) = 0.67
        prof.ctype = 42
        prof.cfrac = 0.3

        do i = 1, head.nchan
          prof.rcalc(i) = 0.05
        end do
        prof.rcalc(head.nchan) = 0.99

        prof.pnote = 'tuna from luna'//char(0)
        prof.udef(1) = 42
        prof.udef(2) = 43
c
c create the file, write the attributes, write the (mostly unfilled) 
c profile structure, and close the file
c
        mode = 'c'
        fname = 'ftest1.hdf'

        status = rtpopen(fname, mode, head, hatt, patt, rchan)
        print *, 'open status = ', status

        status = rtpwrite(rchan, prof)
        print *, 'write status = ', status

        status = rtpclose(rchan)
        print *, 'close status = ', status

        print *, 'MAXCALF = ', MAXCALF
        print *, 'MAXPN4 = ', MAXPN4 

        stop
        end
