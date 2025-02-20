#
# RTP utilities Makefile
#
# this file builds some sample C and Fortran utilities
#
# the site-specific parameters HDFHOME, RTPHOME, and the C and 
# Fortran compiler options can be set in the top level makefile

# --------------
# RTP parameters
# --------------

# RTP definitions, relative to RTPHOME
RTPINC  = -I$(RTPHOME)/include		# RTP C and F77 includes
RTPLIB  = -L$(RTPHOME)/lib -lrtp	# RTP C and F77 library

# -------------
# Build Targets
# -------------

all: rtpinfo rtpdump install

# ----------------
# rtp dump utility
# ----------------

rtpdump: rtpdump.o
	$(CC) $(CFLAGS) -o rtpdump rtpdump.o $(RTPLIB) $(HDFLIB)
rtpdump.o: rtpdump.c
	$(CC) $(CFLAGS) -c rtpdump.c $(RTPINC) $(HDFINC)

# ------------------
# rtp subset utility
# ------------------

# SUBSRC = subset_rtp.f \
# 	rdinfo_subset.f rdlist_subset.f util_subset.f sublst.f
# SUBOBS = subset_rtp.o \
# 	rdinfo_subset.o rdlist_subset.o util_subset.o sublst.o
# 
# rtpsubset: $(SUBOBS)
# 	$(F77) $(FFLAGS) -o rtpsubset $(SUBOBS) $(RTPLIB) $(HDFLIB) $(FLIB)
# 
# $(SUBOBS): $(SUBSRC) incSUB.f
# 	$(F77) $(FFLAGS) -c $(SUBSRC) $(RTPINC)

# ------------------
# rtp info utility
# ------------------

INFOSRC = info_rtp.f  util_subset.f
INFOOBS = info_rtp.o  util_subset.o

rtpinfo: $(INFOOBS)
	$(F77) $(FFLAGS) -o rtpinfo $(INFOOBS) $(RTPLIB) $(HDFLIB) $(FLIB)

$(INFOOBS): $(INFOSRC)
	$(F77) $(FFLAGS) -c $(INFOSRC) $(RTPINC)

# -------
# install
# -------

install: rtpdump rtpinfo
	cp -a rtpdump rtpinfo $(RTPHOME)/bin

# ---------------
# cleanup
# ---------------

clean:
	rm *.o rtpdump rtpinfo 2> /dev/null || true

