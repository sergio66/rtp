#
# RTP top level Makefile for rtpV201
#
# this file sets some common parameters and then calls makefiles
# in the src, test, and utils subdirectories
#
# parameters should normally be set here in preference to setting
# them in the subdirectories.  In most cases, all you need to set
# are HDFHOME, RTPHOME, and your local compiler options
#

# --------------
# HDF parameters
# --------------

# set HDFHOME to point to the local HDF installation
HDFHOME = /asl/opt/hdf/4.2r4
#HDFHOME = /asl/opt/lib_old

# --------------
# RTP parameters
# --------------

# set RTPHOME to point to the local RTP installation
RTPHOME = ..

# ------------------
# C compiler options
# -------------------

# -64 for 64-bit IRIX
# CFLAGS = -g -64
# CFLAGS = -g
CFLAGS = -O 
CC = icc    # F77 = ifort "icc does not adds underscores to symbols"
#CC = cc       # F77 = /asl/opt/absoft/absoft10.0/bin/af77

# ------------------------
# Fortran compiler options
# ------------------------

# Absoft Fortran
# --------------
# -N109  fold all names to upper case
# -C     check array bounds
# -O     some optimizations
# -N3    add record info to unformatted files
# -s     static allocation
# -N15   so that LAPACK lib can be included; adds trailing underscores
# FFLAGS = -N109 -C -O -N3 -s
# FFLAGS = -g
#FFLAGS = -O
#FLIB = -lU77	# Linux Absoft needs -lU77, # comment this out if using ifort
#F77 = /asl/opt/absoft/absoft10.0/bin/af77   # comment this out if using ifort

# SunOS options
# -------------
# FFLAGS = -e -fast -w     
# F77 = f77

# SGI options
# -----------
# include -cpp option if not default 
# FFLAGS = -O
# -64 for 64-bit IRIX
# FFLAGS  = -O -64
# F77 = f77

# ifort options
# -----------
# -us : puts trailing underscore
F77 = /cm/shared/apps/intel/composer_xe/2015.3.187/bin/intel64/ifort
F77 = /cm/shared/apps/intel/composer_xe/2015.5.223/bin/intel64/ifort
#FLAGS = -O3 -names uppercase -us
FLAGS = -O -assume nounderscore

# pass the variables set above to the subdirectory makefiles
#
EXPORTS = HDFHOME="$(HDFHOME)" RTPHOME="$(RTPHOME)" \
	CC="$(CC)" CFLAGS="$(CFLAGS)" CLIB="$(CLIB)" \
	F77="$(F77)" FFLAGS="$(FFLAGS)" FLIB="$(FLIB)"

# -------------
# Make Targets
# -------------

all: SRC TEST UTILS

SRC:
	cd src && make $(EXPORTS)
TEST:
	cd test && make $(EXPORTS)
UTILS:
	cd utils && make $(EXPORTS)

clean:
	cd src && make clean
	cd test && make clean
	cd utils && make clean

# -------------------------
# make an RTP distribution
# -------------------------
#
# "make dist" makes a distribution named by the current working
# directory.  For example, if we are in the subdirectory rtpV201
# "make dist" will clean things up and then create an rtpV201.tar 
# in the parent directory that unpacks to rtpV201/<etc>.
# 
dist: clean
	rm rtp.tar bin/* lib/* 2> /dev/null || true
	rbase=`/bin/pwd`                && \
	    rbase=`basename $${rbase}`  && \
	    cd ..                       && \
	    tar -cf $${rbase}.tar          \
		$${rbase}/bin              \
		$${rbase}/doc              \
		$${rbase}/include          \
		$${rbase}/lib              \
		$${rbase}/Makefile         \
		$${rbase}/README           \
		$${rbase}/src              \
		$${rbase}/test             \
		$${rbase}/utils
	@echo created `/bin/pwd`.tar
	@echo "\"make all\" to rebuild the local distribution"

