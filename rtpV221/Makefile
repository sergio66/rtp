#
# RTP top level Makefile
#
# this file sets some common parameters and then calls makefiles
# in the src, test, and utils subdirectories
#
# parameters should normally be set here in preference to setting
# them in the subdirectories.  In most cases, all you need to set
# are HDFHOME, RTPHOME, and your local compiler options
#

# get HDF libs and Intel compilers
# module load HDF/4.2.14-GCCcore-8.3.0
# module load intel-compilers/2021.4.0

# --------------
# HDF parameters
# --------------

HDFHOME = /usr/ebuild/software/HDF/4.2.14-GCCcore-8.3.0
HDFINC = -I$(HDFHOME)/include/hdf
HDFLIB = -L$(HDFHOME)/lib -ldf

# --------------
# RTP parameters
# --------------

RTPHOME = ..

# ------------------
# C compiler options
# -------------------

CC = icc
# CFLAGS = -O2 -fpack-struct
CFLAGS = -O2

# CC = gcc
# CFLAGS = -O2

# ------------------------
# Fortran compiler options
# ------------------------

F77 = ifort
# FFLAGS = -extend-source 132 -check bounds -O2 -align rec8byte
FFLAGS = -extend-source 132 -check bounds -O2

# F77 = gfortran
# FFLAGS = -ffixed-line-length-none -fdec-structure

# -------------------------------------------
# pass options to the subdirectory makefiles
# -------------------------------------------

EXPORTS = RTPHOME="$(RTPHOME)" \
	HDFINC="$(HDFINC)" HDFLIB="$(HDFLIB)" \
	CC="$(CC)" CFLAGS="$(CFLAGS)" \
	F77="$(F77)" FFLAGS="$(FFLAGS)" \

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

