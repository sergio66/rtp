
I've put a new draft of the RTP format specification, a trial
implementation of the Fortran and C RTP API, and some test data
on my AIRS work page, at http://asl.umbc.edu/pub/motteler .  There
is also a link to Matlab API there, where it is part of a set of
Matlab HDF tools.

The Fortran and C interface source is bundled in a tar file, the
README included there there gives info on building things and
pointers to further information.  This is just an "alpha" release;
it passes the tests included with it, but has not been extensively
tested.
  
While the fortran interface is (hopefully) relatively friendly,
the Fortran implementation is complicated in part because it tries
to be sensible about variable field sets and field sizes; this is
very convenient for the user, but the sort of thing that is much
easier from Matlab than Fortran.  A simpler implementation would
map Fortran structures directly to the HDF vdatas, but would require
recompilation for every format or field variation, would need to
include space for radiances even when they weren't needed, and so
on.



