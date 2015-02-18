
I did some new profiling tests with the latest versions of
rtpread.m and rtpwrite.m, with Matlab 6.0, and they turned out
somewhat differently from earlier tests.

The new versions of the reader and writer are doing more field
manipulation and checking; for example, a lot of time is taken
moving between the "gamnt" 2D gas array, and the 1D "gas_n" fields
needed by the vdatas, which even without the "rmfield" involves a
lot of structure manipulation.  This structure manipulation is the
bottleneck now, with maybe 10% or less (depending on data set and
record sizes) of the time spend in the lower-level HDF read or
write calls.

The HDF calls seem a little faster in Matlab 6.0, compared to the
prerelease, but I can't check that accurately as I deleted the
prerelease, and have been changing the Matlab code code.

While things could certainly be faster, after commenting out the
"rmfields" I'm not seeing the very slow times you were getting;
using the 48 fitting profiles as a test, I was reading about 12
profiles/second, and writing about 24 profiles/second.  This could
be partly because there's less in the profiles, but I think you
were running into a problem with the code slowing down with big
profile sets.

What might be happening is some of the Matlab field operations, or
maybe the HDF calls, are O(n^2) in the number of structure records,
so that things just grind to a halt with big data sets.  The fix
for that is to make the Matlab code closer to the Fortran, with one
read or write per call, rather than doing the whole array at once,
and I'd intended to do that anyway so we could more easily process
really big data sets.


On a related matter, I am worried about JPL trying to write their
own readers; while the RTP format is based on plain HDF 4 vdatas,
the varying field sets might cause them problems.  If a main goal
was to allow for a simple reader and writer, I would have done
something a lot simpler.

(Note that the first draft of the RTP format, along with a Matlab
implementation, was released way back last Nov 9, and the first
C/Fortran implementation Jan 23, so it has been out for JPL to
look at for some time now, so they have had time to look at it if
they were so inclined.)

The simplest possible way to implement RTP would be to just read
and write out the Fortran structures and attribute arrays directly
out to a file.  This is very tempting, but there are some tradeoffs
and drawbacks.

Benefits of such a simple implementation include:

  - our Fortran API does not change

  - it is easy for us to implement

  - it's easy for other people to write their own readers and writers

  - the Matlab reader might be a little faster, though the HDF 
    version can also be made pretty fast

But there are also some drawbacks:

  - a lot of space is wasted. The header (including attribute text)
    takes about 40K bytes, and each profile currently takes about
    50K, regardless of what's in it; so 1000 profiles will always take 
    up about 50Mbytes

  - any change whatsoever in the file format breaks compatibility
    with earlier file versions, so you can't casually add or move
    fields or change any array sizes.

  - the problem of managing attributes remains; these can be saved
    as a char array in the header, but then something as trivial as
    increasing the number of possible attributes or the max length
    of the attribute text breaks compatibility.

  - a convention (e.g., IEEE little-endian floats and int's) is 
    needed to deal with various machine word-order and data types

  - any archival value for the data is reduced significantly

It's because of such concerns that HDF was invented in the first
place.

The fixed record format would be somewhat similar to the setup with
the kcarta binary output format, but a little cleaner, as we have
an explicitly defined fixed set of fields, rather than fields
defined implicitly from loops and print statements, and no embedded
Fortran record lengths.

There are also various compromise designs, such as implementing
variable field sets without HDF, that have some advantages but do
not lead to simple readers for other users.  Since we already have
a working HDF implementation I would be inclined to stick with it,
unless it looks like other people writing their own readers will 
have problems.




