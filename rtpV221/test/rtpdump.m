
function [head, hattr, prof, pattr] = rtpdump(hfile)

% simplified version of rtpread for tests only

addpath /home/motteler/shome/mot2008/hdf/h4tools

fprintf(1, 'RTP TEST\n')

% open the HDF file
file_id = hdfh('open', hfile, 'read', 0);
if file_id == -1
  error('HDF hopen failed');
end

% initialize the V interface
status = hdfv('start',file_id);
if status == -1
  error('HDF vstart failed');
end

% read the header and profiles
[head, hattr] = vsfid2mat(file_id, 'header');
[prof, pattr] = vsfid2mat(file_id, 'profiles');

% end vgroup interface access
status = hdfv('end',file_id);
if status == -1
  error('HDF vend failed')
end

% close the HDF file
status = hdfh('close',file_id);
if status == -1
  error('HDF hclose failed')
end

