%
% matlab code for the read tests in rtptest3
%

% npro and nlevs must match paramters the file was created with
npro = 36000;
nlevs = 102;

% a precaution
hdfml('closeall')

[h, ha, p, pa] = rtpdump('rtptest3.hdf');

for k = 1 : npro

  for j = 1 : nlevs
    if (p.plevs(j,k) ~= (j-1) * 10 + 1 ) || (p.ptemp(j,k) ~= 200 + (j-1))
      fprintf(1, 'plevs or ptemp error prof %d lev %d\n', k, j);
      return
    end
  end

  for j = 1 : nlevs
    if (p.gas_1(j,k) ~=  100 + (j-1) + 1)
        fprintf(1, 'gas_1 error prof %d lev %d\n', k, j);
        return
      end
    end

  for i = 1 : h.nchan
    if p.robs1(i) ~= mod((i-1), 97)
      fprintf(1, 'robs1 error prof %d chan %d\n', k, i);
      return
    end
  end

  for i = 1 : h.nchan
    if p.rcalc(i) ~= mod((i-1), 101)
      fprintf(1, 'rcalc error prof %d chan %d\n', k, i);
      return
    end
  end
end

