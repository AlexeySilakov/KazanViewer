function res = nfreq2fld(nfreq, nuc)

if nargin < 2, nuc = '1H'; end

[I, gn] = nucdata(nuc);
res = nfreq * planck / gn / nmagn;