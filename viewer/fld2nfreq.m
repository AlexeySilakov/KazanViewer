function res = fld2freq(fld, nuc)

if(nargin==0), error('Usage: res = fld2freq(fld, nuc). nuc is e.g. ''1H''.'); end
if nargin < 2, nuc = '1H'; end

[I, gn] = nucdata(nuc);
res = fld*gn*nmagn/planck;