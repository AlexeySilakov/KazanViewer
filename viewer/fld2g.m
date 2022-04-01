% function res = fld2g(fld, freq)
% calculates g-factor at specified field position 
% fld on spectrometer with frequency freq

function res = fld2g(fld, freq)

if(nargin~=2), error('Usage: res = fld2g(fld, freq).'); end
res = freq * planck ./ fld / bmagn;