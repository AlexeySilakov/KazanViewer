% G2FLD converts g-factors to field
%
% res = g2fld(g, mw_freq)

function res = g2fld(g, freq)
res = freq * 6.6261e-034./ g / 9.2740e-024;