% KVBESTUNIT read string value in ci-standard units
% KAZAN viewer routine
% 
% [str_value, str_koefficient, val_koefficient] = kvgetvalue(str)

% boep 14-mar-04, MPI for Bioinorganic Chemistry 


function [str, strunit, kf] = kvbestunit(val, unit)
prefix = {'p', 'n', 'u', 'm', '', 'k', 'M', 'G', 'T', 'Y'};
koeff  = [1E-12, 1E-9, 1E-6, 1E-3, 1, 1E3, 1E6, 1E9, 1E12, 1E345];
val = abs(val);
for kk=2:length(prefix)
    if val < koeff(kk)
        strunit = prefix{kk-1}; 
        kf = koeff(kk-1);
        str = [num2str(val/kf), ' ', strunit, unit];
        break;
    end
end
