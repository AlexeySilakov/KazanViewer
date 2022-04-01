% KVGETVALUE read string value in ci-standard units
% KAZAN viewer routine
% 
% [val, str_unit, str_koefficient] = kvgetvalue(str)

% boep 14-mar-04, MPI for Bioinorganic Chemistry 


function [val, unit, pref, pref_val] = kvgetvalue(str, varargin)
prefix = ['p','n', 'u', 'm', 'k', 'M', 'G', 'T'];
koeff  = [1E-12, 1E-9, 1E-6, 1E-3, 1E3, 1E6, 1E9, 1E12];
idx = (str >= '0' & str <='9') | str == '.' | ...
    upper(str) == 'E' | str == '+' | str == '-';
pref = '';
pref_val = 1;
val = str2num(str(idx));
unit = str(~idx);
unit = unit(unit~=' ');
if nargin>1
    if varargin{1}==0
        pref = '';
        pref_val = 1;
        return
    end
end
if length(unit) > 1
    if ~isempty(unit)
        kk = findstr(prefix, unit(1));
        if ~isempty(kk)
            val = val * koeff(kk);
            unit = unit(2:end);
            pref = prefix(kk);
            pref_val = koeff(kk);
        end
    end
end
