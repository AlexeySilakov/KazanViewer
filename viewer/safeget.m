% function res = safeget(strct, fld, deflt)
% Returns the field of the structure or default 
% value if field is absent.
% if default value is array then return value
% has not less elements than in this array 
% (not in char case)

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retain all rights. 
% Contact: epel@mpi-muelheim.mpg.de

function res = safeget(strct, fld, deflt)
if isfield(strct, fld)
    res = getfield(strct, fld);
    sd = size(deflt, 2);
    sr = size(res, 2);
    if sr < sd
        if iscell(deflt)
            [res{sr+1:sd}] = deal(deflt{sr+1:sd});
        elseif ~ischar(deflt)
            res(sr+1:sd) = deflt(sr+1:sd);
        end
        
    end
else
    res = deflt;
end