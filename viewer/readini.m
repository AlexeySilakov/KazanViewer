function outstruct = readini(varargin)
% outstruct = readini
% outstruct = readini(desired_dir)
%   Reads 'ini' files.
%   Don't like spaces in a parameter name.
%   Returns struct 'outstruct' which has substructs calls by names of
%   the'chapters' in 'ini'files. They are have their own substructs calls
%   by names of parameters inside chapter.
%   All values are 'string'.
% AlSi 14.12.2004 for Kazan Viewer
outstruct = [];
if nargin >0
    dird = varargin{1};
else
    dird = '';
end
[fname, pname]= uigetfile([dird, '*.ini']);
if ~ischar(fname) return; end
fid = fopen([pname, fname], 'r');
currstruct = 'noname';
outstruct.noname = [];
while ~feof(fid)
    a = fgets(fid);
    tstr = trim(sscanf(a, '%s'));
    if isempty(tstr)
        % do nothing
    elseif strcmp(tstr(1), '[')
        currstruct = strtok(tstr(2:end), ']');
    elseif strcmp(tstr(1), '%')|strcmp(tstr(1), ';')
        % do nothing
    else
        [varname, varval] = strtok(tstr, '=');
        varname = trim(varname);
        if ~isletter(varname(1)), flag = 'a'; 
        else    flag = ''; end
        eval(['outstruct.',currstruct, '.', flag, trim(varname), '= ''', trim(varval(2:end)), ''';']);
    end
end
fclose(fid);

