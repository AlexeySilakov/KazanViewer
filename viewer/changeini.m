function changeini(varargin)
if nargin >0
    fname = varargin{1};
    structname = varargin{2};
    varname = varargin{3};
    if ischar(varargin{4})
        varvalue = varargin{4};
    else
        error('varvalue is not a char string');
        return
    end
else
    return;
end
fid = fopen([fname], 'r');
if fid == -1, return; end
strings = {};
count = 0;
this_is_it = 0;
job_done = 0;
while ~feof(fid)
    a = fgets(fid);
    if ~ischar(a), break; end
    tstr = trim(sscanf(a, '%s'));
    strings{end+1} = tstr;
    count = count + 1;
    if isempty(tstr)
        % do nothing
    elseif strcmp(tstr(1), '[')
        currstruct = strtok(tstr(2:end), ']');
        if strcmp(currstruct, structname), 
            this_is_it = 1;
        else 
            if this_is_it == 1&~job_done 
                % end of the chapter, but desired string was not found
                strings{end} = [varname, ' = ', varvalue]; 
                strings{end+1} = tstr;                
            end
        end
    elseif strcmp(tstr(1), '%')|strcmp(tstr(1), ';')
        % do nothing
    else
        [cvarname, cvarval] = strtok(tstr, '=');
        cvarname = trim(cvarname);
        if strcmp(cvarname, varname)&this_is_it
            strings{end} = [varname, ' = ', varvalue];
            job_done = 1;
        end
    end

end
fclose(fid);
if this_is_it == 0
    strings{end+1} = ['[',structname, ']'];
    strings{end+1} = [varname, ' = ', varvalue];
end

fid = fopen([fname], 'w');
if fid == -1, return; end
for cc = 1:length(strings)
    fprintf(fid, '%s\r\n', strings{cc});
end
fclose(fid);