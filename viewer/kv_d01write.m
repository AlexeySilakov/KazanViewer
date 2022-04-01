% function varargout=kv_d01write(filename, src, varargin)
% arguments 
%    format    double/float 
function varargout=kv_d01write(filename, src, varargin)

opt = [];
if nargin > 2,
    if ~mod(nargin-2,2)
        for kk=1:2:nargin-3
            opt=setfield(src.ax, lower(varargin{kk}), varargin{kk+1});    
        end
    else, error('Wrong amount of arguments')
    end
end

[path,name,ext] = fileparts(filename);
fname = fullfile(path,[name,'.d01']);
dscname = fullfile(path,[name,'.exp']);

if isfield(src, 'fname'), 
    orfname = src.fname;
else
    orfname = 'unknown';
end
if ~isfield(src.ax, 'xlabel'), src.ax.xlabel = ''; end
if ~isfield(src.ax, 'ylabel'), src.ax.ylabel = ''; end

ind = findstr(src.ax.ylabel, '''');
src.ax.ylabel(ind) = '';
ind = findstr(src.ax.xlabel, '''');
src.ax.xlabel(ind) = '';

dim = size(src.y);
dim(end+1:4)=1; dim = dim(1:4);
total = prod(dim);
im = sum(sum(imag(src.y)));
tmpy = reshape(src.y, total,1);

ndim1 = (im~=0) + 1;

fid=fopen(char(fname),'w', 'ieee-le');
fwrite(fid, ndim1,'uint32');       % number of headers, re/im etc.
if strcmp(safeget(opt, 'format', 'float'), 'double'),
   fwrite(fid,0,'uint32');            % format:0-double,1-float
   sformat='double';
else
   fwrite(fid,1,'uint32');            % format:0-double,1-float
   sformat='float';
end

for k=1:ndim1
  fwrite(fid, sum(dim > 1),'int32');
  fwrite(fid,dim.',  'int32');
  fwrite(fid,total,'int32');
end
fwrite(fid,real(tmpy),sformat);
if im,
  fwrite(fid,imag(tmpy),sformat);
end

fclose(fid);

fid = fopen(char(dscname), 'w');
fprintf(fid, '[general]\r\n');
fprintf(fid, ['name = ', safeget(src.ax, 'title', ''),'\r\n']);
if isfield(src.ax,'freq1'), fprintf(fid, 'freq1 = %f GHz\r\n\r\n', src.ax.freq1*1E-9); end
fprintf(fid, '[text]\r\nGenerated by Kazan Viewer\r\n');
fprintf(fid, 'from ''%s''\r\n[sweep]\r\n\r\n', orfname);
% sweep
if im,
    fprintf(fid, 'transient = I,1,1,a,b\r\n');
else
    fprintf(fid, 'transient = I,1,1,a\r\n');
end

xlabel=trim(strtok(src.ax.xlabel, ','));
if strcmp(xlabel, '?'), xlabel = 'no'; end
if isempty(xlabel), xlabel = 'x'; end
ylabel=trim(strtok(src.ax.ylabel, ','));
if strcmp(ylabel, '?'), ylabel = 'no'; end
if isempty(ylabel), ylabel = 'y'; end
if strcmp(xlabel, ylabel)
    ylabel = [ylabel,'_2'];  
end

if dim(2)>1,
    fprintf(fid, 'sweep0 = X,%d,%d,%s\r\n', dim(1), 1, xlabel);
    fprintf(fid, 'sweep1 = Y,%d,%d,%s\r\n\r\n', dim(2), 1, ylabel);
else
    fprintf(fid, 'sweep0 = X,%d,%d,%s\r\n\r\n', dim(1), 1, xlabel);
end
fprintf(fid, '[params]\r\n');

[tokstr, reststr]=strtok(src.ax.xlabel, ',');
unix = trim(reststr(2:end));
[tokstr, reststr]=strtok(src.ax.ylabel, ',');
uniy = trim(reststr(2:end));

if dim(2)>1
    fprintf(fid, '%s = %d %s to %d %s; %s@none\r\n', xlabel,...
        src.ax.x(1, 1), unix, src.ax.x(end, 1), unix, xlabel);
    ll = length(src.ax.y);
    dx = (src.ax.y(end, 1) - src.ax.y(1, 1))/(ll-1);
    if max(abs(src.ax.y - src.ax.y(1, 1) - [0:ll-1]'*dx))>dx/2 % AlSi 27.12.2005 (old version was wrong ... + [0:ll-1]'*dx... )
       str = '';
       for ii=1:ll
           str = [str, num2str(src.ax.y(ii)), ' ', uniy]; 
           if ii~=ll, str = [str, ', ']; end;
       end
        fprintf(fid, '%s = %s; %s@none\r\n\r\n', ylabel, str, ylabel); 
    else
        fprintf(fid, '%s = %d %s to %d %s; %s@none\r\n\r\n', ylabel,...
            src.ax.y(1, 1), uniy, src.ax.y(end, 1), uniy, ylabel); 
    end
else
    fprintf(fid, '%s = %d %s to %d %s; %s@none\r\n\r\n', trim(strtok(src.ax.xlabel, ',')),...
        src.ax.x(1, 1), unix, src.ax.x(end, 1), unix, trim(strtok(src.ax.ylabel, ',')));
end
fprintf(fid, '[aquisition]\r\n');
fprintf(fid, 'a = \r\n');
if im,
   fprintf(fid, 'b = \r\n');
end

% non-systematic but...
if isfield(src, 'dsc')
    fprintf(fid, '[dsc]\r\n');
    flds = fieldnames(src.dsc);
    for ii=1:length(flds)
      fprintf(fid, '%s=%s\r\n', flds{ii}, getfield(src.dsc,flds{ii}));
    end
end
fclose(fid);