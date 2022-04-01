% cmbl Read data from ascii files
%
%   data = kv_cmblread(filename)
%   [ax,data] = asciiread(filename, ...);
%   [ax,data,dsc] = asciiread(filename, ...);
%
%   Reads data and axis information from dat
%   files. Aditional arguments are always 
%   coming in pairs field, value. 
%   ax contains fields x,y,z depends on data 
%   dimension and may contain other fields 
%   as well. dsc is a cell array of description 
%   strings. 

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

% Silakov Alexey, ?-nov-03, MPI  
% bme,   4-dec-03, MPI  


function varargout = kv_cmblread(filename)
fid = fopen(filename, 'r');
if fid == -1, error('sorry, file not found'); return; end

c = fread(fid, 'char');
fclose(fid);

str = char(c');

idx1 = strfind(str, '<ColumnCells>')+numel('<ColumnCells>')+1;
idx2 = strfind(str, '</ColumnCells>')-1;
ax.x = str2num(str(idx1(1):idx2(1)));
Y = str2num(str(idx1(2):idx2(2)));

dsc.size_x = num2str(size(Y,1));
dsc.size_y = num2str(size(Y,2));

relFields = {'FileIDString', 'CreationDateTime', 'StartTime', 'ModifiedDateTime', 'DataObjectName'};
for ii = 1:numel(relFields)
    sss = ['<', relFields{ii}, '>'];
    nnn = ['</', relFields{ii}, '>'];
    idx1 = strfind(str, sss)+numel(sss);
    idx2 = strfind(str, nnn)-1;
    for ss = 1:min(numel(idx1), numel(idx2))
        dsc = setfield(dsc, relFields{ii}, str(idx1(ss):idx2(ss)));
    end
end

ax.type = 'data';
ax.xlabel = 'time, s';
ax.xlabel = 'time, s';
ax.title = '';
dsc.ASCII='cmbl loader for KAZAN viewer';


% assign output depending on number of output arguments
% requested
switch nargout
  case 1,
    varargout = {Y};
  case 2,
    varargout = {ax, Y};
  case 3,
    varargout = {ax, Y, dsc};
end