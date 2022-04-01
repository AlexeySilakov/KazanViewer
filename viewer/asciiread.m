% ASCIIREAD Read data from ascii files
%
%   data = asciiread(filename, delimiter, fstcol)
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


function varargout = asciiread(filename, delimiter, fstcol);
fid = fopen(filename, 'r');
if fid == -1, error('sorry, file not found'); return; end

k = 1; a='';
strings = {};
while ~feof(fid)
  try
    a = fgets(fid);
    F = sscanf(a, ['%f' delimiter]); %%%%% alsi 27.03.2008
    
    if ischar(F)
        strings{end + 1, 1} = a;
    else
        In(k,:) = F';
        k = k+1;
    end
  catch
    strings{end + 1, 1} = a;
  end
end
fclose(fid);

if fstcol, 
  if size(In, 2)>1
    ax.x = In(:, 1);
    Y = In(:, 2:end);
  else
    ax.x = (1:size(In, 1))';
    Y = In;
  end
else
  ax.x = (1:size(In, 1))';
  Y = In;
end

ax.type = 'data';
ax.xlabel = '?';
ax.title = '';
dsc.ASCII='ASCII loader for KAZAN viewer';
dsc.Delimiter=delimiter;
if fstcol,dsc.FirstColumn='X';else dsc.FirstColumn='Y';end
dsc.FirstColumn = ['Treated as ', dsc.FirstColumn];
dsc.size_x = num2str(size(Y,1));
dsc.size_y = num2str(size(Y,2));
if ~isempty(strings)
    for n = 1:size(strings, 1)
      dsc = setfield(dsc, ['c', num2str(n)], strings{n});
    end
end

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