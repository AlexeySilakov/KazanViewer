% D00READ Read data from .d00/.exp data files
%
%   data = d00read(filename)
%   [ax,data] = d00read(filename);
%   [ax,data,dsc] = d00read(filename);
%
%   Reads data and axis information from .d00
%   and .exp spectrometer files. File extension 
%   is insignificant. ax contains fields x,y,z 
%   depends on data dimension and may contain
%   other fields  as well
%   dsc is a cell array of description strings

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

% stst, 14-mar-01, WIS
% bme,   4-dec-03, MPI  

function varargout = d00read(filename);

% .d00 reading
%----------------------------------------------
[path,name,ext] = fileparts(filename);

[h,msg] = fopen(fullfile(path, [name '.d00']));
error(msg);

% read in first three 16bit integers
Dims = fread(h,3,'int16').';
nDims = sum(Dims>1);

% read in data, complex
data = fread(h,[2,inf],'double');
data = data(1,:) + i*data(2,:);

% and shape into correct array size
data = reshape(data,Dims);

% close data file
fclose(h);

%----------------------------------------------
str = [];
if nargout>1
  daxis = {[],[],[]};
  dscaxis = {'','',''};
  GroupTagList = {'[x]','[y]','[z]'};
  [h,msg] = fopen(fullfile(path, [name '.exp']));
  error(msg);
  k = 1;
  while ~feof(h)
    l = fgetl(h);
    str = setfield(str, ['f' num2str(k)],l); k = k + 1;
    dimfound = find(strcmp(l,GroupTagList));
    if ~isempty(dimfound)
      txt = fgetl(h);
      str = setfield(str, ['f' num2str(k)],txt); k = k + 1;
      [dscaxis{dimfound}, txt] = strtok(txt, '=');
      t = taketokens(txt(2:end));
      dscaxis{dimfound} = [dscaxis{dimfound},', ',t{end}];
      if strcmp(t{3},'to'),
        start = sscanf(t{1},'%f');
        thend = sscanf(t{4},'%f');
        daxis{dimfound} = linspace(start,thend,Dims(dimfound));
      elseif strcmp(t{3},'step')
        start = sscanf(t{1},'%f');
        step = sscanf(t{4},'%f');      
        daxis{dimfound} = start + (0:Dims(dimfound)-1)*step;
      elseif strcmp(t{3},'around')
        center = sscanf(t{4},'%f');
        around = sscanf(t{1},'%f');      
        daxis{dimfound} = linspace(center-around/2,center+around/2,Dims(dimfound));
      end
    end
  end
  fclose(h);
  % filling unknown axis
  for i=1:3
    if Dims(i) > 1 & isempty(daxis{i})
      daxis{i}=1:Dims(i);
    end
  end
end

% TODO: 2and3D
if ~isempty(daxis), ax.x = daxis{1}; end
ax.xlabel = dscaxis{1};
ax.type = 'data';

% assign output depending on number of output arguments
% requested
switch nargout
case 1,
  varargout = {data};
case 2,
  varargout = {ax, data};
case 3,
  varargout = {ax, data, str};
end

return


function t = taketokens(line);
% returns a list of white-space delimited tokens contained in the line
r = deblank(line); % remove white-space
k = 1;
while ~isempty(r)
  [t{k},r] = strtok(r);
  k = k+1;
end
return
