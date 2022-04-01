% BRUKERREAD Read data from .dsc/.dta and .par/.spc data files
%
%   data = brukerreadJSS(filename)
%   [ax,data] = brukerreadJSS(filename);
%   [ax,data,dsc] = brukerreadJSS(filename);
%
%   Reads data and axis information from Bruker
%   files DSC/DTA and PAR/SPC which contains 'JSS' parameter. 
%   ax contains fields x,y,z 
%   depends on data dimension and may contain
%   other fields  as well
%   dsc is a cell array of description strings

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003-2004
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

% bme,   4-dec-03, MPI
% allsi, 20-jan-04, MPI

function varargout = brukerreadJSS(filename);

[path,name,ext] = fileparts(filename);

switch upper(ext)
case {'.DSC', '.DTA'}
  fname = fullfile(path,[name,'.dta']);
  dscname = fullfile(path,[name,'.dsc']);
case {'.PAR', '.SPC'}
  fname = fullfile(path,[name,'.spc']);
  dscname = fullfile(path,[name,'.par']);
end

  dsc = XEPRdsc(dscname);
  ax = XEPRparJSS(dsc);
  if ~isfield(dsc, 'DOS')
      dsc.DOS = '0';
  end
try 
%   % from EasySpin by Stefan Stoll www.esr.ethz.ch 
%   [y] = eprload(upper(filename));
% catch
  hdl = fopen(fname);
  [fname, mode, mformat] = fopen(hdl);
  
  
  if strcmp(dsc.DOS, 'Format')
      y = fread(hdl, 'float32', 'ieee-le');
  else
      y = fread(hdl, 'long', 'ieee-be');
  end
  fclose(hdl);              
end

if nargout > 1

  if isempty(ax.x) | size(ax.x, 1)~=size(y, 1)
    ax.x = [1:size(y, 1)].';    
  end
  if isempty(ax.y) | size(ax.y, 1)~=size(y, 2)
    ax.y = [1:size(y, 2)].';    
  end

end

% assign output depending on number of output arguments
% requested
switch nargout
 case 1,
   varargout = {y};
 case 2,
   varargout = {ax, y};
 case 3,
   varargout = {ax, y, dsc};
end
