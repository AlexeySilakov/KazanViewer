% XWINNMR Read data from XWINNMR data files (1D)
%
%   data = XWINNMR(filename)
%   [ax,data] = XWINNMR(filename)
%   [ax,data,dsc] = XWINNMR(filename)
%
%   Reads data and axis information from fid./acqus
%   and 1r./1i./proc spectrometer files. File extension 
%   is insignificant. ax contains fields x,y,z 
%   depends on data dimension and may contain
%   other fields as well
%   dsc is a cell array of description strings

% bme,   4-dec-03, MPI  
% bme,   12-mar-04, MPI  

function varargout = WINNMR(filename);

[fpath,name,ext] = fileparts(filename);

SizeTD2 = 1;

% if strcmp(name,'fid')
dsc = XWINNMRpar(filename);
AQ_mod = str2num(dsc.AQ_mod);
seq = bitand(AQ_mod, 2) > 0;

SizeTD1 = str2num(dsc.TD)/(1+~seq);
SW=str2num(dsc.SW_h); % s
offset = str2num(safeget(dsc, 'O1', '0'));
ax.x = (SW*[-SizeTD1/2:SizeTD1/2-1]/SizeTD1).' + offset;
ax.xlabel = 'Frequency, Hz';
ax.reference = str2num(safeget(dsc, 'SFO1', '0'))*1E6;
ax.freq1 = str2num(safeget(dsc, 'BF1', '0'))*1E6;
ax.dx = 0;
if str2num(dsc.BYTORDA) > 0.5, byteorder = 'b';
else, byteorder = 'l';
end  
sz = SizeTD2*SizeTD1;
id = fopen(fullfile(fpath, [name,'.1r']), 'r', byteorder);
[p1,count] = fread(id, sz, 'float');
fclose(id);
id = fopen(fullfile(fpath, [name,'.1i']), 'r', byteorder);
[p2,count] = fread(id, sz, 'float');
fclose(id);
data = p1+j*p2;
%   
% else
%   dsc = XWINNMRpar(filename); 
%   SizeTD1 = str2num(dsc.SI);
%   
%   SW=str2num(dsc.SW_p); % Hz
%   SF=str2num(dsc.SF);   % spec freq MHz
%   off = str2num(dsc.OFFSET)*SF;
% 
%   ax.x = off-SW*[0:SizeTD1-1].'/(SizeTD1-1);
%   ax.xlabel = 'Frequncy, Hz';
% 
%   sz = SizeTD2*SizeTD1;
%   if str2num(dsc.BYTORDP) > 0.5, byteorder = 'b';
%   else, byteorder = 'l';
%   end
%   id = fopen(fullfile(fpath, '1r'), 'r', byteorder);
%   [Return1,count] = fread(id, sz, 'int32');
%   fclose(id);
%   
%   id = fopen(fullfile(fpath, '1i'), 'r', byteorder);
%   [Return2,count] = fread(id, sz, 'int32');
%   fclose(id);
%   
%   data = Return1 + j*Return2;
% end
ax.type = 'data';

% assign output depending on number of output arguments
% requested
switch nargout
  case 1,
    varargout = {data};
  case 2,
    varargout = {ax, data};
  case 3,
    varargout = {ax, data, dsc};
end