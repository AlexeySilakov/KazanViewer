% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: silakov@mpi-muelheim.mpg.de

% Alexey Silakov, MPI-BAC 2005
% alsi, 13.06.05

function varargout=moss_webbin(fname, varargin)

[path,name,ext] = fileparts(fname);
switch lower(ext)
    case '.a'
        [y, ical, rcal, TXT_dsc] = readascii(fname);
        dsc.type = 'ascii WEB data';
    case '.b'
        [y, ical, rcal, TXT_dsc] = readbin(fname);
        dsc.type = 'binary WEB data';
    otherwise
        error('Extention ''', ext, ''' is not supported by "moss_webbin"');
end

ax.x = [(1:ical(1))-rcal(6)].'*rcal(7);
ax.y = 0;

ax.xlabel = 'Velosity, mm/s';
ax.type = 'data';
dsc.MagnField = [num2str(rcal(1)), ' kG'];
dsc.Temperature = [num2str(rcal(2)), ' K'];
dsc.Baseline = [num2str(rcal(3))];
dsc.some1 = [num2str(rcal(5))];
dsc.ZeroChanel = [num2str(rcal(6))];
dsc.slope = [num2str(rcal(7))];
dsc.some2 = [num2str(rcal(8))];
dsc.some3 = [num2str(rcal(9))];
dsc.some4 = [num2str(rcal(10))];
dsc.Description = TXT_dsc;


Gamma_pol = {'??? (-3)',...
        'powder (-2)',...
        'Perpendicular / split coil (-1)',...
        'Unpolarized (0)',...
        'Circular (1)',...
        'Linear (2)'};
dsc.Gamma_polarisation = Gamma_pol{ical(3)+4};
% IPOL=-3 : INTENSITY TENSOR AND INTENSITY VECTOR ARE TO BE
%           CALCULATED. NO INPUT FOR EX,EY,EZ NECESSARY
% IPOL=-2 : POWDER MEAN VALUE; NO INPUT FOR EX,EY,EZ NECESSARY
% IPOL=-1 : PERPENDICULAR CONFIGURATION (SPLIT COIL)
%           EX,EY,EZ IS TO BE PARALLEL TO THE EXTERNAL FIELD
% IPOL=0  : UNPOLARIZED GAMMA-RAY PARALLEL TO EX,EY,EZ
% IPOL=1  : CIRCULAR POLARISATION,WAVE VECTOR PARALLEL EX,EY,EZ
% IPOL=2  : LINEAR POLARISATION; IN 57FE (EX,EY,EZ) IS THE VECTOR
%           PERPENDICULAR TO THE WAVE VECTOR AND PERPENDICULAR TO
%           THE VECTOR POTENTIAL A


% some analysis (as in original WEB-Moss program)
if rcal(3)~=0
    dsc.Effect = num2str(-100*min(y)*rcal(4)/rcal(3));
else
    dsc.Effect = '0';
end
dsc.Area = num2str(rcal(4));%num2str(sum(y));


ax.title = trim(dsc.Description);

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


function [y, ical, rcal, TXT_dsc] = readbin(fname)
% 
fid=fopen(char(fname),'r');
if fid<0, error(['File ''', char(fname), ''' does not exist.']); end
TXT_dsc = char(fread(fid, 64, 'char'))';
precision = 'real*4';

fseek(fid, 432, 'bof');
[ical,count] = fread(fid, 10, 'integer*4');
[rcal,count] = fread(fid, 10, 'real*4');

[yval,count] = fread(fid, ical(1)*4, precision);
fclose(fid);

y = yval; %*sqrt(rcal(4))/abs(sum(yval));
% up to now no noise bars in the spectrum
% ynoi(i) = sqrt( rcal(3)+rcal(4)*yval(i) ) / abs (rcal(4))

function [y, ical, rcal, TXT_dsc] = readascii(fname)
fid=fopen(char(fname),'r');

a = fgets(fid);
TXT_dsc = a;
for i = 1:5
    a = fgets(fid);    
end

ical = [];
rcal = [];
for i = 1:2
    a = fgets(fid);
    ical = [ical; sscanf(a, '%i')];
end
for i = 1:2
    a = fgets(fid);
    rcal = [rcal; sscanf(a, '%f')];    
end
yt = [];
strings = '';
while ~feof(fid)
  try
    a = fgets(fid);
    F = sscanf(a, ['%f \t']);
    yt(end+1,:) = F';
  catch
    strings{end + 1, 1} = a;
  end
end
y = yt(1:ical(1), 2);
fclose(fid);