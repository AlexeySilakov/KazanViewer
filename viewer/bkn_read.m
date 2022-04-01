function [ax,y,dsc]=bkn_read(varargin)
filename = varargin{1};
if nargin>1, datapart = varargin{2}; else datapart = 1; end
warning off MATLAB:nonIntegerTruncatedInConversionToChar;

fid=fopen(filename,'r', 'ieee-le');
CC = fread(fid, 'char'); 
fclose(fid); 

SS = char(CC');
idx1 = strfind(SS, 'Sample');
idx2 = strfind(SS, 'BCB');
idx31 = strfind(SS, 'TContinuumStore');
idx32 = strfind(SS, 'TGraphStore');

count = 0;
cha = '%$';
for ii = 1:18
    tidx1 = strfind(SS, ['Sample ', num2str(ii), cha(1)]);
    tidx2 = strfind(SS, ['Sample ', num2str(ii), cha(2)]);
    if isempty(tidx1) &  isempty(tidx2) 
        break;
    end
    count = count+1;
end
nData = count; %round(numel(idx2)/3); 

if datapart>nData, 
    datapart = nData;
    warning(['dataset # > number of datasets in the batch file (', num2str(nData), ')']);
end
idx = find(idx1>idx2(datapart)); 
dat1 = (idx2(datapart))+3;
dat2 = round(idx1(idx(1)))-2;

ttidx = [idx31, idx32];
Sidx = find(ttidx>dat2);

Dsc = SS((dat2+11):ttidx(Sidx(1)));
cDsc = CC((dat2+11):ttidx(Sidx(1)));


fid = fopen(filename, 'rb', 'L');
fseek(fid, dat1, 'bof');
A = fread(fid, (dat2-dat1)/4, 'float32');
fclose(fid) ;

y(:, 1) = A(2:2:end);

ax.x(:, 1) = A(1:2:(end-1));
ax.y = 1;

dsc.c1 = filename;
dsc.nDataSets = num2str(nData);
ax.title = ['Dataset #', num2str(datapart), ' out of ', num2str(nData)];

% 
% disp(SS(dat2:ttidx(Sidx(1))));
lbrake = char(17);
stridx = [1, strfind(Dsc, '')];
for ii = 2:numel(stridx)
    dsc = setfield(dsc, ['d', num2str(ii+1)], Dsc(stridx(ii-1):stridx(ii)) );
end


% for ii = 1:numel(idx1)
%     dsc = setfield(dsc, ['c', num2str(ii+1)], SS([0:10]+idx1(ii)) );
% end
