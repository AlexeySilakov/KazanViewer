function varargout = SPAread(filename, formShow)
PosTitle = 30;
PosNPoints = 564; %hex2dec('234')
PosWNunmbers = 576;
PosDataHeader = 288;

AppWin = 620;%? not really

SomeNumber1 = 632;
PosLaserFrequency = 640;

PosInfo = 896;
PosPosData = 386;
LenBufText = 1728-40; %padding of text between processed data and the interferogram.
fid = fopen(filename, 'r');
if fid<0, error(sprintf('File "%s" not found.', filename)); end
fseek(fid, 0, 'eof');
FSize = ftell(fid);

fseek(fid,PosTitle,'bof');
Title = char(nonzeros(fread(fid,255,'uint8'))');

fseek(fid,PosNPoints,'bof');
nPoints = fread(fid,1,'int32');
fseek(fid,PosWNunmbers,'bof');

maxX = fread(fid,1,'single');
minX = fread(fid,1,'single');
un1 = fread(fid,1,'single');
NumberOfScanPoints = fread(fid,1,'int32');
ResolutionPoints = fread(fid,1,'int32');
NumberOfSampleScans = fread(fid,1,'int32');
un5 = fread(fid,1,'single');
NumberOfFFTpoints = fread(fid,1,'int32');

fseek(fid,PosLaserFrequency,'bof'); 
LaserFrequency = fread(fid,1,'single');

% The starting byte location of the absorbance data is stored in the
% header. It immediately follows a flag value of 3

fseek(fid,PosPosData,'bof'); 
PosData = fread(fid,1,'int32');

% Flag=0;
% fseek(fid,PosDataHeader,'bof');
% IDS = fread(fid,PosNPoints-PosDataHeader,'uint16')
% id3 = find(IDS==3);
% td = IDS(id3+1); PosData = td(td~=0);
% id4 = find(IDS==4);
% td = IDS(id4+1); PosInfo = td(td~=0);
% 
% while Flag ~= 3
%     Flag = fread(fid,1,'uint16');
% end
% PosData=fread(fid,1,'uint16')';
% 
% Flag=0;
% fseek(fid,PosDataHeader,'bof');
% while Flag ~= 4
%     Flag = fread(fid,1,'uint16');
% end
% PosInfo = fread(fid,1,'uint16')';

fseek(fid, PosInfo, 'bof');
Info = char(fread(fid,PosData-PosInfo, 'uint8')');


fseek(fid,PosData,'bof');
Data = fread(fid,nPoints,'single'); % Absorbance Data

SomeInfo = char(fread(fid,LenBufText, 'uint8')');

tD = fread(fid,FSize-ftell(fid),'single'); % Other datasets, e.g. interference data


fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
Li = [1, strfind(Info, newline)];

idx = strfind(Info, 'Final format:');
if ~isempty(idx)
    eidx = find(Li>idx(end));
    modes = {'Interferogram', 'Single Beam', 'Transmittance', 'Absorbance', 'Kubelka-Munk', 'Photoacoustic', 'Reflectance', 'log(1/R'};
    
    StrMode = strtrim(Info(idx(end):Li(eidx(1))));
    for ii = 1:numel(modes)
        if contains(StrMode, modes{ii})
            NumMode = ii;
            break;
        end
    end
else
    NumMode = 2;
    StrMode = 'Single Beam';
end
dsc = [];
dsc.Mode = StrMode;
for ii = 1:(numel(Li)-1)
    dsc = setfield(dsc, ['c', num2str(ii)], strtrim(Info(Li(ii):Li(ii+1))));
end
Li = [1, strfind(SomeInfo, newline)];
for ii = 1:(numel(Li)-1)
    dsc = setfield(dsc, ['s', num2str(ii)], strtrim(SomeInfo(Li(ii):Li(ii+1))));
end

if formShow == 1
    y = Data;
    x = flipud(linspace(minX, maxX,nPoints).');
    xy = 0;
    xlab = 'Wavenumbers, cm-1';
else
    if NumMode~=2, 
        nP = floor(numel(tD)/2)-3;
        InterfData = tD([1:nP]+nP+3, 1); %tD(1:nP, 1); %, tD([1:nP]+nP+3)
    else
        nP = floor(numel(tD))-3;
        InterfData = tD(1:nP, 1);
    end

    df = (maxX-minX)/nPoints;
    
    y = InterfData;
    x = [1:size(y, 1)]'/df;
    xy = ones(size(y, 2), 1);
    xlab = 'Mirror position, cm';
    
    
end

ax.x = x;
ax.title = Title;
ax.xlabel = xlab;
ax.ylabel = '';
dsc.BackwardPeakLocation = '0 '; % this is not recorded in OMNIC but needed for FTIRBox
dsc.HighFoldingLimit = sprintf('%f cm-1', LaserFrequency); % this is not recorded in OMNIC but needed for FTIRBox

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
