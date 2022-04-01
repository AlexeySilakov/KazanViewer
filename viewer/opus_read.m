function [ax,y,dsc]=opus_read(varargin)
% OPUS FT-IR data file loader
%   [ax,y,dsc]=opus_read(filename)
%   Loads one selected data section (see format description below )
% Made to be part of KAZAN Viewer
% alsi 24.11.2005
%
% FORMAT description
%   format is terrible. original data has no X-axis
%   Structure of the file is following:
%    first comes the difraction patern
%    second, FFT(baseline) cutted for some sertain region
%    and third FFT(data)-FFT(baseline)
%   ---- in the new version the order is reverse --------
%   If during experiment baseline is recordered (extention is usually *.0),
%   the third part is not present
%   Data is in the "float32" binary format. 
% HEADERS
%   Data sections are separated by headers, starting usually with "CSF" or
%   "FXV" and finish with END
%   There is also some common header with description of
%   experimental setup Parameter, which I was able to find by "High/Low
%   folding limit". Headers contain variables names in ascii code followed 
%   by values in binary coded (all sorts)
%   "float64" binary for HFL/LFL - upper and lower frequency for X-axis, in
%   FFT(data) :-) 

filename = varargin{1};
if nargin>1, datapart = varargin{2}; else datapart = 2; end
warning off MATLAB:nonIntegerTruncatedInConversionToChar;
fid=fopen(filename,'r', 'ieee-le');

%ntotal = 83236;
fseek(fid, 0, 'eof');
ntotal = ftell(fid);

frewind(fid);
tmpstr=char(fread(fid,ntotal,'schar')).'; % whole file in one string
frewind(fid);
tmpdat=fread(fid,ntotal/4,'float32'); % whole file in one array
frewind(fid);
tmpint=fread(fid,ntotal/4,'int32'); % whole file in one array

hfl = strfind(tmpstr, 'HFL');
fseek(fid, hfl+7, 'bof');
HighFold = fread(fid, 1, 'float64');
lfl = strfind(tmpstr, 'LFL');
fseek(fid, lfl+7, 'bof');
LowFold = fread(fid, 1, 'float64');

fxv = strfind(tmpstr, 'FXV'); % first X point
lxv = strfind(tmpstr, 'LXV'); % last X point
npts = strfind(tmpstr, 'NPT'); % Number of data points 
ndatas = min([length(fxv), length(lxv)]);
for ii = 1:ndatas
        fseek(fid, fxv(ii)+7, 'bof');    
    FreqStPoint(ii) = fread(fid, 1, 'float64');
        fseek(fid, lxv(ii)+7, 'bof');
    FreqEnPoint(ii) = fread(fid, 1, 'float64');
        fseek(fid, npts(ii)+7, 'bof');
    NDataPoints(ii) = fread(fid, 1, 'int32');
end

npts = strfind(tmpstr, 'RSN'); % Running sample number (size(data)?)
if ~isempty(npts)
    fseek(fid, npts+7, 'bof');
    NPoints = fread(fid, 1, 'int32');
else
    NPoints = 0;
end

%
%APT..2.5 mm 
%BMS KBr 
%CHN  Front 
%DTC  ...MCT 
%LPF  1  
%SRC MIR-Source 
%VEL  7   
%SGN -1  
%RGN -1  
%........END
%APF B3  
%HFQ  @?@
%LFQ ±@ 
%PHR @
%PHZ ...ML  
%SPZ... NO  
%ZFF  2   
%........END     
%AQM  DN  
%COR  NO  
%DEL     
%DLY     
%HFW ±@
%LFW     y@
%NSS   è  
%PLF AB  
%RES @ (bin 000001)
%TDL  &   
%.....END     
%CNM default 
% 
%EXP cryostat_mct4000-0_c.XPM ...Pun
%%% sample form
%SFM Cryostat 40K /2,5 mm  ..PM
%SNM Creinhardtii FeFeHase, HoxCO, rec, 200uM, 15uL ..K
%...END     
%CSF      ð?
%MXY    âÏ?
%MNY `vû¿
%DAT 17/10/2007 
%TIM 03:09:36 ech
%DXU WN  
%.....END
npts = strfind(tmpstr, 'ARS'); % number of bg. scans
if ~isempty(npts) 
    fseek(fid, npts(1)+7, 'bof');
    NBGScans = fread(fid, 1, 'int32');
else
    NBGScans = 0;
end    
npts = strfind(tmpstr, 'ASS'); % number of sample scans
if ~isempty(npts) 
    fseek(fid, npts+7, 'bof'); 
    NScans = fread(fid, 1, 'int32');
else
    NScans = 0;
end
npts = strfind(tmpstr, 'PKL'); % Peak Location
if ~isempty(npts) 
    fseek(fid, npts+7, 'bof'); 
    PeakLocation = fread(fid, 1, 'int32');
else
    PeakLocation = 0;
end
npts = strfind(tmpstr, 'PRL'); % Backward Peak Location
if ~isempty(npts) 
    fseek(fid, npts+7, 'bof'); 
    BackwardPeakLocation = fread(fid, 1, 'int32');
else
    BackwardPeakLocation = 0;
end

% res = strfind(tmpstr, 'RES'); % Running sample number (size(data)?)
% fseek(fid, res+3, 'bof');
% Resol = fread(fid, 1, 'int32');

%aqm = strfind(tmpstr, 'AQM');

dd = strfind(tmpstr, 'DAT');
tt = strfind(tmpstr, 'TIM');
dsc.Date = [tmpstr((dd(1)+8):(dd(1)+17))];
dsc.Time = [tmpstr((tt(1)+8):(tt(1)+15))];

nptse = strfind(tmpstr, 'EXP'); 
nptsf = strfind(tmpstr, 'SFM'); 
nptsn = strfind(tmpstr, 'SNM'); 
nptsen  = strfind(tmpstr, 'END'); 
eie= find(nptsen>nptsn(1));
nptsen = nptsen(eie(1));
% nptsen = 
if ~isempty(nptsn)
    SampleName = [tmpstr((nptsn+8):(nptsen(1)-2))];
%     ee = isletter(SampleName);
%     iee = find(ee);
%     SampleName = SampleName(iee(1):iee(end));
else
    SampleName = '';
end
dsc.SampleName = SampleName;
if ~isempty(nptse)
    dsc.Experiment = [tmpstr((nptse+8):(nptsf-4))];
%     ee = isletter(dsc.Experiment);
%     iee = find(ee);
%     dsc.Experiment = dsc.Experiment(iee(1):iee(end));
else
    dsc.Experiment = '';
end
if ~isempty(nptsf)
    dsc.SampleForm = [tmpstr((nptsf+8):(nptsn))];
%     ee = isletter(dsc.SampleForm);
%     
%     iee = find(ee);
%     if ~isempty(iee)
%         dsc.SampleForm = dsc.SampleForm(iee(1):iee(end));
%     end
else
    dsc.Experiment = '';
end

fclose(fid);

%stcom = strfind(tmpstr, 'CSF');
encom = strfind(tmpstr, 'END') + 3;

%%%%%%% alsi - new opus software hack 26.01.2012
if isempty(encom(encom<fxv(1)))
    encom = [hex2dec('1F0')/4 + 2, encom];
%     encom = [hex2dec('1F0')/4 + 2, encom];
end
inx = find(diff(fxv)>1500);
idxx = [1, inx+1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:length(fxv(idxx))
    tencom = encom(encom<fxv(idxx(ii)));
    tencom = sort(tencom);
    stcom(ii) = tencom(end);
end
ax.title = SampleName;
% 1-Interferogram, 2, 3 - FFT
if length(fxv)>2 % full set of datas
    datapart = min(length(fxv), datapart);
    if datapart==1,
        ax.xlabel = 'Mirror steps, pnts';
        ax.title = [ax.title, ' <S_IFG>'];
    elseif datapart==2,
        ax.xlabel = 'Wavenumber, cm^{-1}';
        ax.title = [ax.title, ' <S_SC>'];
    elseif datapart==3,
        ax.xlabel = 'Wavenumber, cm^{-1}';
        ax.title = [ax.title, ' <_AB>'];        
    end
elseif length(fxv)==2 % baseline
    % somehow in file with baseline first comming fft(baseline) and only
    % then interferogram (sado_maso inc. :)
    datapart = min(length(fxv), datapart);
    if datapart==2,
        ax.xlabel = 'Wavenumber, cm^{-1}';
        ax.title = [ax.title, ' <S_SC>'];
       datapart = 1;
    else
        ax.xlabel = 'Mirror steps, pnts';
        ax.title = [ax.title, ' <S_IFG>'];
        datapart = 2;        
    end
elseif length(fxv)==1
    % only the AB
    datapart = min(length(fxv), datapart);
        ax.xlabel = 'Wavenumber, cm^{-1}';
        ax.title = [ax.title, ' <S_SC>'];
end


y = tmpdat(floor((stcom(datapart)+3)/4)+1+(1 : NDataPoints(datapart)));

ax.x = linspace(FreqStPoint(datapart), FreqEnPoint(datapart), NDataPoints(datapart)).';
ax.y = 1;

ax.HFL = HighFold;
ax.LFL = LowFold;

dsc.HighFoldingLimit = [num2str(HighFold), ' cm-1'];
dsc.LowFoldingLimit = [num2str(LowFold), ' cm-1'];
dsc.NumberOfBGScans = [num2str(NBGScans)];
dsc.NumberOfSampleScans = [num2str(NScans)];
dsc.NumberOfDataSets = [num2str(length(fxv))];
dsc.PeakLocation = num2str(PeakLocation);
dsc.BackwardPeakLocation = num2str(BackwardPeakLocation);

dsc.c1 = '???';