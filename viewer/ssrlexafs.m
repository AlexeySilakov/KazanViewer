function varargout = ssrlexafs(filename, varargin)
%*************
% SSRL format
%*************
% 1. 800 bytes = 40 bytes (title MUST start with 'SSRL -') '\n' '0'
%              + 40 bytes (recording date:time) '\n' '0'
%              + 40 bytes (scaler_file) '\n' '0'
%              + 40 bytes (region_file) '\n' '0'
%              + 80 bytes (mirror_info) '\n' '0'
%              + 40 bytes (mirror param) '\n' '0'
%              + 80 bytes x 6 (User comments) '\n' '0'
%
% 2. NCOL*4 bytes (offsets)
%    NCOL*4 bytes (weights)
%
% 3. NCOL*20 bytes (labels for each column)
%
% 4. 4 bytes integer (npts)
%
% 5. mysterious 12 bytes starting w/ char(8),char(0), char(1), char(0)
%    probably it's not important.
%
% 6. NCOL*NPTS*4 bytes (data)
%
%
% Tsu-Chien Weng, 2002-07-24
% Copyleft(c), by the Penner-Hahn research group

[fpath,name,ext] = fileparts(filename);

fid=fopen(filename,'r', 'ieee-le');
Pts = 0;
Col = 0;
tt = [];
dtt=[];
try
    for ii = 1:12
        ff = ftell(fid);
        text = fgetl(fid);
        tt{end+1} = text;
        if ~isempty(strfind(tt{end}, 'PTS'))
            [nums] = sscanf(tt{end}, '%*s%d%*s%d');
            Pts = nums(1);
            Col = nums(2);
        end
    end
    if ~Pts
       error('Header was not found');
    end
    
    fseek(fid, 798, 'bof');
    ff = ftell(fid);
    
    Offsets = fread(fid, Col, 'real*4');
    Weights = fread(fid, Col, 'real*4');
    
    for ii= 1:Col
   
        text = fgetl(fid);
        dtt{end+1} = text;
    end
%      cc = fread(fid, 2, 'char');

ff = ftell(fid);
cc = fread(fid, 1, 'char');
     Dig = fread(fid, 4, 'int16');
     ccdd = fread(fid, 6, 'char');
ff = ftell(fid);
     DAta = fread(fid, Col*Pts, 'float');
     
catch
    disp(lasterr);
end

fclose(fid);
dsc.FileName = filename;
dsc.nCols = num2str(Col);
dsc.nPts  = num2str(Pts);
for ii = 1:length(tt)
    idx = isstrprop(tt{ii}, 'alphanum')+isstrprop(tt{ii}, 'punct')+isstrprop(tt{ii}, 'wspace');
    icc = find(idx~=0);

    dsc = setfield(dsc, ['c', num2str(ii)], tt{ii}(icc));
end
idxI0=4;
idxI1=5;
idxI2=6;

idxeV=3;

idxrtc=1;
idxFF = [];
idxSCA1 = [];
idxSCA2 = [];
idxICR = [];
for ii = 1:length(dtt)
    idx = isstrprop(dtt{ii}, 'alphanum')+isstrprop(dtt{ii}, 'punct')+isstrprop(dtt{ii}, 'wspace');
    icc = find(idx~=0);
        
    dsc = setfield(dsc, ['data_', num2str(ii)], dtt{ii}(icc));
    if ~isempty(strfind(dtt{ii}, 'rtc'))
        idxrtc = ii;
    end
    if ~isempty(strfind(dtt{ii}, 'I0'))
        idxI0 = ii;
    end
    if ~isempty(strfind(dtt{ii}, 'I1'))
        idxI1 = ii;
    end
    if ~isempty(strfind(dtt{ii}, 'I2'))
        idxI2 = ii;
    end
    if ~isempty(strfind(dtt{ii}, 'LYTLE'))
        idxLYTLE = ii;
    end
    if ~isempty(strfind(dtt{ii}, 'Achieved Energy'))
        idxeV = ii;
    elseif ~isempty(strfind(dtt{ii}, 'FF'))
        idxFF(end+1) = ii;
    end
    if ~isempty(strfind(dtt{ii}, 'SCA1'))
        idxSCA1(end+1) = ii;
    end
    if ~isempty(strfind(dtt{ii}, 'SCA2'))
        idxSCA2(end+1) = ii;
    end
    if ~isempty(strfind(dtt{ii}, 'ICR'))
        idxICR(end+1) = ii;
    end
end
y = reshape(DAta, Col, Pts).';

xax = (y(:, idxeV))/4;
[xax, idxSort]=sort(xax);
y = y(idxSort, :); %.*Weights(:, ones(Pts, 1)).';

if nargin>1
    func = varargin{1};
    ns = strip_names(func);
    for ii = 1:length(ns)
        if exist(['idx', ns{ii}])
            eval([ns{ii}, '= y(:, idx', ns{ii},');']);
        else
            error(sprintf('input formula %s contains unknown channel name %s', func, ns{ii}));
        end
    end
    eval(['y= ', func, ';']);
    
end
% I0 = (y(:, idxI0*ones(Col, 1)));
% rtc = (y(:, idxrtc));
ax.x = xax; %[1:Pts].';
ax.y = [1:Col].';
ax.title = trim(dsc.c8);

% y = y(:, idxFF);
% y = (y./I0);
% y = y(:, idxSCA);
% y = (y./I0);
% y = y(:, idxSCA);
switch nargout
 case 1,
   varargout = {y};
 case 2,
   varargout = {ax, y};
 case 3,
   varargout = {ax, y, dsc};
end

function names = strip_names(str)
tokill = {'ones', 'size', 'log', 'sum', 'exp', '*', '/', '+', '-', '(', ')', '.', ':'};
tokill_len = [3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1]; % just to speed thing up a bit

fstr = [' ', str, ' '];
for ii=1:length(tokill)
    fstr = strrep(fstr, tokill{ii}, ' ');
end

idx = isstrprop(fstr, 'alphanum');
id1 = find(diff(idx)>0)+1;
id2 = find(diff(idx)<0);

names = {};
for ii = 1:length(id1)
    tr = isstrprop(fstr(id1(ii):id2(ii)), 'alpha');
    if sum(tr)
        names{end+1} = fstr(id1(ii):id2(ii));
    end
end

