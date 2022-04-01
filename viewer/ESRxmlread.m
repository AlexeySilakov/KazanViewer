function varargout = ESRxmlread(filename, datapart)
fid = fopen(filename, 'r');
if fid == -1, error('sorry, file not found'); end


% fname = 'c:\My documents\PSU\Data\MS5000\test\20180503_185640345_Chromium_E1703.xml';
% fid = fopen(fname, 'r');
% cc = char(fread(fid, 'char')).';
ss = {};
par = [];
lin = [];
meas = [];
rec = [];
while ~feof(fid)
    ss{end+1} = fgetl(fid);
    if ~isempty(strfind(ss{end}, '<Curve'));
        lin(end+1) = numel(ss);
    end
    if ~isempty(strfind(ss{end}, '<Param Name'));
        par(end+1) = numel(ss);
    end
    if ~isempty(strfind(ss{end}, '<Measurement Name'));
        meas(end+1) = numel(ss);
    end
    if ~isempty(strfind(ss{end}, '<Recipe Name'));
        rec(end+1) = numel(ss);
    end
end
fclose(fid);

Pars = [];
for ii = 1:numel(par)
    
    str = trim(ss{par(ii)});
        tStr = readParams(str, []);
%     idx = strfind(str, '"');
    ParName = tStr.Name; %str((idx(1)+1):(idx(2)-1));
    if isfield(tStr, 'Unit')
        idss = isstrprop(tStr.Unit, 'alphanum');
        U = tStr.Unit(idss);
        idss = find(uint8(U)<123);
        ParName = [ParName, '_', U(idss)];
    end
    
    enst = strfind(str, '</Param>')-1;
    best = strfind(str, '>')+1;
    Val = str(best(1):enst);
    Pars = setfield(Pars, ParName, Val);
end
for ii = 1:numel(rec)
    str = trim(ss{par(ii)});
    Pars = readParams(str, Pars);
end
base64 = org.apache.commons.codec.binary.Base64;
Curves =  {};
Bvals = 0;
SinVal = 0;
CosVal = 0;
Proc  = 0;
for ii = 1:numel(lin)
    str = trim(ss{lin(ii)});
    enst = strfind(str, '</Curve>')-1;
    if isempty(enst)
        Dat{ii} = 0;
        continue;
    end
    best = strfind(str, '>')+1;
    Val = str(best(1):enst);
    
    
    y = base64.decode(uint8(Val));

    nPoints = floor(numel(y)/9);
    sa = reshape(y(1:(nPoints*9)), 9, nPoints);
    sa = sa(1:8, :);
    Data = typecast(sa(:), 'double');
%     Curves(ii).Data = Data;
    Stru.Data = Data;
    Stru = readParams(str, Stru);
%     in = strfind(str, '>');
%     str = str(1:in(1));
%   
%     iv = strfind(str, '"');
%     ie = strfind(str, '=');
%     is = strfind(str, ' ');
%     
%     Stru.Data = Data;
%     for jj = 1:numel(ie)
%         ic = find(is-ie(jj)<0);
%         Nam = str((is(ic(end))+1):(ie(jj)-1));
%         idx1 = iv(2*jj-1)+1;
%         idx2 = iv(2*jj)-1;
%         Val = str(idx1:idx2);
%         Stru = setfield(Stru, Nam, Val);
%     end
    
    Curves{ii} = Stru;
    switch Stru.Name
        case 'MWAbsorption'
            Proc = ii;
        case 'MWAbsorption sinus'
            SinVal = ii;
        case 'MWAbsorption cosinus'
            CosVal = ii;
        case 'BField'
            Bvals = ii;
    end
    
end

if ~isempty(meas)
    str = trim(ss{meas(1)});
    Pars = readParams(str, Pars);
%     in = strfind(str, '>');
%     str = str(1:in(1));
%   
%     iv = strfind(str, '"');
%     ie = strfind(str, '=');
%     is = strfind(str, ' ');
%     
%     for jj = 1:numel(ie)
%         ic = find(is-ie(jj)<0);
%         Nam = str((is(ic(end))+1):(ie(jj)-1));
%         idx1 = iv(2*jj-1)+1;
%         idx2 = iv(2*jj)-1;
%         Val = str(idx1:idx2);
%         Pars = setfield(Pars, Nam, Val);
%     end
    
end

ax.freq1 = str2num(Pars.MwFreq)*1e9;

if lower(datapart(1))=='p'
    Y = Curves{Proc}.Data;
    dt = str2double(Curves{Proc}.XSlope);
    T = 0:dt:(numel(Y)-1)*dt;
    ax.title = 'Processed Data';
elseif lower(datapart(1))=='s'
    % 1) shrink data to 4096 points with appropriate filtering 
    % 2) subtract linear slope
    noRawData = 1;
    if SinVal >0
       if ~isempty(Curves{SinVal}.Data)
           noRawData = 0;
       end 
    end
    if noRawData==0 
        ty = Curves{SinVal}.Data; % + 1i*Curves{CosVal}.Data;
        dt = str2double(Curves{SinVal}.XSlope);
        T = 0:dt:(numel(ty)-1)*dt;

        nPoints = 2^12;

        nP = floor(numel(ty)/nPoints);
        sP = floor(numel(ty)/nP); % closest number to nPoints

        T = mean(reshape(T(1:(nP*sP)), nP, sP), 1).' - dt*nP/2;
        ccy = reshape(ty(1:(nP*sP)), nP, sP);
        Y = mean(ccy, 1).';
        Y = kv_mvgavg(Y, numel(Y)/nPoints, 'savgol', 1);
%         Y = loc_rcfilt(Y, str2num(Pars.SweepTime_s), convtime);
        %%%%%%% quick extraction of fisrt order polinomial for baseline
        idx = [1:sP]; 
        Sl = mean(diff(Y));
        Y = Y-[1:sP]'*Sl;
        Y = Y-mean(Y);
        
        Y = spline([1:sP], Y, linspace(1, sP, nPoints));
        T = spline([1:sP], T, linspace(1, sP, nPoints));
        ax.title = 'Reprocessed Data';
    else
        Y = Curves{Proc}.Data;
        dt = str2double(Curves{Proc}.XSlope);
        T = 0:dt:(numel(Y)-1)*dt;
        ax.title = 'Processed Data';        
    end
else
    if ~isempty(Curves{SinVal}.Data)
        Y = Curves{SinVal}.Data + 1i*Curves{CosVal}.Data;
        dt = str2double(Curves{SinVal}.XSlope);
        T = 0:dt:(numel(Y)-1)*dt;
        ax.title = 'Raw Data';
    else
        Y = Curves{Proc}.Data;
        dt = str2double(Curves{Proc}.XSlope);
        T = 0:dt:(numel(Y)-1)*dt;
        ax.title = 'Processed Data';        
    end
end

Bs = Curves{Bvals}.Data;
if numel(Bs) == 1
    
    B = linspace(Bs, Bs+str2double(Curves{Bvals}.XSlope), numel(Y)).';
else
    dt = str2double(Curves{Bvals}.XSlope) ;
    tt = 0:dt:(numel(Bs)-1)*dt + str2double(Curves{Bvals}.XOffset);
%     x = 1:numel(Bs);
%     xx = linspace(1, numel(Bs), numel(Y));
    B = spline(tt, Bs, T).';
end

% ax.x = B;
Bf = str2double(Pars.Bfrom_mT);
Bt = str2double(Pars.Bto_mT);
ax.x = linspace(Bf, Bt, numel(Y)).';
Y = spline(B, Y, ax.x);

ax.xlabel = 'mT';
ax.y = 1;

dsc = Pars;
dsc.size_x = num2str(size(Y,1));
dsc.size_y = num2str(size(Y,2));

% assign output depending on number of output arguments
% requested
switch nargout
  case 0, 
  figure(53); clf;
  plot(ax.x, Y);
  disp(dsc);
  case 1,
    varargout = {Y};
  case 2,
    varargout = {ax, Y};
  case 3,
    varargout = {ax, Y, dsc};
end

function Pars = readParams(str, Pars)
in = strfind(str, '>');
str = str(1:in(1));

iv = strfind(str, '"');
ie = strfind(str, '=');
is = strfind(str, ' ');

for jj = 1:numel(ie)
    ic = find(is-ie(jj)<0);
    Nam = str((is(ic(end))+1):(ie(jj)-1));
    idx1 = iv(2*jj-1)+1;
    idx2 = iv(2*jj)-1;
    Val = str(idx1:idx2);
    Pars = setfield(Pars, Nam, Val);
end

function ry = loc_rcfilt(y, scantime, convtime)
nP = numel(y);
dx = scantime/nP;
ff = linspace(-1/dx/2, 1/dx/2, nP).';
tc = convtime;
func = 1./(1+1i*2*pi*(ff)*tc);
fy = fftshift(fft(y));
ry = real(ifft(fftshift(fy.*func)));



