function [varargout] = fftprocessing(ax, y, pars)
% FFT processing. Written for fftplugin of Kazan Viewer
%     [rax, ry]= fftprocessing(ax, y, pars)
%     [rax, ry, apwin]= fftprocessing(ax, y, pars)
%     
%  Depending from field 'pars' can prepare data for FFT
%  like cutting data, zerofilling, creating appodization window 
%  or/and make 1D FFT of 'y'.
%   'y' sould be column  or 2D matrix.
%   'ax' has to be structure ans include field ax.x (X axis).
%   'pars' has to be structure as well. Possible fields:
%      pars.opt = 'real'/'cta'/'qseq'
%      pars.rshift, pars.lshift / cutting data
%      pars.awin /if isfield, apodization window is applied to data 
%           'bla','bar','con', 'cos','ham','han','wel'
%           'exp','gau', 'kai'
%           'lor-gau'
%      pars.awidth
%      pars.ashift
%      pars.aalpha
%      pars.zerofill / zerofilling
%      pars.fft / if =1, FFT is performed
%      pars.phase0 / phase correction
%      pars.xscale / X axis scalling
%
% See also FFT, FFTSHIFT, APOWIN, SAFEGET
% NOTE:
%     This function uses 'apowin' from EasySpin, and 'safeget' from
%     KazanViewer package
% AlSi 26.09.2004

x = ax.x;
x1 = safeget(ax, 'y', 0);
rax =ax;
ndim = size(x,1);

opt = safeget(pars, 'opt', 'real');
if strcmp(opt, 'real') | strcmp(opt, 'cta')
   y = real(y);
end

% cut data arrays from the right
rshift = floor(safeget(pars, 'rshift', 0)+0.5);
if rshift > 0
  y1(rshift+1:ndim, :) = y(1:ndim-rshift, :);
  y1(1:rshift, :) = 0;
else
  %script{end+1} = sprintf('y1 = y;', rshift);
  y1 = y;
end

% cut data arrays from the left
lshift = floor(safeget(pars, 'lshift', 0)+.5);
if lshift > 0
  y1(1:ndim - lshift, :) = y1(lshift+1:end, :);
  y1(ndim - lshift:end, :) = 0;
end

dx = abs(x(2)-x(1));
ndx = abs(x(end)-x(1));

% apodization window 
if isfield(pars, 'awin'),
  width = abs(safeget(pars, 'awidth', .01))*ndx;
  npoints = floor(width/dx+.5);
  nshift = ndim - safeget(pars, 'ashift', 0)/dx + npoints;
  switch pars.awin
    case {'bla','bar','con', 'cos','ham','han','wel'}
      aw = apowin(pars.awin, npoints*2);
      lastv = aw(end);
      awin = [ones(ndim, 1)*lastv; aw; zeros(ndim, 1)*lastv];
    case 'exp'
      ttx = abs(linspace(-1, 1, ndim*2+npoints*2)).';
      %ttx0 = safeget(pars, 'ashift', 0);
      ttG = safeget(pars, 'awidth', 1);
      taa = safeget(pars, 'aalpha', 2);  
      awin = exp(ttx/ttG/ttaa);
      nshift = ndim+npoints;
    case 'gau'
      ttx = abs(linspace(-1, 1, ndim*2+npoints*2)).';
      ttx0 = safeget(pars, 'ashift', 0);
      ttG = safeget(pars, 'awidth', 1);
      ttaa = safeget(pars, 'aalpha', 2);
      awin = exp(-2*((ttx-ttx0)/ttG).^2);
      nshift = ndim+npoints;
    case 'kai'
      aw = apowin(pars.awin, npoints*2,safeget(pars, 'aalpha', 2));
      lastv = aw(end);
      awin = [ones(ndim, 1)*lastv; aw; zeros(ndim, 1)*lastv];
    case {'lor-gau'}
        ttx = abs(linspace(-1, 1, ndim*2+npoints*2)).';
        ttx0 = safeget(pars, 'ashift', 0);
        ttG = safeget(pars, 'awidth', 1);
        ttaa = safeget(pars, 'aalpha', 2);
      aw = exp(-2*((ttx-ttx0)/ttG).^2).*exp(ttx/ttG/ttaa);
      lastv = aw(end);
      awin = aw; %[ones(ndim, 1)*lastv; aw; zeros(ndim, 1)*lastv];
      nshift = ndim+npoints;
  end
  awin = awin([nshift:nshift+ndim-1]);
  y1 = y1.*awin(:, ones(size(y1, 2), 1)); 
else
awin = ones(size(y1, 1), 1);  
end
y2 = y1;
apwin = awin(:, ones(size(y1, 2), 1));
    
% zerofilling
ndim1 = ndim;
zfill = floor(safeget(pars, 'zerofill', 0)+0.5);
if zfill,
    npow = fix(log2(ndim)+.5); %AlSi (ceil -> fix)
    npow = 2^(abs(zfill)+npow);
    if zfill < 0
        zfill = abs(zfill);
        y1(npow/2+1:npow/2+ndim,:) = y1;
        y1(1:npow/2-1 , :)=0;
        y1(end+1:npow, :)=0;
    else
        y1(end+1:npow, :)=0;
    end
    ndim1 = npow;
end
rx = x;
ry = y2;
pars.fft = safeget(pars, 'fft',1);
if pars.fft == 1
    switch opt
    case 'qseq'  
        y1(2:2:end) = -y1(2:2:end);
        ytmp=zeros(1,2*ndim1);
        ytmp(1:2:2*ndim1)=real(y1);
        ytmp(2:2:2*ndim1)=imag(y1);
        ry = fft(real(ytmp));
        rx = [-ndim1:ndim1-1].'*1/dx/ndim1/2;
    case 'real'
        rx = [0:ndim1/2-1].'*1/dx/ndim1;
        ry = fftshift(fft(y1), 1); 
        ry = ry(ndim1/2+1:end, :);
    case 'cta', % cross-term averaging
        avpar = safeget(pars, 'cta', ndim/4);
        rx = [0:ndim1/2-1].'*1/dx/ndim1;
        ry = 0;
        for p=1:avpar,yy=y1(p:end,:);    
            yy(end+1:ndim1,:)=0;    
            ry = ry + abs(fftshift(fft(yy), 1)).^2;
        end
        ry = sqrt(ry(ndim1/2+1:end, :));
    otherwise
        rx = [-ndim1/2:ndim1/2-1].'*1/dx/ndim1;
        ry = fftshift(fft(y1), 1);
    end
    if isfield(pars, 'xscale'), 
        rx = rx * pars.xscale;
    end
    if isfield(pars, 'phase0')
        ry = ry.*exp(j*pars.phase0*pi/180);
    end
    [rax.xlabel,rax.x] = findunit(ax.xlabel, rx);
end
varargout{1} = rax;
varargout{2} = ry;
if nargin ==3
    varargout{3} = apwin;
end

%-----------------------------------------------------
function [newlabel,newaxx] = findunit(label, axx)
newlabel = label;
newaxx = axx;
p = findstr(label, ',');
if isempty(p), return; end
[un, unit, c1, c2] = kvgetvalue(['1',label(p+1:end)]);
[a,b,c] = kvbestunit(max(abs(axx))/un, '');
newaxx = axx/(c*c2);

switch unit
   case 's',  unit = 'Hz'; name = 'Frequency,';
   case 'Hz',  unit = 's'; name = 'Time,';
   otherwise, unit = ''; name = ['1/',label(1:p)];
end
newlabel = [name, ' ', b, unit];

