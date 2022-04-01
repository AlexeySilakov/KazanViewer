function varargout = kv_pseumod(varargin)
%   kv_pseumod  Pseudo modulation for Kazan Viewer
%
%     modY = kv_pseumod(x,y,ModAmpl);
%     modY = kv_pseumod(x,y,ModAmpl,[Harmonic], [Dimention]);
%  
%   Pseudomodulation - Convolution of the spectrum with bessel function
%  
%     Input:
%     - x: axis data
%     - y: Spectrum
%     - ModAmpl: peak-to-peak modulation amplitude [same units as x]
%     - Harmonic: harmonic, greater than 0, default is 1
%     - Dimention: if 2D, which dimention to pseudomodulate
%     Output;
%     - modSpec: pseudomodulated spectrum

x = varargin{1};
if sum (size(x)==sum(size(varargin{2})))
    error(['Dimentions of X [', num2str(size(x)) , '] do not fit dimentions of Y [',num2str(size(x)) , ']']);
end
y = varargin{2};
ModAmpl = varargin{3};

if nargin >=4, harm = varargin{4}; else harm = 1; end
if nargin ==5, dim = varargin{5}; else dim = 2; end

isvect = find(size(y)==1); % is empty == it is a 2D data set

if ~isempty(isvect)
    if isvect ~= 2 
        y = y.';
    end
end

if dim == 1 % do rows
    y = y.';
end

fs = fft(y);
mS = 1/(x(2)-x(1))/2;
%%%%%%%%%%%%%%%%%%%%%% alsi 30.08.2009 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hh = fftshift(besselj(harm, linspace(-1, 1, size(fs, 1)).'*mS*ModAmpl*pi));
% res = real((i)^harm*ifft(fs.*[hh(:, ones(size(y, 2),1)); y*0]));
res = real((i)^harm*ifft(fs.*hh(:, ones(size(y, 2),1)))/2);
if isvect ~= 2
     varargout{1} = res(1:size(y, 1), :).'; % get back to original column/row orientation
else
    varargout{1} = res(1:size(y, 1), :);
end

if nargout>=2
    disp('pseumod: too much output arguments');
end
% y = gaussian(x, 60, 10) + gaussian(x, 30, 10);
