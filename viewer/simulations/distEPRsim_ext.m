function varargout = distEPRsim_ext(varargin)
Sys = varargin{1};
Exp = varargin{2};
if nargin>2
    Opt = varargin{3};
else 
    Opt = [];
end

distE = 0;
distD = 0;
N = 0;
spread = 3;
if isfield(Sys, 'Ddist')
    N = safeget(Sys, 'nDist', 11);
    rr = linspace(-spread, spread, N).'*Sys.Ddist;
    
    distD = [-1*rr, -1*rr, 2*rr]/3;
    amp = exp(-rr.^2/(2*Sys.Ddist)^2)/sqrt(2*pi)/Sys.Ddist;
end
if isfield(Sys, 'Edist')
    
    N = safeget(Sys, 'nDist', 11);
    rr = linspace(-spread, spread, N).'*Sys.Edist;
    
    distE = [-1*rr, 1*rr, 0*rr];
    amp = exp(-rr.^2/(2*Sys.Edist)^2)/sqrt(2*pi)/Sys.Edist;
end

if N==0
    [x, y] = EPRsim_ext(Sys, Exp, Opt);
else
    dist = distE+distD;
%     amp = ampE*ampD;
    tSys = Sys;
    arr = zeros(Exp.nPoints, N);
    for ii = 1:N
        tSys.D  = Sys.D+dist(ii, :);
        disp(tSys.D)
        [x, ty] = EPRsim_ext(tSys, Exp, Opt);
        arr(:, ii) = ty*amp(ii);
    end
    y = sum(arr, 2);
    if safeget(Opt, 'distPlot', 0)
        figure;
        plot(x, arr, 'r', x, y, '.b');
    end
    
end

varargout{1} = x;% G
varargout{2} = y;