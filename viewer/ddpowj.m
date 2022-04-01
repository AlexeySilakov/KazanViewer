function varargout = ddpowj(varargin)
Sys = [];
Exp = [];
Opt = [];
if nargin>0, Sys =varargin{1}; end
if nargin>1, Exp =varargin{2}; end
if nargin>2, Opt =varargin{3}; end

fpath = 'd:\Install\Science\ddpowj\';

str = {};
% something               TITLE (<= 72 CHARACTERS)
str{end+1} = 'kv generated';
% 2.0						S(1) ENTER SPIN OF SITE 1
str{end+1} = sprintf('%-2.1f', Sys.S(1));
% 2.0						S(2) ENTER SPIN OF SITE 2
str{end+1} = sprintf('%-2.1f', Sys.S(2));
% 2.0 2.0 2.0						GXYZ(1) ENTER g VALUES FOR ELECTRONIC SPIN 1
str{end+1} = sprintf('%-8.5f %-8.5f %-8.5f', Sys.g(1, 1),Sys.g(1, 2), Sys.g(1, 3));

D = safeget(Sys, 'D', [0, 0]);
E = safeget(Sys, 'E', [0, 0]);
mu = safeget(Sys, 'mu', [0, 0]);
F = safeget(Sys, 'F', [0, 0]);
A = safeget(Sys, 'A', [0, 0]);
B = safeget(Sys, 'B', [0, 0]);
% 1.0						XD(1) ENTER 2nd ORDER AXIAL ZFS PARAMETER FOR SPIN 1 IN cm-1
str{end+1} = sprintf('%-8.5f', D(1));
% 0.0						XE(1) ENTER 2nd ORDER RHOMBIC ZFS PARAMETER FOR SPIN 1 IN cm-1
str{end+1} = sprintf('%-8.5f', E(1));
% 0.0						XMU(1) ENTER 3rd ORDER ZEEMAN PARAMETER FOR SPIN 1 ISOTROPIC, UNITLESS
str{end+1} = sprintf('%-8.5f', mu(1));
% 0.0						XF(1) ENTER 4th ORDER AXIAL ZFS PARAMETER FOR SPIN 1
str{end+1} = sprintf('%-8.5f', F(1));
% 0.0						XA(1) ENTER 4th ORDER CUBIC ZFS PARAMETER FOR SPIN 1
str{end+1} = sprintf('%-8.5f', A(1));
% 0.0						XB(1) ENTER 4th ORDER RHOMBIC ZFS PARAMETER FOR SPIN 1
str{end+1} = sprintf('%-8.5f', B(1));
% 2.0 2.0 2.0						GXYZ(2) ENTER g VALUES FOR ELECTRONIC SPIN 2
str{end+1} = sprintf('%-8.5f %-8.5f %-8.5f', Sys.g(1, 1),Sys.g(1, 2), Sys.g(1, 3));
% 0.0 0.0 0.0						ALPHAS(2), BETAS(2), GAMMAS(2) ENTER ROTATION ANGLES FOR ELECTRONIC SPIN 2
gpa = safeget(Sys, 'gpa', zeros(2, 3));
str{end+1} = sprintf('%-8.5f %-8.5f %-8.5f', gpa(2, 1),gpa(2, 2), gpa(2, 3));
% 1.0						XD(2)
str{end+1} = sprintf('%-8.5f', D(2));
% 0.0						XE(2)
str{end+1} = sprintf('%-8.5f', E(2));
% 0.0						XMU(2)
str{end+1} = sprintf('%-8.5f', mu(2));
% 0.0						XF(2)
str{end+1} = sprintf('%-8.5f', F(2));
% 0.0						XA(2)
str{end+1} = sprintf('%-8.5f', A(2));
% 0.0						XB(2)
str{end+1} = sprintf('%-8.5f', B(2));
% 0.1						XJI ENTER ISOTROPIC SPIN-EXCHANGE COUPLING (J, IN cm-1)
str{end+1} = sprintf('%-8.5f', safeget(Sys, 'J', 0));
% 0.0						XJJ ENTER 2nd ORDER SPIN-EXCHANGE COUPLING (J(2), ISOTROPIC, IN cm-1
str{end+1} = sprintf('%-8.5f', safeget(Sys, 'J2', 0));
% 0.0						XJD ENTER AXIAL PART OF DIPOLAR SPIN COUPLING (Ddip, IN cm-1)
str{end+1} = sprintf('%-8.5f', safeget(Sys, 'Jax', 0));
% 0.0						XJE NTER RHOMBIC PART OF DIPOLAR SPIN COUPLING (Edip, IN cm-1)
str{end+1} = sprintf('%-8.5f', safeget(Sys, 'Jrh', 0));
% 0.0	0.0 0.0				ALPHADIP, BETADIP, GAMMADIP ENTER ROTATION ANGLES FOR DIPOLAR COUPLING WITH RESPECT TO SPIN 1 (IN DEGREES)
str{end+1} = sprintf('%-8.5f ', safeget(Sys, 'Jpa', [0,0,0 ]));
% 0.0	0.0 0.0				XJA(1, 2, 3) ENTER ANISOTROPIC EXCHANGE COUPLING (JA (aka D), IN cm-1)
str{end+1} = sprintf('%-8.5f ', safeget(Sys, 'Ja', [0,0,0 ]));
% 0.0						XJQ ENTER ELECTRONIC QUADRUPOLE COUPLING (JQ (aka A), ISOTROPIC, IN cm-1)
str{end+1} = sprintf('%-8.5f ', safeget(Sys, 'Jq', 0));
% 9.4						FREQ ENTER MICROWAVE FREQUENCY (IN GHz)
str{end+1} = sprintf('%-8.5f', safeget(Exp, 'mwFreq', 9.4));
% 0						MODEH1 ENTER 0 FOR H1 PERPENDICULAR, NONZERO FOR H1 PARALLEL
str{end+1} = sprintf('%-d', strcmp(safeget(Exp, 'Mode', 'perp'), 'par'));
% 4.0						TEMP ENTER TEMPERATURE (IN K, DEFAULT IS 2)
str{end+1} = sprintf('%-8.5f', safeget(Exp, 'Temperature', 4.0));
% 1024					NPTS ENTER NUMBER OF SPECTRAL POINTS (MAX= 1024)
Exp.nPoints = safeget(Exp, 'nPoints', 1024);
str{end+1} = sprintf('%d', Exp.nPoints);
% 0						IGLOO ENTER IGLOO (H0) GRID INTERVAL FOR POWDER PATTERN
str{end+1} = sprintf('%d', safeget(Opt, 'nKnots', 10));
% 0						NRHO ENTER RHO (H1, H2) GRID INTERVAL (DEFAULT= 8)
str{end+1} = sprintf('%d', safeget(Opt, 'Nrho', 16));
% N  						ANGL DO YOU WANT TO ENTER H0 ANGLE RANGE (Y/N)
str{end+1} = 'N';
% 0						MSIMIN MINIMUM STATE LEVEL from WHICH TO CALCULATE EPR (0=min)
str{end+1} = '0';
% 100						MSIMAX MAXIMUM STATE LEVEL from WHICH TO CALCULATE (1000=max)
str{end+1} = '100';
% 0						MSJMIN MINIMUM STATE LEVEL to WHICH TO CALCULATE (0=min)
str{end+1} = '0';
% 100						MSJMAX MAXIMUM STATE LEVEL to WHICH TO CALCULATE (1000=max)
str{end+1} = '100';
% G						GORT ENTER UNITS FOR MAGNETIC FIELD (G/T)
str{end+1} = 'G';
% 0.0 4000.0				HMIN, HMAX ENTER MAGNETIC FIELD MINIMUM AND MAXIMUM
str{end+1} = sprintf('%-8.5f ', safeget(Exp, 'Range', [100, 400])*10);
% 1						IDERIV EPR LINE FORM (0 FOR ABSORPTION, 1 FOR FIRST DERIVATIVE)
str{end+1} = sprintf('%d', safeget(Exp, 'Harmonic', 0));
% G						ELNTYP ENTER EPR LINE TYPE (L FOR LORENTZIAN, G FOR GAUSSIAN)
str{end+1} = sprintf('%s', safeget(Exp, 'lshape', 'G'));
% 1.0 1.0 1.0				WXYZ(1, 2, 3) ENTER ANISOTROPIC EPR LINEWIDTH IN ORDER Wx, Wy, Wz; HWHM (MHz)
lw = safeget(Sys, 'lw', 1);
if length(lw)==3
    str{end+1} = sprintf('%-8.5f ', lw);
else
    str{end+1} = sprintf('%-8.5f ', lw([1, 1, 1]));
end
% 3						LECUTOF ENTER LINEWIDTH CUTOFF FOR EPR (# OF LINEWIDTHS)
str{end+1} = sprintf('%d', safeget(Opt, 'Threshold', 10));
% N   					SETWMS ... just say {n}
str{end+1} = 'N';
% test.out
str{end+1} = 'kv_ddpowj.out';
% test.dat
str{end+1} = 'kv_ddpowj.dat';
% test.txt
str{end+1} = 'kv_ddpowj.txt';

fid = fopen([fpath, 'kv_ddpowj.inp'], 'w');
for ii = 1:length(str)
    fprintf(fid, '%s \r\n', str{ii});
end
fclose(fid);

aaa = cd;
cd(fpath)
% tic
system(['del kv_ddpowj.dat']);
[status,status2]= system(['ddpowj.exe kv_ddpowj.inp']);

% toc
cd(aaa);
if ~isempty(strfind(status2, 'abnormal'))
    disp(status2);
    error('PROGRAM CRASHED');
end

fid = fopen([fpath, 'kv_ddpowj.dat'], 'r');
count = 0;
while fid<0 && count<10
%     pause(0.0001);
    fid = fopen([fpath, 'kv_ddpowj.dat'], 'r');
    count = count+1;
end

if fid<0, error('can not read the output file'); end
try
    A = textscan(fid, '%f', min(Exp.nPoints, 1024)*2);
end
fclose(fid);

x  =A{1}(1:2:end);
y  =A{1}(2:2:end);

if ~nargout
    figure(12); clf; plot(x, y);
else
    varargout{1} = x;
    varargout{2} = y;
end
