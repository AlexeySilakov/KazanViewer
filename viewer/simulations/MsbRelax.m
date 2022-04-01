function varargout = MsbRelax(varargin)

S = 1;
Aunit = 'G';
un.G = 1;
un.kG = 1e3;
un.T = 1e4;
un.MHz = 1e4*0.7264;
un.mm_sFe = 11.625*un.MHz;

Omega = 0.0;
nPoints = 256;
nKnots = 10;
RelaxType = 1;
Range = [-10, 10];
IS = 0.1;
lw = 0.3;
kType = 'perp';
Field = 40e3;
Temp = 1.0;
RelaxPar = 0.0000165;
TDeb = 200; 
g = [2.1, 2.04, 2.00];
A = [-200 -200 -65]*1e3; %G
Deq = 2.5;
eta = 0.3;
D = -8.0; %K
E = -1.6; % K
Conv = 1.0; %
Qpa = [0, 0, 0];
ZFS4 = [0, 0];
Apa = [0, 0, 0];
DirectDim = 0.0;

fpath = 'D:\Matlab_work\Mossbauer\Schutz\gfort_comp\';

% kind of global safeget
if nargin>=1
    for kk = 1:nargin
        if isstruct(varargin{kk})
            Datas = struct2cell(varargin{kk});
            Names = fieldnames(varargin{kk});
            for ii = 1:length(Datas)
                if isnumeric(Datas{ii})
                    eval([Names{ii}, '=[',num2str(Datas{ii}), '];']);
                else
                    eval([Names{ii}, '= ''',Datas{ii}, ''';']);
                end    
            end
        end
    end
end

if isfield(un, Aunit)
    scale = getfield(un, Aunit);
else
    scale = 1;
    disp(sprintf('WARNING (MsbRelax), unit "%s" for Sys.A is not supported. Check Sys.Aunit.', Aunit));
end

A = A*scale;
dV = abs(diff(Range))/(nPoints-1);
fid = fopen([fpath, 'RELAX.PAR'], 'w');

fprintf(fid, 'RELAX             FILENAME WITHOUT EXTENSION\n');
fprintf(fid, 'DFX               TITLE\n');
fprintf(fid, '%2.1f \t SPIN\n', S);
fprintf(fid, '%2.1f \t SPIN SPIN FLUCTUATION PARAMETER\n',Omega); 
fprintf(fid, '%i \t NUMBER OF CHANNELS\n', nPoints);
fprintf(fid, '%i \t THETA INTEGRATION STEPS\n', nKnots);
fprintf(fid, '%i \t PHI   INTEGRATION STEPS\n', nKnots);
fprintf(fid, '%i \t TYPE OF RELAXATION [SPIN SPIN (1),SPIN LATTICE (0), W (-1)]\n', RelaxType);
fprintf(fid, '0                IFPUN\n');
fprintf(fid, '0                IPRT\n');
fprintf(fid, '0                IIX\n');
fprintf(fid, '%5f \t DELTA V [mm/s PER CHANNEL]\n', dV);
fprintf(fid, '%5f \t ZERO VELOCITY CHANNEL [mm/s PER CHANNEL]\n', abs(Range(1))/dV+1);
fprintf(fid, '%5f \t ISOMERSHIFT [mm/s]\n', IS);
fprintf(fid, '%5f \t LINEWIDTH [mm/s]\n', lw);
if isnumeric(kType)
    Th = kType;
else
    if strfind(kType, 'par')
        Th = 0.0;
    else
        Th = 90.0;
    end
end
fprintf(fid, '%5f \t ANGLE BETWEEN H-FIELD AND GAMMA BEAM [DEGREES]\n', Th);
fprintf(fid, '%7f \t MAGNETIC FIELD [G]\n', Field);
fprintf(fid, '%5f \t TEMPERATURE [K]\n', Temp);
fprintf(fid, '%5f \t RELAXATIONPARAMETER\n', RelaxPar);
fprintf(fid, '%5f \t DEBYE TEMPERATURE [K]\n', TDeb);
fprintf(fid, '%5f \t GX\n', g(1));
fprintf(fid, '%5f \t GY\n', g(2));
fprintf(fid, '%5f \t GZ\n', g(3));
fprintf(fid, '%5f \t AX [G]\n', A(1));
fprintf(fid, '%5f \t AY [G]\n', A(2));
fprintf(fid, '%5f \t AZ [G]\n', A(3));
fprintf(fid, '%5f \t QUADRUPOL SPLITTING [mm/s]\n', Deq);
fprintf(fid, '%5f \t ASYMMETRY PARAMETER\n', eta);
fprintf(fid, '%5f \t D [K]\n', D);
fprintf(fid, '%5f \t  E [K]\n', E);
fprintf(fid, '%5f \t CONSTANT MULTIPLYING S-L MATRIX ELEMENTS\n', Conv);
fprintf(fid, '%5f \t ALPHA   EULER ANGLES BETWEEN G-TENSOR AND EFG\n', Qpa(1));
fprintf(fid, '%5f \t BETA\n', Qpa(2));
fprintf(fid, '%5f \t GAMMA\n', Qpa(3));
fprintf(fid, '%5f \t 4TH ORDER CUBIC ZFS TERM [K]\n', ZFS4(1));
fprintf(fid, '%5f \t 4TH ORDER AXIAL ZFS TERM [K]\n', ZFS4(2));
fprintf(fid, '0.0             INITIAL THETAVALUE\n');
fprintf(fid, '0.0             INITIAL PHIVALUE\n');
fprintf(fid, '90.0             FINAL THETAVALUE\n');
fprintf(fid, '90.0             FINAL PHIVALUE\n');
fprintf(fid, '%5f \t ALPHA EULER ANGLES BETWEEN G-TENSOR AND A-TENSOR\n', Apa(1));
fprintf(fid, '%5f \t BETA\n', Apa(2));
fprintf(fid, '%5f \t GAMMA\n', Apa(3));
fprintf(fid, '%5f \t DIMENSIONALITY OF PHONON DENSITY FOR DIRECT PROCESS\n', DirectDim);

fclose(fid);
aaa = cd;
cd(fpath)
% tic
system(['del RELAX.DAT']);
system(['RELAX.EXE']);
% toc
cd(aaa);

fid = fopen([fpath, 'RELAX.DAT'], 'r');
count = 0;
while fid<0 && count<10
%     pause(0.0001);
    fid = fopen([fpath, 'RELAX.DAT'], 'r');
    count = count+1;
end

if fid<0, error('can not read the output file'); end
try
    A = textscan(fid, '%f', nPoints*2);
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




