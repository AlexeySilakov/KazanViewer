function varargout = autofit_moss(varargin)
% if nargin==0
clear global;
global noeff temp  DATA intData inpAll relpar W1 Amp DE;
[ax,y,dsc]=kv_d01read('D:\PSU\Data\Tempstudies_for_alexey\Cyp119\1_Cyp119-I\Cyp119_CompI-AllTemp.d01');

noeff.nPoints = length(ax.x);
noeff.Range = [min(ax.x), max(ax.x)];
noeff.g = [2, 1.9, 1.8];

DATA = y(:, 1) - mean(y(1:10, 1)); temp = 4.2;
DATA = DATA/sum(DATA);
% intData = sum( DATA);
x0(1) = 0.08;% minp.IS = 
x0(2) = 0.2735;% minp.lw = 
x0(3:5) = [-29.45, -32.10, -3.90]*1e4;% minp.A(1:3) = 
x0(6) = -10.88;% minp.Apa2 = 
x0(7) = 0.916;% minp.Deq = 
x0(8) = 0.0029;% minp.eta = 

DE = 60;
Amp = 4;
options = optimset('Display','off','MaxIter',10);
% for ii = 1:10
res = fminsearch(@lsq_pars, x0, options);
inpAll = res;
temp = ax.y(1);
plot_stuff;

pause(1);
x0 = [];
inpAll(2) = 0.1735; % lw
x0(1) = 4; % Amp
x0(2) = 10; % DE

for ii = 2:size(y, 2)
    DATA(:, ii) = y(:, ii) - mean(y(1:10, ii)); temp(ii) = ax.y(ii);
    DATA(:, ii) = DATA(:, ii)/sum(DATA(:, ii));
end
options = optimset('Display','off','MaxIter',500);
res = fminsearch(@lsq_rel, x0, options);
% disp(sprintf('Temperature, %3f', temp));
Amp = abs(res(1));
DE = abs(res(2));

     plot_stuff;
     
%     title(sprintf('Temperature, %3f', temp));
    disp(res);
%     pause(1);
%     x0 = res;
%     W1(ii) = res(2);
% 
% 
% assignin('base', 'W1', W1);
% assignin('base', 'T', ax.y);
% temp = ax.y';
% 
% % fit exponential 
% disp('making exponential fit')
% options = optimset('Display','off','MaxIter',5000);
% 
% x0 = [20, 60];
% res = fminsearch(@lsq_exp, x0, options);
% figure(444); plot(ax.y, W1); hold on;
% yt = res(1)*res(2)^3./(exp(res(2)/temp)-1);
% plot(ax.y, yt); 
% end
% else
%         [sol, val] = lsqlsq1(varargin{1},varargin{2});
%     varargout{1} = sol;
%     varargout{2} = val;
% end
function res = lsq_exp(inp)
global temp W1;
y = W1-inp(1)*inp(2)^3./(exp(inp(2)/temp)-1);
res= sum(y.^2);

function res = lsq_pars(inp)
global  DATA intData;

minp.IS = inp(1);
minp.lw = inp(2);
minp.A(1:3) = inp(3:5);
minp.Apa2 = inp(6);
minp.Deq = inp(7);
minp.eta = inp(8);

[x, y] = LocalMoss(minp, 0);
y = (y-y(1))/sum(y);

res = sum(  (DATA-y+y(1)).^2  );

function res = lsq_rel(inp)
global  DATA inpAll temp;

minp.IS = inpAll(1);
minp.lw = inpAll(2);
minp.A(1:3) = inpAll(3:5);
minp.Apa2 = inpAll(6);
minp.Deq = inpAll(7);
minp.eta = inpAll(8);
res = sum(abs(inp - abs(inp))); % punishment for negative values
inp = abs(inp);
for ii = 1:length(temp)
    W = inp(1)./( exp(inp(2)./temp(ii))-1 );
    [x, y] = LocalMoss(minp, W);
    res = res+sum(  (DATA(:, ii)-y+y(1)).^2  );
end

function plot_stuff

global  DATA inpAll relpar Amp DE temp;

minp.IS = inpAll(1);
minp.lw = inpAll(2);
minp.A(1:3) = inpAll(3:5);
minp.Apa2 = inpAll(6);
minp.Deq = inpAll(7);
minp.eta = inpAll(8);
disp(minp);

step = 0.03;
figure(142); clf; hold on;
tt = temp;
for ii = 1:size(DATA, 2)
    W = Amp./( exp(DE./tt(ii))-1 );
    temp = tt(ii);
    [x, y] = LocalMoss(minp, W);
    y = (y-y(1))/sum(y);
    
    plot(x, DATA(:, ii)+step*(ii-1), '.b'); hold on; plot(x, y+step*(ii-1), '-r');
end

function [x, y] = LocalMoss(inp, relpar)
global noeff temp;

S = 0.5;
Omega = 0.0;
nPoints = noeff.nPoints;
nKnots = 10;
RelaxType = -1;
Range = noeff.Range;
IS = inp.IS;
lw = inp.lw;
kType = 'par';
Field = 540;
Temp = temp;
RelaxPar = relpar;
TDeb = 200; 
g = noeff.g;
A = inp.A; %G
Apa = [0, inp.Apa2, 0];
Deq = inp.Deq;
eta = inp.eta;
D = 0; %K
E = 0; % K
Conv = 1.0; %
Qpa = [0, 0, 0];
ZFS4 = [0, 0];
DirectDim = 0.0;


path = 'd:\MATLABR2006a\work\Mossbauer\Schutz\gfort_comp\';

dV = abs(diff(Range))/nPoints;
fid = fopen([path, 'RELAX.PAR'], 'w');

fprintf(fid, 'RELAX             FILENAME WITHOUT EXTENSION\n');
fprintf(fid, 'DFX               TITLE\n');
fprintf(fid, '%2.1f \t SPIN\n', S);
fprintf(fid, '%2.1f \t SPIN SPIN FLUCTUATION PARAMETER\n',Omega); 
fprintf(fid, '%i \t NUMBER OF CHANNELS\n', nPoints);
fprintf(fid, '%i \t THETA INTEGRATION STEPS\n', nKnots);
fprintf(fid, '%i \t PHI   INTEGRATION STEPS\n', nKnots);
fprintf(fid, '%i \t TYPE OF RELAXATION [SPIN SPIN (1) OR SPIN LATTICE (0)]\n', RelaxType);
fprintf(fid, '0                IFPUN\n');
fprintf(fid, '0                IPRT\n');
fprintf(fid, '0                IIX\n');
fprintf(fid, '%5f \t DELTA V [mm/s PER CHANNEL]\n', dV);
fprintf(fid, '%5f \t ZERO VELOCITY CHANNEL [mm/s PER CHANNEL]\n', abs(Range(1))/dV );
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
cd(path);
% tic
system(['del RELAX.DAT']);
system(['RELAX.EXE']);
% toc
cd(aaa);

fid = fopen([path, 'RELAX.DAT'], 'r');
count = 0;
while fid<0 && count<10
%     pause(0.0001);
    fid = fopen([path, 'RELAX.DAT'], 'r');
    count = count+1;
end

if fid<0, error('can not read the output file'); end
try
    A = textscan(fid, '%f', nPoints*2);
end
fclose(fid);

x  =A{1}(1:2:end);
y  =A{1}(2:2:end);