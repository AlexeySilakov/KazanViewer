function varargout = ODF3(varargin)

ODFpath = 'D:\Programs\ODF3\kv_interface';
SysFile = 'ODF3_sys.par';
OPTFile = 'ODF3_opt.par';

if nargin==0
    Sys.S = 1/2;
    Sys.Nucs = '1H, 1H';
    Sys.A(1, :) = [1, 1,2 ];
end

nVar = 0;
numave = 1;
Ctrl = 0;
SimType = 1; % Type of spectrum simulation 1 -  coinciding axes of A- and g-tensors 
% 2  -  tilted  A-tensor,  nuclear  Zeeman,  forbidden 
% transitions,  3  -  similar  to  type  2,  but  restricted  by 
% the  number  of  transitions  according  to  level 
% probability,  4  -  similar  to  2,  probability  of 
% transitions is taken from precalculated file *.trn; 
% the second value is not used 
gau = 1;
lor = 0;
nNucs = 1;
 
keyODF = 1; % Type of axial ODF:  
% keyODF=1  -  general  form  is  described  by  expression 
% (1),  
% keyODF=2  -  hidden  uniaxiality  of  R  is  described  by 
% expressions (2), 
%  keyODF=3  -  ordered  uniaxiality  is  described  by 
% expressions (3) 

ErrMode = 1; % 1 - using Hessian and Jacobian, 2 - using Hessian only, 
% 3 - using Jacobian only. 
% The second value is not used. 

fid = fopen('nameODF3.nam', 'w');
fprintf(fid, '1\n')
fprintf(fid, '%s,%s, no, 1\n', SysFile, OPTFile);
fprintf(fid, 'out_%s,out_%s\n', SysFile, OPTFile);
fprintf(fid, 'exp.esr,res.dat,0,1\n', SysFile, OPTFile);
fprintf(fid, 'end\n', SysFile, OPTFile);
fclose(fid)




str{1} = sprintf('%d', Ctrl); %0 - continue, 1 - stop iterations 
str{2} = sprintf('%d, %d', nVar, numel(Sys.S)); % 
str{3} = sprintf('%d, %d', Exp.nPoints, Opt.nKnots); % 
str{4} = sprintf('%d, %d', numave, safeget(Opt, 'Lib', 0)); % Numave, 0-3 Libration (?) model
str{5} = '1, 1'; % phase of experimental and calculated EPR spectra  
str{6} = sprintf('%d, 0', SimType);
str{7} = sprintf('%d, %d', gau, lor); 
str{8} = sprintf('%d, 0', nNucs);
str{9} = '2, 0';%2 - use FFT for convolution; 
% 0  -  key  for  use  of  weighted  points  of  spectra,  the 
% weights are taken from a separate file. 

t = sprintf('%2.1f, ',Sys.I); 
str{10} = t(1:(end-2));

str{11} = sprintf('%7.2f, %7.2f, 1.0000', Exp.Range(1), Exp.Range(2)); 
str{12} = '1.0D-8, 1.0D-6, 1.0D-6, 1.0D-5, 0.0';%Step  for  calculation  of  gradients,  RFCTOL  - 
% function  convergence  level,  XCTOL  -  arguments 
% convergence level, probability level for transitions, 
% last value is not used 
str{13} = '  flag';
str{14} = '66.0   0  0.31   0  1.0   0 '% ?????
str{15} = sprintf('%7f 0 %7f 0 %7f 0', gg(1), gg(2), gg(3)); 
str{16} = sprintf('0.0 0 0.0 0 0.0 0'); %%% rotate gg to g-tensor... useless
str{17} = sprintf('%7f 0 %7f 0 %7f 0', ll(1), ll(2), ll(3)); 
str{18} = sprintf('0.0 0 0.0 0 0.0 0'); %%% rotate gg to g-tensor... useless
str{19} = sprintf('%7f 0 %7f 0 %7f 0', Sys.g(1), Sys.g(2), Sys.g(3)); 
str{20} = '0.0     0  0.0     0  0.0   0'; % Amplitudes of librations around gx, gy, gz axes 
str{21} = '0.0     0  0.0     0  0.0   0'; % Euler  angles  that  transform  libration  frame  to  g-frame, in degrees 
str{22} = '0.0     0  0.0     0  0.0   0';% betta  coefficients  (before  m)  for  Gaussian  line width 
str{23} = '0.0     0  0.0     0  0.0   0';% gamma  coefficients  (before  m2)  for  Gaussian  line width 
str{24} = '0.0     0  0.0     0  0.0   0';% betta  coefficients  (before  m)  for  Lorentzian  line width 
str{25} = '0.0     0  0.0     0  0.0   0';% gamma  coefficients  (before  m2)  for  Lorentzian  line width 
for ii = 1:nNucs
    str{end+1} = sprintf('%7f 0 %7f 0 %7f 0', Sys.A(ii, 1), Sys.A(ii, 2), Sys.A(ii, 3)); 
    str{end+1} = sprintf('%7f 0 %7f 0 %7f 0', Sys.Apa(ii, 1), Sys.Apa(ii, 2), Sys.Apa(ii, 3)); 
    str{end+1} = sprintf('%7f 0', Sys.gn(ii)); 
end

fid = fopen(SysFile, 'w');
for ii = 1:length(str)
    fprintf(fid, '%s\n', str{ii})
end
fclose(fid);

strp{1} = '0'; %Rank of expansion 
strp{2} = '0 1.000'; %Penalty flag: 0 - penalty off, ? 0 penalty on; Penalty coefficient 
strp{3} = sprintf('%d', keyODF);   
strp{4} = sprintf('%d 0.0', ErrMode);   
strp{5} = 'Parameters: ';  
strp{6} = '0.0       0 ';  
strp{7} = 'Uniaxial: ';  
strp{8} = sprintf('%7f 0 %7f 0 %7f 0', modelAng(1), modelAng(2), modelAng(3)); % ???? not explained !, only used with keyODF=2or3
strp{9} = '  An0          Fl  An2         Fl   An4        Fl  An6         Fl  An8       Fl';
strp{10} = ' 1.0  0';
strp{11} = 'General: ';  
strp{12} = sprintf('%7f 0 %7f 0 %7f 0', modelAng(1), modelAng(2), modelAng(3)); % ???? not explained !, only used with keyODF=2or3
strp{13} = '      An0          Fl  An1        Fl   Bn1        Fl   An2       Fl   Bn2    Fl';
strp{14} = '0 1.0  0';

fid = fopen(OPTFile, 'w');
for ii = 1:length(strp)
    fprintf(fid, '%s\n', strp{ii})
end
fclose(fid);


aaa = cd;
cd(ODFpath)
% tic
system(['del res.dat']);
system(['ODF3_1.exe']);
% toc
cd(aaa);

fid = fopen([ODFpath, 'res.dat'], 'r');
count = 0;
while fid<0 && count<10
%     pause(0.0001);
    fid = fopen([ODFpath, 'res.dat'], 'r');
    count = count+1;
end

if fid<0, error('can not read the output file'); end
try
    A = textscan(fid, '%f', nPoints*2);
end
fclose(fid);




