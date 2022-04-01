
function  [x, fSpec] = simpleEPR(Sys,Exp,Opt)
% [x, Spec] = simpleEPR(Sys,Exp,Opt)
%    routine for fast simulation of S=1/2 EPR spectra 
%    with or withough HStrain
%    with or without HF splitting (only relatively small HFC will be calculated correctly).
%    X - in [mT]
%
% How it works:
%   The core program is "cwepr_an".
%   It is written in C and therefore need to be compiled (mex cwepr_an.c)
%   It calculates EPR spectrum according to the provided magnetic field
%   values of principal components of the EPR spectrum using analytical
%   solution. This program also takes care of gaussian or pseudovoigtian HStrain. 
%   
%   HF coupling spliting is calculated in the first order approxiamtion
%   by calculating EPR spectrum for each MI separatly
%
%   !!! Some cosy routins of easyspin are used, so be sure it is
%   installed prior using this progaram
%
% Input parameters (mostly compatible with Easyspin ones):
%   Sys.S - electronic spin (not esential, anyway only S=1/2 is taken)
%   Sys.g - el. g-values Sys.g = [gx, gy, gz]
%   Sys.HStrain - strain for field components of the spectrum
%   Sys.gStrain - strain for g-values, converted to Sys.HStrain
%   Sys.lw - linewith
%      could be a number (gaussian), or two numbers (voight, i.e convolution of gaussian and lorentzian)
%   Sys.psVoight - pseudoVoightian option, provided by the core program
%      lshape = psVoight*lorenzian + (1-psVoight)*gaussian
%   Sys.Nucs - specify nuclei e.g. "Sys.Nucs = '1H, 57Fe' "
%      alternatively, one can specify explicitely Sys.I 
%      (only taken if Sys.Nucs is not present)
%      Sys.I = [1/2, 1/2]
%   Sys.A  - hyperfine coupling constants 
%      (row of three values for each nucleus)   
%   Sys.Apa - angles for HF thensors. Building the projections for gx, gy, gz
%      either in degree [default] or in radians (specify Opt.AngUn = 'rad')
%      if not specified, colinear case is taken
%
%   Exp.mwFreq - MW Frequency [GHz], default 9.7
%   Exp.Range - range of field to calculate [mT], default [300, 400]
%   Exp.nPoints - number of experimental points, default 1000 
%   
%   Opt.AngUn - which units are specified for Sys.Apa
%       Opt.AngUn = 'deg' [default] - degree
%       Opt.AngUn = 'rad'  - radians
%
%   Alexey Silakov, Max Planck Institute for Bioinorganic Chemistry, 2011

% bohrM = 9.27400915e-24;
tS = safeget(Sys, 'S', 1/2);
tg = safeget(Sys, 'g', [2.1, 2.04, 2.00]);
tHStrain = safeget(Sys, 'HStrain', [0.0, 0.0, 0.0]);
taHStrain = safeget(Sys, 'aHStrain', [0.0, 0.0, 0.0]);
tgStrain = safeget(Sys, 'gStrain', [0.0, 0.0, 0.0]);
tSRatio = safeget(Sys, 'SRatio', tS*0+1);

Range = safeget(Exp, 'Range', [300, 400]); % mT
nPoints = safeget(Exp, 'nPoints', 1000);
mwFreq = safeget(Exp, 'mwFreq', 9.7); % GHz

x = linspace(Range(1), Range(2), nPoints).';
fSpec = x*0; 
for ss = 1:numel(tS)
    S = tS(ss);
    if size(tg, 1)==numel(tS)
        g = tg(ss, :);
    elseif size(tg, 1)==1
        g = tg;
    else
        error('size of g-matrix does not match number of spins specified');
    end
    
    if size(tHStrain, 1)==numel(tS)
        HStrain = tHStrain(ss, :);
    elseif size(tHStrain, 1)==1
        HStrain = tHStrain;
    else
        error('size of HStrain-matrix does not match number of spins specified');
    end
    
    if size(tgStrain, 1)==numel(tS)
        gStrain = tgStrain(ss, :);
    elseif size(tHStrain, 1)==1
        gStrain = tgStrain;
    else
        error('size of HStrain-matrix does not match number of spins specified');
    end
    
    if size(taHStrain, 1)==numel(tS)
        aHStrain = taHStrain(ss, :);
    elseif size(taHStrain, 1)==1
        aHStrain = taHStrain;
    else
        error('size of aHStrain-matrix does not match number of spins specified');
    end
    if sum(abs(aHStrain))
        nASubs = 11;
    else
        nASubs = 1;
    end
    
    if S>0.5, error('This program can only work with Sys.S=1/2'); end
    
    % g = safeget(Sys, 'g', [2.1, 2.04, 2.00]);
    %
    % HStrain = safeget(Sys, 'HStrain', [0.0, 0.0, 0.0]);
    % gStrain = safeget(Sys, 'gStrain', [0.0, 0.0, 0.0]);
    
    lw = safeget(Sys, 'lw', [1.0, 0.0]);
    
    const = 71.44773042495561; % mT/GHz
    Bizz = const*mwFreq./g;
    nSubs = 1;
    dB = x(2)-x(1);
    
    % ---------------- Nuclei begin
    Nucs = safeget(Sys, 'Nucs', '');
    HF_AngUn = safeget(Opt, 'AngUn', 'deg'); % or "rad"
    
    if strcmp(HF_AngUn, 'rad')
        AngFact = 1;
    else
        AngFact = pi/180;
    end
    
    if isempty(Nucs)
        I = safeget(Sys, 'I', []);
        gn = safeget(Sys, 'gn', []); % not really necessary
    else
        [I, gn] = nucdata(Nucs); %%% to do ... make independent from easyspin
    end
    
    if ~isempty(I)
        HF_A = safeget(Sys, 'A', [0, 0, 0]);
        
        if size(HF_A, 1)~=length(I),
            error(['N. sets of A = ', num2str(size(HF_A, 1)), ' while length of Sys.I or Sys.Nucs = ' num2str(length(I))]);
        end
        
        HF_Apa = safeget(Sys, 'Apa', HF_A*0);
        
        
        Bs = Bizz;
        for ii = 1:length(I)
            
            %%%%%%%%%%%%% not sure that this is completely right
            RM = local_erot(HF_Apa(ii, :)*AngFact);
            hfc = diag(HF_A(ii, :))*RM;
            hfc = sqrt(hfc(1, :).^2+hfc(2, :).^2+ hfc(3, :).^2);
            HF_F = hfc*const./g*1e-3;
            
            Bs = local_splitBs(Bs, I(ii), HF_F);
        end
        Bizz = Bs;
        nSubs = size(Bizz, 1);
    end
    
    % ---------------- Nuclei end
    
    extra(1) = safeget(Opt, 'nHStrain', 100); % hHStrain
    if sum(abs(HStrain))==0
        extra(1) = 1;
    end
    extra(2) = safeget(Exp, 'Harmonic', 0.0); % derivative
    extra(3) = safeget(Sys, 'psVoight', 0.5); % lor/gau
    
    Spec = x*0;
    for mm = 1:nSubs
        
        if nASubs>1
            for ii = 1:nASubs
                tA = aHStrain(ones(size(Bizz, 1), 1), :)*ii/nASubs;
                tAmps = exp(-(ii-1).^2/nASubs.^2*9);
                Spec = Spec + cwepr_an(x, sort(Bizz(mm, :))+tA, dB, HStrain, lw(1), extra)*tAmps;
            end
        else
            Spec = Spec + cwepr_an(x, sort(Bizz(mm, :)), dB, HStrain, lw(1), extra);
        end
    end
    if length(lw)>1
        if lw(2)~=0
            lw = lorentzian(x, mean(x), Sys.lw(2), 0);
            Spec = real(fft( ifft(Spec).*ifft( fftshift(lw))));
        end
    end
    
    fSpec = (2*Spec.*x/const/mwFreq)*tSRatio(ss) + fSpec; % scale due to operation in B0 space, not frequency space
end


%%% calculate splitting
function Bs = local_splitBs(Bs, I, HF_F)

tBs = [];
for ii = 1:size(Bs, 1)
    for proj = (-I):I
        tBs(end+1, 1:3) = Bs(ii, 1:3) + proj*HF_F;
    end
end
Bs = tBs;

function R = local_erot(ang)
ph = ang(1);
th = ang(2);
ps = ang(3);

D = [-sin(ph), cos(ph), 0; -cos(ph), -sin(ph), 0; 0, 0, 1];
C = [1 0, 0; 0, cos(th), sin(th); 0, -sin(th), cos(th)];
B = [sin(ps), -cos(ps), 0; cos(ps), sin(ps), 0; 0, 0, 1];

R = B*C*D;

