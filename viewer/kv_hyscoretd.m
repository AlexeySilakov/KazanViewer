function [y, ax] = kv_hyscorefd(Sys, Exp, Opt, Shift)
% KV_HYSCOREFD powder HYSCORE simulation program.
% Is used as a part of the 'hyscorebox' plugin for the 'Kazan Viewer',
% However, can be used as a standalone script for executing in MatLab
% Some EasySpin routines are required
%
% [y, ax] = kv_hyscorefd(Sys, Exp, Opt, Shift)
% 
% input:
% Sys.A/Q = [MHz], Exp.mwFreq = [GHz], Exp.Field = [mT], Exp.tau = [ns]
% Sys.Apa/Qpa = [degree]
% output
% y - 2D matrix
% ax.x,y - axes in [Hz]
% Example: 
%   clear Sys Exp Opt Shift;
%   Sys = struct('S', 0.5, 'g', [2.0006, 2.006, 2.0654]);
%   Sys.I = [0.5; 0.5];  % there are two nuclei in the system (writen as a colomn)
%   Sys.gn = 0.1806; % both nuclei has the same g-factor (this equal two [0.1806; 0.1806]);
%   Sys.lw = 0.12;  % Unitary linewith (always the same for all nuclei)
%   Sys.A(1, :) = [-2, -2, 1.71]; Sys.A(2, :) = -1*[4 6.2 2.0]; % HF coupling [MHz]
%   Sys.Apa(1, :) = [0,28, 30];   Sys.Apa(2, :) =  [0,50, 145];  % Orientation of the HF tensor along g-tensor [degree]
%   Exp.gField = 2.0042; % the Magn. field automaticaly calculated as: 
%                        % Exp.Field = (Exp.mwFreq*1e+9*planck)/(Exp.gField*bmagn*1e-3)
%   Exp.mwFreq = 33.8;    % MW frequency [GHz]
%   Exp.ExciteWidth = 40; % same as in EasySpin, exitation width [MHz]
%   Exp.tau = 340;        % delay between first and second MW pulses [ns]
%   Exp.nPoints = 256;    % number of points
%   Opt.Symmetry = 'Ci'; Opt.nKnots = 20; Opt.Verbosity = 1;
%   Opt.KillNeg = 1;      % flag for suppression of negative time in timedomain spectrum
%   Opt.ProdRule = 1;      % use product rule, producing combination frequencies (memory/time consuming)
%   Opt.AmpRatios = [1; 2]; % relative amplitudes of the peaks from different HF couplings
%   Shift.ding = 0.2;
%   [y, ax] = kv_hyscorefd(Sys, Exp, Opt, Shift);
% 
%   % plot
%   ax.contour = [0.1:0.1:1];
%   ax.x = ax.x/1e6; % Hz ->MHz
%   ax.y = ax.y/1e6; % Hz ->MHz
%   ax.xlabel = 'Frequency, MHz';ax.ylabel = 'Frequency, MHz';
%   ax.Color = [0, 0, 0]; % color of lines in the "skyline" projections
%   kvskyline(100, ax, y, 'plot', 'image', 'hold', 'off', 'style', 1);

tstart = cputime;
tau = safeget(Exp, 'tau', 0)*1e-9; %ns->s

%%% simulation for a few nuclei
Parameters ={'A', 'Apa', 'Q', 'Qpa', 'Nucs', 'gn', 'I', 'nNucs'};

Sys.Apa = safeget(Sys, 'Apa', [0, 0, 0])*pi/180; % degree -> rad
Sys.Qpa = safeget(Sys, 'Qpa', [0, 0, 0])*pi/180; % degree -> rad
% lets find number of simmulaitons (parameters can be common, it means
% size(XX, 1)==1)
Sys.nNucs = safeget(Sys, 'nNucs', 1); % number of equivalent nuclei
sim_size = 1;
isfield_num = [];
for pk = 1:length(Parameters)
    if isfield(Sys, Parameters{pk})        
        tPar = getfield(Sys, Parameters{pk});
        isfield_num(end+1) = pk;        
        if size(tPar, 1)>1
            sim_size = size(tPar, 1);
        end
    end
end
if sim_size>1
    if size(Sys.I, 2) > 1,
        Sys.I = Sys.I.';
    end
    if size(Sys.gn, 2) > 1,
        Sys.gn = Sys.gn.';
    end
    for kk = 1:length(isfield_num)
        pk = isfield_num(kk); % kind of optimisation
        tPar = getfield(Sys, Parameters{pk});
        if size(tPar, 1)==1
            Sys = setfield(Sys, Parameters{pk}, tPar(ones(sim_size, 1), :));
        elseif size(tPar, 1)~=sim_size
%             error({['Parameter Sys.', Parameters{pk}, ' has wrong number of rows.']; ['Correct number is ', num2str(sim_size)]});
              error(['Parameter Sys.', Parameters{pk}, ' has wrong number of rows.'], ['Correct number is ', num2str(sim_size)]);
        end
    end
end

if isfield(Sys, 'Aiso'),
    Sys.A = safeget(Sys, 'A', 0) + Sys.Aiso;
end
% 
%%% orientation selection

%%%%%%%%%%%%%%%%%%%% new version 18.02.2006 %%%%%%%%%%%
Exp.ExciteWidth = safeget(Exp, 'ExciteWidth', 0);
Exp.mwFreq = safeget(Exp, 'mwFreq', 0);
if isfield(Exp, 'gField'),
   if ~Exp.mwFreq & ~isfield(Exp, 'Field'), ...
           error('ESEEM: there is no Exp.mwFreq to calculate B0');
   end
   Exp.Field = (Exp.mwFreq*1e+9*planck)/(Exp.gField*bmagn*1e-3);  % mT
else 
   Exp.Field = safeget(Exp, 'Field', 0);
end

[phi, theta, weigths] = Local_orisel(Sys, Exp, Opt);

% Weights = orisel(Sys,Par,Ori)

if isfield(Sys, 'Nucs'), [Sys.I,Sys.gn] = nucdata(Sys.Nucs); end % EasySpin function

Opt.Treshold = safeget(Opt, 'Treshold', 1e-3);

oriind = find(weigths >Opt.Treshold);
Exp.phi = phi(oriind, 1);
Exp.theta = theta(oriind, 1);
Exp.Weights = weigths(oriind, 1)/max(weigths(oriind, 1));

% draw sphere showing orientation selection
if safeget(Opt, 'ShowOri', 0)
    figure(126);
    [xx, yy, zz] = sphere(20);
    surf(xx, yy, zz, 'Linestyle', 'none'); %zz*0)
    colormap([0.5, 0.5, 1])
    shading interp;
    hold on;
    X =ang2vec(phi(oriind, 1), theta(oriind, 1));
    www = weigths(oriind, 1);
    for ii = 1:length(X(1, :))
        plot3(X(1, ii), X(2, ii), X(3, ii), 'Marker', '.', 'LineStyle', 'none', 'Color', [www(ii), 0, 0]);
        hold on;
    end
  %  X1 =ang2vec(Exp.phi, Exp.theta);
  %  plot3(X1(1, :), X1(2, :), X1(3, :), '.r');
    hold off;
end
% 
% MatrA = zeros(safeget(Exp, 'nPoints', 500));
% MatrB = MatrA;
% start calculation
AmpRation = safeget(Opt, 'AmpRatios', ones(sim_size, 1));
if (length(AmpRation)~=sim_size ),
    if length(AmpRation)~=1
        error('length(AmpRatios)~=size(A, 1)'); 
    else
        AmpRation = ones(sim_size, 1)*AmpRation;
    end
end

Ex = Exp;
Freqs = [];
Amps = [];
nNucs = Sys.nNucs;
make_prod = safeget(Opt, 'ProdRule', 1);
if sum(nNucs)==1, make_prod = 0; end
% nNucs - number of equivalent nuclei (column)
% all possible combinations for the product rule
nsim = sum(nNucs);

idxrange = [];
if ~isfield(Opt, 'Nuclei')
    Opt.Nuclei = 1:sim_size;
elseif max(Opt.Nuclei)>sim_size
    error('max(Opt.Nuclei) > number of nuclei');
end

% for ci =1:sim_size
%     if ~Opt.Nuclei(ci), continue; end
%     idxrange = [idxrange; ones(nNucs(ci), 1)*ci];    
%     
% end
nidxrange = Opt.Nuclei.'; % has to be a column

tnsim = length(nidxrange);
comb = [];
 for i=1:tnsim-1
%     if (~Opt.Nuclei(i))|(~Opt.Nuclei(i+1)), continue; end
    tempout=nidxrange(i+1:tnsim);
    comb=[comb; nidxrange(i)*ones(size(tempout,1),1)	tempout];
end
if isempty(comb), make_prod = 0; end

%%%%%%%%%%%%%%%%%%% optimisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. generation of Iij inside the cycle makes simulation twice longer
% 2. creating indexes here for crosspeaks calculations should make it 40% faster
% so, take it out of the cycle
global sIx sIy sIz
sIx = {};
sIy = {};
sIz = {};
maxmult = 2*max(Sys.I)+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Opt.Nuclei %%%%%%%%%%%%%%%%%%

for ii  =1:sum(nNucs)
    mult = 2*Sys.I(ii) + 1;

    sIx{mult} = local_sop([Sys.I(ii)], 'x');
    sIy{mult} = local_sop([Sys.I(ii)], 'y');
    sIz{mult} = local_sop([Sys.I(ii)], 'z');

    sIxx = sIx{mult}*sIx{mult};
    sIxy = sIx{mult}*sIy{mult};
    sIxz = sIx{mult}*sIz{mult};

    sIyx = sIy{mult}*sIx{mult};
    sIyy = sIy{mult}*sIy{mult};
    sIyz = sIy{mult}*sIz{mult};

    sIzx = sIz{mult}*sIx{mult};
    sIzy = sIz{mult}*sIy{mult};
    sIzz = sIz{mult}*sIz{mult};

    %         HF coupling in angle-dependent components
    A = diag(Sys.A(ii, :))*1e6; % [Hz]
    RM = erot(Sys.Apa(ii, :));
    A = RM*A*RM.';

    HFIx{ii} = (sIx{mult}*A(1, 1) + sIy{mult}*A(2, 1) + sIz{mult}*A(3, 1));
    HFIy{ii} = (sIx{mult}*A(1, 2) + sIy{mult}*A(2, 2) + sIz{mult}*A(3, 2));
    HFIz{ii} = (sIx{mult}*A(1, 3) + sIy{mult}*A(2, 3) + sIz{mult}*A(3, 3));

%       Quadrupole coupling 
    Quad{ii} = 0;
    if isfield(Sys, 'Q')&~(Sys.I==1/2);
        if sum(abs(Sys.Q(ii, :))) ~= 0

            Q = diag(Sys.Q(ii, :))*1e6; %[Hz]
            if isfield(Sys, 'Qpa'),
                RMQ = erot(Sys.Qpa(ii, :));
            else
                RMQ = diag(ones(3, 1));
            end
            Q = RMQ*Q*RMQ.';
            Quad{ii} = Q(1,1)*sIxx + Q(1, 2)*sIxy + Q(1, 3)*sIxz +...
                       Q(2,1)*sIyx + Q(2, 2)*sIyy + Q(2, 3)*sIyz +...
                       Q(3,1)*sIzx + Q(3, 2)*sIzy + Q(3, 3)*sIzz;
        end
    end

    % generating all possible combination of 6 numbers in range [1:nst]
    % made to avoid multiple "for" cycles (another optimisation)
    idxrange = 1:mult;

    b1(:, 1, 1, 1) = idxrange;
    b2(1, :, 1, 1) = idxrange;
    b3(1, 1, :, 1) = idxrange;
    b4(1, 1, 1, :) = idxrange;
    b5(:, 1) = idxrange;
    b6(1, :) = idxrange;

    nstones = ones(mult, 1);
    tii = reshape(b1(:, nstones, nstones, nstones), mult^4, 1);
    tkk = reshape(b2(nstones, :, nstones, nstones), mult^4, 1);
    tll = reshape(b3(nstones, nstones, :, nstones), mult^4, 1);
    tnn = reshape(b4(nstones, nstones, nstones, :), mult^4, 1);
    
    tjj = reshape(b5(:, nstones), mult^2, 1);
    tmm = reshape(b6(nstones, :), mult^2, 1);
    
    cii{ii} = tii;
    ckk{ii} = tkk;
    cll{ii} = tll;
    cnn{ii} = tnn;
    cjj{ii} = tjj;
    cmm{ii} = tmm;
    
    
    clear b1 b2 b3 b4 b5 b6 tii tkk tll tnn tjj tmm;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Amps = [];  Freqs = [];
Exp.nPoints = safeget(Exp, 'nPoints', 500);
Exp.MaxFreq = safeget(Exp, 'MaxFreq', 10)*1e6;
% SubtractLater = sqrt(-1);
MatrA =  zeros(Exp.nPoints);
dx = 1/Exp.MaxFreq/2;
Tax = linspace(0, dx*(Exp.nPoints-1), Exp.nPoints);
Ax1 = Tax(ones(Exp.nPoints, 1), :);
Ax2 = Ax1.';
for ang = 1:length(Exp.phi) %% check
    Ex.phi = Exp.phi(ang);
    Ex.theta = Exp.theta(ang);
%     AmpsA = []; AmpsB = [];
%     Freq1 = []; Freq2 = [];
    aAmps = [];  aFreqs = [];
%     for ci = 1:sim_size  
    for ci = Opt.Nuclei
        %%%%%%%%if Opt.Nuclei(ci)=0, do not include in the calculation
             
        Sy = Sys;
        for kk = 1:length(isfield_num)
            pk = isfield_num(kk);
            tPar = getfield(Sy, Parameters{pk});
            Sy = setfield(Sy, Parameters{pk}, tPar(ci, :));
        end  
         [HamA, HamB] = Local_spinham_opt(Sy, Ex, HFIx{ci}, HFIy{ci},HFIz{ci}, Quad{ci});

        [V_a, tD_a] = eig(HamA);  D_a(:, 1) = diag(tD_a); 
        [V_b, tD_b] = eig(HamB);  D_b(:, 1) = diag(tD_b);
        Str.V_a = V_a;
        Str.V_b = V_b;
        Str.E_a = D_a;
        Str.E_b = D_b;
        Str.vi = nmagn*Sy.gn*Exp.Field/planck;
        Str.I = Sy.I;

% %%%%%%%%%%%%%%%%%%%% -------------- test
% wa = real(D_a(:, ones(size(D_a, 1), 1)) - D_a(:, ones(size(D_a, 1), 1)).');
% wb = real(D_b(:, ones(size(D_b, 1), 1)) - D_b(:, ones(size(D_b, 1), 1)).');
% 
% W = [wa(:), wb(:)];
% Amp = W*0+1;
% %%%%%%%%%%%%%%%%%%%% -------------- test
        [Amp, W] = Local_crosspeak(Str, tau, cii{ci}, ckk{ci}, cll{ci}, cnn{ci}, cjj{ci}, cmm{ci});
        ampl{ci} = Amp*Exp.Weights(ang)*AmpRation(ci);
        frcs{ci} = W;
%         aAmps = [aAmps; Amp*Exp.Weights(ang)*AmpRation(ci)];
%         aFreqs = [aFreqs; W];
%          AmpsA = [AmpsA; Amp(:, 1)];% alhpa monifold
%          AmpsB = [AmpsB; Amp(:, 2)];% alhpa monifold         
%          Freq1 = [Freq1; W(:, 1)];
%          Freq2 = [Freq2; W(:, 2)];
    for jj =  1:size(Amp, 1)
        MatrA = MatrA + Amp(jj, 1)*exp(-2i*pi*(Ax1*W(jj, 1) + Ax2*W(jj, 2)))*Exp.Weights(ang)*AmpRation(ci)...
                      + Amp(jj, 2)*exp(-2i*pi*(Ax1*W(jj, 2) + Ax2*W(jj, 1)))*Exp.Weights(ang)*AmpRation(ci);
    end
    end
    if make_prod
        % product rule, only first product is counted, nothing like
        % exp(w1a*t1)*exp(w2a*t1)*exp(w3a*t1)... 
        disp('Product rule is not implemented yet');
%         for ind_c = 1:length(comb(:, 1))        
%             ci = comb(ind_c, 1);
%             cj = comb(ind_c, 2);
%                 sz1 = size(frcs{ci}, 1); 
%                 sz2 = size(frcs{cj}, 1);
%                 combination frequencies
%                 tDiff = frcs{ci}(:, ones(sz2, 1)) - frcs{cj}(:, ones(sz1, 1)).'; % w_1
%                 tDiff1 = frcs{ci}(:, ones(sz2, 1)+1) - frcs{cj}(:, ones(sz1, 1)+1).'; % w_2
%                 tSumm = frcs{ci}(:, ones(sz2, 1)) + frcs{cj}(:, ones(sz1, 1)).';
%                 tSumm1 = frcs{ci}(:, ones(sz2, 1)+1) + frcs{cj}(:, ones(sz1, 1)+1).';                
%                 amplitudes
%                 tAmps = ampl{ci}(:, ones(sz2, 1)).*ampl{cj}(:, ones(sz1, 1)).'*make_prod; % alpha
%                 tAmps1 = ampl{ci}(:, ones(sz2, 1)+1).*ampl{cj}(:, ones(sz1, 1)+1).'*make_prod; % beta 
%                 output
%                 aAmps = [aAmps; [tAmps(:), tAmps1(:)]; [tAmps(:), tAmps1(:)]];
%                 aFreqs = [aFreqs; [tDiff(:), tDiff1(:)]; [tSumm(:), tSumm1(:)]];
%         end
    end


end    
% disp(cputime - tstart);

% MatrA = make_fftcomb_bin(Freqs, Amps, Exp, kill_neg, safeget(Sys, 'lw', 0)*1e6/dx);  %alsi 08.06.06

MatrA = make_fft(real(MatrA), Exp, safeget(Sys, 'lw', 0)*1e6);

ax.x = linspace(-Exp.MaxFreq, Exp.MaxFreq, Exp.nPoints).'; % Hz
ax.y = linspace(-Exp.MaxFreq, Exp.MaxFreq, Exp.nPoints).'; % Hz

y = MatrA/max(max(MatrA)); 

switch safeget(Opt, 'Verbosity', 0)
    case 1
        disp('----kv_hyscorefd----')
        disp(['Calucation time (s): ',num2str(cputime-tstart)]);
    case 2
        disp('----kv_hyscorefd----')
        disp(['Calucation time (s): ', num2str(cputime-tstart)]);
        disp(['Maximum amplitude: ', num2str(max(max(real(y))))]);
    case 0
end

if safeget(Shift, 'ding', 0)
    x = [1:1000]/1.5;
    sound(min(Shift.ding, 1)*sin(x).*sin(x/2).*exp(-x/1000));
end

%----------------------------------------------------
function [MatrA, MatrB] = make_fftcomb(tMatrA, tMatrB, kill_neg, make_prod, lw)
% tMatrA, tMatrB  = [Exp.Dim, Exp.Dim, sim_size] (see main function)
% kill_neg = 1 or 0. If 1, negative part of the timedomain "spectra" will
%   be == 0
% make_prod amplitude of the combination freq. peaks. if 0, product rull
%   will be screwed
% lw  unitary linewidth IN POINTS, not in a MHZ !!!

fftMatrAt = (ifft2(fftshift(tMatrA/max(max(max(tMatrA))))));
fftMatrBt = (ifft2(fftshift(tMatrB/max(max(max(tMatrB))))));
% making linewidth
if lw
    fftGaga = real(ifft2(fftshift(Local_gaga2(1:size(tMatrA, 1), size(tMatrA, 1)/2+1, 1:size(tMatrA, 1), size(tMatrA, 1)/2+1, [1, 1]*lw))));
    fftMatrAt = fftMatrAt.*fftGaga(:, :, ones(size(tMatrA, 3), 1))/max(max(fftGaga));
    fftMatrBt = fftMatrBt.*fftGaga(:, :, ones(size(tMatrB, 3), 1))/max(max(fftGaga));
end
% killing revers time domains
if kill_neg
    MS = floor(size(tMatrA, 1)/2+0.5);
    fftMatrAt(MS:end, MS:end, :) = 0;
    fftMatrAt(1:MS, MS:end, :) = 0;
    fftMatrAt(MS:end, 1:MS,  :) = 0;
    
    fftMatrBt(MS:end, MS:end, :) = 0;
    fftMatrBt(1:MS, MS:end, :) = 0;
    fftMatrBt(MS:end, 1:MS,  :) = 0;

    tMatrA = (fftshift(fft2(fftMatrAt)));
    tMatrB = (fftshift(fft2(fftMatrBt)));    
end
% product rule
if make_prod
    MatrA  = abs(fftshift(fft2(prod(fftMatrAt, 3))))*make_prod + abs(fftshift(fft2(sum(fftMatrAt, 3))));
    MatrB  = abs(fftshift(fft2(prod(fftMatrBt, 3))))*make_prod + abs(fftshift(fft2(sum(fftMatrBt, 3))));
else
    MatrA  =abs(fftshift(fft2(sum(fftMatrAt, 3))));
    MatrB  =abs(fftshift(fft2(sum(fftMatrBt, 3))));
end

%----------------------------------------------------
function MatrA = make_fftcomb_bin(Freqs, Amps, Exp, kill_neg, lw)
% tMatrA, tMatrB  = [Exp.Dim, Exp.Dim, sim_size] (see main function)
% kill_neg = 1 or 0. If 1, negative part of the timedomain "spectra" will
%   be == 0
% make_prod amplitude of the combination freq. peaks. if 0, product rull
%   will be screwed
% lw  unitary linewidth IN POINTS, not in a MHZ !!!

ind1 = 1+floor( (Freqs(:, 1)+Exp.MaxFreq)/(abs(2*Exp.MaxFreq/(Exp.nPoints-1))) );
ind2 = 1+floor( (Freqs(:, 2)+Exp.MaxFreq)/(abs(2*Exp.MaxFreq/(Exp.nPoints-1))) );   

find1 = find(ind1(:)>0 & ind1(:) <= Exp.nPoints);
ind1 = ind1(find1);
ind2 = ind2(find1);   
find2 = find(ind2(:)>0 & ind2(:) <= Exp.nPoints);
MatrA = zeros(Exp.nPoints);
kill_zeros = 1;
if ~isempty(find2),
  ind1 = ind1(find2);
  ind2 = ind2(find2); 
  f_ind2 = find1(find2);
%        MatrA = bin2D(MatrA, ind1(:), ind2(:), (abs(real(Aikln(find2))) + abs(imag(Aikln(find2)))));
%        MatrB = bin2D(MatrB, ind2(:), ind1(:), (abs(real(Bikln(find2))) + abs(imag(Bikln(find2)))));   
   MatrA = bin2D(MatrA, ind1(:), ind2(:), (Amps(f_ind2, 1) + Amps(f_ind2, 2)));
%    MatrA = bin2D(MatrA, ind2(:), ind1(:), (abs(Amps(f_ind2, 1)) + abs(Amps(f_ind2, 2))));   
   MatrA = bin2D(MatrA, Exp.nPoints - ind1(:), Exp.nPoints - ind2(:), (Amps(f_ind2, 1) + Amps(f_ind2, 2)));
%    MatrA = bin2D(MatrA, Exp.nPoints - ind2(:), Exp.nPoints - ind1(:), (abs(Amps(f_ind2, 1)) + abs(Amps(f_ind2, 2))));   
end
MatrA = MatrA + MatrA.';
if kill_zeros
        gauss = gaussian(1:(Exp.nPoints), Exp.nPoints/2, Exp.nPoints/100);
        gauss = 1-gauss/max(gauss);
        Killing_mashine = gauss(ones(Exp.nPoints, 1), :).*gauss(ones(Exp.nPoints, 1), :).';
        MatrA = MatrA.*Killing_mashine;
end
if lw | kill_neg
    fftMatrAt = (ifft2(fftshift(MatrA/max(max(max(MatrA))))));
    % making linewidth using 2D gaussian
    if lw
        fftGaga = real(ifft2(fftshift(Local_gaga2(1:Exp.nPoints, Exp.nPoints/2+1, 1:Exp.nPoints, Exp.nPoints/2+1, [1, 1]*lw))));
        fftMatrAt = fftMatrAt.*fftGaga/max(max(fftGaga));
    end
    % killing revers time domains
    if kill_neg
        MS = floor(Exp.nPoints/2);
        fftMatrAt(MS:end, MS:end, :) = 0;
        fftMatrAt(1:MS, MS:end, :) = 0;
        fftMatrAt(MS:end, 1:MS,  :) = 0;
    end
    MatrA = abs(fftshift(fft2(fftMatrAt)));

end
% %----------------------------------------------------
% function [rMatrA, iMatrA] = make_bin(rMatrA, iMatrA, Freqs, Amps, Exp);
% % kill_neg = 1 or 0. If 1, negative part of the timedomain "spectra" will
% %   be == 0
% % make_prod amplitude of the combination freq. peaks. if 0, product rull
% %   will be screwed
% % lw  unitary linewidth IN POINTS, not in a MHZ !!!
% 
% nPoints = Exp.nPoints*2; 
% ind1 = 1+floor( (Freqs(:, 1)+Exp.MaxFreq)/(abs(2*Exp.MaxFreq/(Exp.nPoints-1))) );
% ind2 = 1+floor( (Freqs(:, 2)+Exp.MaxFreq)/(abs(2*Exp.MaxFreq/(Exp.nPoints-1))) );   
% 
% find1 = find(ind1(:)>0 & ind1(:) <= Exp.nPoints);
% ind1 = ind1(find1);
% ind2 = ind2(find1);   
% find2 = find(ind2(:)>0 & ind2(:) <= Exp.nPoints);
% 
% if ~isempty(find2),
%   ind1 = [ind1(find2), Exp.nPoints - ind1(find2)];
%   ind2 = [ind2(find2), Exp.nPoints - ind2(find2)]; %ind2(find2); 
%   f_ind2 = [1, 1]*find1(find2), find1(find2);
% 
%    rMatrA = bin2D(rMatrA, [ind1(:); ind2(:)], [ind2(:); ind1(:)], [real(Amps(f_ind2, 1)+ Amps(f_ind2, 2)); real(Amps(f_ind2, 1)+ Amps(f_ind2, 2))]);
% %    rMatrA = bin2D(rMatrA, ind2(:), ind1(:), real(Amps(f_ind2, 1)+ Amps(f_ind2, 2)));
%    
%    iMatrA = bin2D(iMatrA, [ind1(:); ind2(:)], [ind2(:); ind1(:)], [imag(Amps(f_ind2, 1)+ Amps(f_ind2, 2)); imag(Amps(f_ind2, 1)+ Amps(f_ind2, 2))]);
% %    iMatrA = bin2D(iMatrA, ind2(:), ind1(:), imag(Amps(f_ind2, 1)+ Amps(f_ind2, 2)));
%    
% %    MatrA = bin2D(MatrA, Exp.nPoints - ind1(:), Exp.nPoints - ind2(:), (Amps(f_ind2, 1) + Amps(f_ind2, 2)));
% %    MatrA = bin2D(MatrA, Exp.nPoints - ind2(:), Exp.nPoints - ind1(:), (Amps(f_ind2, 1) + Amps(f_ind2, 2)));
% 
% end
% % MatrA = sqrt(rMatrA.^2 + iMatrA.^2);
%----------------------------------------------------
function MatrA = make_fft(MatrA, Exp, lw)
% MatrA = abs(MatrA) + abs(MatrA.');

use_LoGa = 0;   %%%%%%%%%%% to make lorentz-gauss transofrmation / remove "stars"
[A, B] =size(MatrA) ;

xmean = mean(MatrA, 1);
MatrA = MatrA - xmean(ones(A, 1), :);

xmean = mean(MatrA, 2);
MatrA = MatrA - xmean(:, ones(A, 1));
if lw
%         fftGaga = abs(ifft2(fftshift(Local_gaga2(1:Exp.nPoints, Exp.nPoints/2+0.5, 1:Exp.nPoints, Exp.nPoints/2+0.5, [1, 1]*lw))));
        fftGaga = fftshift(Local_gaga2(1:Exp.nPoints, 0, 1:Exp.nPoints, 0, [1, 1]*Exp.nPoints*Exp.MaxFreq/lw));
        MatrA = MatrA.*fftGaga/max(max(fftGaga));
end
    % killing revers time domains
    
%     if kill_neg
% %         MS = floor(Exp.nPoints/2);
%         fftMatrAt(MS:end, MS:end, :) = 0;
%         fftMatrAt(1:MS, MS:end, :) = 0;
%         fftMatrAt(MS:end, 1:MS,  :) = 0;
%     end
    if use_LoGa
       gg =  Local_gaga2(1:Exp.nPoints, Exp.nPoints/4+0.5, 1:Exp.nPoints, Exp.nPoints/4+0.5, [1, 1]*Exp.nPoints/3);
       t = 1:Exp.nPoints;
       zzzff = t(ones(1, Exp.nPoints), :);
       fftMatrAt = fftMatrAt.*gg.*exp(zzzff/Exp.nPoints*3).*exp(zzzff.'/Exp.nPoints*3);
    end
    MatrA = abs(fftshift(fft2(MatrA)));
% end
%----------------------------------------------------
function [Ham_A, Ham_B] = Local_spinham_opt(Sys, Exp, HFIx, HFIy, HFIz, Quad)
% global sIxx sIxy sIxz sIyx sIyy sIyz sIzx sIzy sIzz sIx sIy sIz
% Generation of spin hamiltonian for 'ms'/'0' electron spin manifolds (S'=1/2)
% 
% [Ham_A, Ham_B] = spinham(Sys, Exp, ms)
% input units: mT, ns, MHz, rad
global sIx sIy sIz
Wplanck = planck;
% angl = length(Exp.phi);
Bo = Exp.Field*1e-3; % T
phi = Exp.phi;
theta = Exp.theta;

gxx = Sys.g(1);
gyy = Sys.g(2);
gzz = Sys.g(3);

mult = 2*Sys.I+1;

E = sop([Sys.I], 'e');

% if ~isfield(Sys, 'Apa'), Sys.Apa = [0, 0, 0]; end
% 
% A = diag(Sys.A)*1e6;
% RM = erot(Sys.Apa);
% A = RM*A*RM.';
% calcQ = 0;
% if isfield(Sys, 'Q')&~(Sys.I==1/2);
%     if sum(abs(Sys.Q)) ~= 0
% 
%         Q = diag(Sys.Q)*1e6;
%         if isfield(Sys, 'Qpa'), RMQ = erot(Sys.Qpa);
%         else RMQ = diag(ones(3, 1));
%         end
%         Q = RMQ*Q*RMQ.';
%         calcQ = 1;
%     end
% end
sn = 2*Sys.I + 1;
snum = ones(1, sn);

costheta = cos(theta);
sintheta = sin(theta);

cosphi = cos(phi);
sinphi = sin(phi);
% 
angle(1,1) = sintheta*cosphi;
angle(2,1) = sintheta*sinphi;
angle(3,1) = costheta;
sangle = sintheta;

Gx=angle(1,1)*Sys.g(1);
Gy=angle(2,1)*Sys.g(2);
Gz=angle(3,1)*Sys.g(3);

geff=sqrt(Gx.*Gx+Gy.*Gy+Gz.*Gz);

we = (1/Wplanck)*bmagn*Bo*geff;
nu_n = 1/Wplanck*nmagn*Sys.gn*Bo;

Zeeman_A = -nu_n*(angle(1,1)*sIx{mult} + angle(2,1)*sIy{mult} + angle(3,1)*sIz{mult});
%%%%%%333OLD
% HFI_A = HFIx*(gxx./geff)*sintheta*cosphi*1/2 + ...
%         HFIy*(gyy./geff)*sintheta*sinphi*1/2 + ...
%         HFIz*(gzz./geff)*costheta*1/2;
HFI_A = HFIx*Gx./geff*1/2 + ...
        HFIy*Gy./geff*1/2 + ...
        HFIz*Gz./geff*1/2;
Ham_A = Zeeman_A - HFI_A + Quad;
Ham_B = Zeeman_A + HFI_A + Quad;

%--------------------------------------------------------
function [Amp, W] = Local_crosspeak(Str, tau, ii, kk, ll, nn, jj, mm)
E_a = Str.E_a;
E_b = Str.E_b;
% sn = size(E_a, 3);
I = Str.I;
nst = 2*I+1;
V_a = Str.V_a;
V_b = Str.V_b;

%%%%%%%%% Frequencies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Hz]
wa = real(E_a(:, ones(size(E_a, 1), 1)) - E_a(:, ones(size(E_a, 1), 1)).');
wb = real(E_b(:, ones(size(E_b, 1), 1)) - E_b(:, ones(size(E_b, 1), 1)).');

%%%%%%%%% Transition moments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = V_a'*V_b; %%% corrected alsi 14.01.2008
cM = M';      %%% corrected alsi 20.01.2008

Amp = [];
W = [];
freqind = 1;
CC = 1;

%%%%%%%%% Calculation of crosspeaks frequencies and amplitudes
maxamp = 0;
nummax = 0;

for cc = 1:length(ii) % 1:nst^2*(cnt-1)^2 %
    Aikln = 0;
    Bikln = 0;
    for bb = 1:length(jj)
        ampar  = M(ii(cc), ll(cc)).*cM(ll(cc), jj(bb)).*M(jj(bb), nn(cc)).*cM(nn(cc), kk(cc)).*M(kk(cc), mm(bb)).*cM(mm(bb), ii(cc));
        exppar = exp(-2i*pi*(wa(ii(cc), jj(bb)) + wb(ll(cc), mm(bb)))*tau) + exp(2i*pi*(wa(kk(cc), jj(bb)) + wb(nn(cc), mm(bb)))*tau);
        Aikln =  (ampar.*exppar)+Aikln;
        Bikln =  (conj(ampar).*exppar)+Bikln;
    end
    Amp(end+1, [1, 2]) = ([Aikln, Bikln]);
      W(end+1, [1, 2]) = [wa(ii(cc), kk(cc)), wb(ll(cc), nn(cc))];
end

function varargout = Local_gaga2(x, xo, y, yo, fwhh, varargin)
% out = gaga(x, xo, y, yo, fwhh)
%  
%     y = Local_gaga2(x,x0, y, yo, fwhh)
%     Returns a 2D Gaussian line shape
%  
%     Input:
%     - 'x': First abscissa vector (for columns)
%     - 'x0': Centre of the lineshape function 
%     - 'y': Second abscissa vector (for rows)
%     - 'y0': Centre of the lineshape function
%     - 'fwhh': Full width at half height for 2 dimentions [fwhh_x fwhh_y]
%     - 'diff': Derivative. 0 is no derivative, 1 first,
%       2 second, -1 the integral with -infinity as lower
%       limit. If omitted, 0 is the default.
%  
%     Output:
%     - 'y': [n, m] vector  of function values for arguments 'x' (size=n)
%     and 'y' (size=m)

G = fwhh/sqrt(2*log(2));
sx = size(x);
sy = size(y);
if sx(1) ==1, xx(:,1) = x(:); else xx(:, 1) = x; end
if sy(1) ==1, yy(1, :) = y; else yy(1, :) = y.'; end
tx = xx(:, ones(length(yy), 1));
ty = yy(ones(length(xx), 1), :);
out = (2/pi)*(1/G(1))*(1/G(2))*exp(-2*((tx-xo)/G(1)).^2).*exp(-2*((ty-yo)/G(2)).^2);
if nargout > 1,
    varargout{1} = tx;
    varargout{2} = ty;
    varargout{3} = out/max(max(out));
elseif nargout == 1
    varargout{1} = out/max(max(out));
else
    disp('something is wrong');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% calculate orientation selection and produce a set of angles
function [phi, theta, weights] = Local_orisel(Sys, Exp, Opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%

if isfield(Exp, 'phi')& isfield(Exp, 'theta'),
    nphi = length(Exp.phi);
    ntheta = length(Exp.theta);

    theta = reshape(Exp.theta(ones(nphi, 1), :), nphi*ntheta, 1)*pi/180;
    phi = reshape(Exp.phi(ones(ntheta, 1), :).', nphi*ntheta, 1)*pi/180;
    ww = ones(nphi*ntheta, 1);
elseif isfield(Opt, 'Symmetry'),
    if ~Opt.Symmetry % if [0], use grid from MR2 (spiral [0:pi/2]-(0:2*pi))
        % taken from FORTRAN program of Ed Reijerse "MR2" 
        ncall = 0;
        nselect= 0;

        krid = safeget(Opt, 'nKnots', 20);

        theta = [];
        phi = [];
        last = 0;
        step = 0;
        while (~last)
            if(step==0 & ncall~=0), error(' PWDPEAL unitialized!'); end
            if(krid<0), error(' PWDPEAL illegal value of krid'); end
            if(ncall==0)
                itheta = krid;
                iphi   = 0;
                step = 90 / krid;
                nphi = 2*krid;
                thetaa = 90;
                dthe = -0.5 * step / nphi;
                dphi = 360 / nphi;
            end
            if(itheta<=0) error(' pdrpeal too many calls'); end
            if(iphi==nphi)
                iphi = 0;
                itheta = itheta -1;
                thetam = (itheta) * step;
                nphi = floor(sin(pi/180*thetam)*(4*krid));
                if(nphi<=0) error( 'BUG pdrpeal 1'); end
                % *              SMALLEST NPHI IS TYPICALLY 6
                thetaa = thetam + 0.5 * step;
                dthe = - step / (nphi);
                dphi = 360 / (nphi);
            end
            thetap = thetaa + (iphi) * dthe;
            phip   = (iphi) * dphi;
            iphi  = iphi + 1;
            ncall = ncall + 1;
            last = (iphi==nphi & itheta==1);
            if(last), step = 0; end

            theta(end+1, 1) = thetap*pi/180;
            phi(end+1, 1)   = phip*pi/180;
        end
        ww  = phi*0 + 1;
    else
        Opt.Symmetry = safeget(Opt, 'Symmetry', 'Ci');
        Opt.nKnots = safeget(Opt, 'nKnots', 20);
        [phi, theta, ww] = sphgrid(Opt.Symmetry, Opt.nKnots);
        phi = phi.'; theta = theta.'; ww = ww.';
    end

else
    t_theta = linspace(0, pi, safeget(Opt, 'nKnots', 20));
    Opt.nKnots = safeget(Opt, 'nKnots', 20);
    phi = [];
    theta = [];
    nnn = [];
    for cc = 1:Opt.nKnots
        num = floor(Opt.nKnots*abs(sin(t_theta(cc))));
        nnn(cc) = num;
        if num ~=0, 
            t_phi = linspace(0, 2*pi, num+1).';
            phi = [phi; t_phi(1:(end-1), 1)];
            theta = [theta; ones(num, 1)*t_theta(cc)];
        end
    end
    ww  = phi*0 + 1;

end

nangles = length(phi);
if isfield(Sys, 'Nucs'), [Sys.I,Sys.gn] = nucdata(Sys.Nucs); end % EasySpin function



% calculating orientation selection
ctheta = cos(theta);
stheta = sin(theta);
cphi = cos(phi);
sphi = sin(phi);
angle(:, 1)=stheta.*cphi;
angle(:, 2)=stheta.*sphi;
angle(:, 3) = ctheta;

Gx=angle(:, 1)*Sys.g(1);
Gy=angle(:, 2)*Sys.g(2);
Gz=angle(:, 3)*Sys.g(3);
geff=sqrt(Gx.*Gx+Gy.*Gy+Gz.*Gz);

%%%%%%%%%%%%%%%%%%%% new %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% include HF to the orientation selection 
useHFsel = safeget(Opt, 'useHFsel', 0);
if useHFsel
    projA = zeros(nangles, sum(Sys.nNucs));
    k = 0;
    for ii = 1:length(Sys.nNucs)
        for jj = 1:Sys.nNucs(ii)
            k = k+1;
            if sum(Sys.Apa(ii, :))
                R = erot(Sys.Apa(ii, :)*pi/180);
                At = R*diag(Sys.A(ii, :))*R';
            else
                At = diag(Sys.A(ii, :));
            end
            LA(:, 1) = angle(:, 1) * At(1, 1) + angle(:, 2) * At(2, 1) + angle(:, 3) * At(3, 1);
            LA(:, 2) = angle(:, 1) * At(1, 2) + angle(:, 2) * At(2, 2) + angle(:, 3) * At(3, 2);
            LA(:, 3) = angle(:, 1) * At(1, 3) + angle(:, 2) * At(2, 3) + angle(:, 3) * At(3, 3);
            projA(:, k) = sqrt(LA(:, 1).^2 + LA(:, 2).^2 + LA(:, 3).^2);
            mult(k) = 2*Sys.I(ii)+1;
            
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Exp.ExciteWidth = safeget(Exp, 'ExciteWidth', 0);
Exp.mwFreq = safeget(Exp, 'mwFreq', 0);
if isfield(Exp, 'gField'),
   if ~Exp.mwFreq & ~isfield(Exp, 'Field'), ...
           error('ESEEM: there is no Exp.mwFreq to calculate B0');
   end
   Exp.Field = (Exp.mwFreq*1e+9*planck)/(Exp.gField*bmagn*1e-3);  % mT
else 
    Exp.Field = safeget(Exp, 'Field', 0);
end

if (Exp.ExciteWidth > 0) & (Exp.mwFreq>0)
%     gcenter = fld2g(Exp.Field*1E-3, Exp.mwFreq*1E9); % killed by alsi 27.12.07

    feff = Exp.Field*1E-3*bmagn*geff/planck *1e-6; % Must be in MHZ
    if ~useHFsel
        
        ak = exp(-2*((feff-Exp.mwFreq*1e3)./Exp.ExciteWidth).^2); % all in MHz
    else
        % dE (epr) = nu_S + Sum(a*MIi)
        nTrans = prod(mult);
        if length(mult)==1,
            MIs(:, 1) = (-(mult(1)-1)/2 : 1 : (mult(1)-1)/2)';
        else
            prod2 = prod(mult(2:end));
            tA = (-(mult(1)-1)/2 : 1 : (mult(1)-1)/2)';
            MIs(:, 1) = kron(tA, ones(prod2, 1));

            for ii = 2:(length(mult)-1)
                tA = (-(mult(ii)-1)/2 : 1 : (mult(ii)-1)/2)';
                prod1 = prod(mult(1:(ii-1)));
                ttA = kron(ones(prod1, 1), tA);
                prod2 = prod(mult((ii+1):end));
                MIs(:, ii) = kron(ttA, ones(prod2, 1));
            end
            tA = (-(mult(end)-1)/2 : 1 : (mult(end)-1)/2)';
            prod2 = prod(mult(1:(end-1)));
            MIs(:, length(mult)) = kron(ones(prod2, 1), tA);
       end
        ak = zeros(nangles, 1);
        for ii = 1:nTrans
            ffre = zeros(nangles, 1);
            for jj = 1:length(mult)
                ffre = ffre + projA(:, jj)*MIs(ii, jj);
            end
            ak = ak + exp(-2*((feff+ffre-Exp.mwFreq*1e3)./Exp.ExciteWidth).^2); % all in MHz
        end
        
    end

else
    ak = ones(nangles, 1);
end
%         Sy.Nucs = '57Fe';
%         Parameters ={'A', 'Apa', 'Q', 'Qpa', 'Nucs', 'gn', 'I'};
%         Sy = Sys;
%         for kk = 1:length(isfield_num)
%             pk = isfield_num(kk);
%             tPar = getfield(Sys, Parameters{pk});
%             rmfield(Sy, Parameters{pk})
%             Sy = setfield(Sy, Parameters{pk}, tPar(1, :));
%         end
%         ak = orisel(Sy, Exp, [phi.', theta.']);    
if max(ak) < 1e-4, error('_orisel: no resonance'); end 
weights = ak.*ww;
weights = weights/max(weights);
% Weights = orisel(Sys,Par,Ori)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% generate spin operators %%%%%%%%%%%%%%%%%%%%%%%%%
function totalout = local_sop(spins, corrs)
% totalout = local_sop(spins, corrs)
%    spins - row of total spins (arbitrary number)
%    corrs - string containing 'x', 'y' , 'z' or 'e';
totalout = 1;

for jj=1:length(spins)
    s = spins(jj);
    mult = (2*s+1);
    sp = zeros(mult);
    sm = zeros(mult);
    %%%%%%%%%%%%%% calculate Sx
    if strcmp(corrs(jj), 'x')
        for ii =1:mult-1
            jj = ii+1;
            ms = -s+(ii-1);
            sp(ii, jj) = sqrt(s*(s+1) - ms*(ms+1));
        end

        for ii =2:mult
            jj = ii-1;
            ms = -s+(ii-1);
            sm(ii, jj) = sqrt(s*(s+1) - ms*(ms-1));
        end

        S_s = (sp+sm)/2;
    %%%%%%%%%%%%% calculate Sy
    elseif strcmp(corrs(jj), 'y')
        for ii =1:mult-1
            jj = ii+1;
            ms = -s+(ii-1);
            sp(ii, jj) = sqrt(s*(s+1) - ms*(ms+1));
        end

        for ii =2:mult
            jj = ii-1;
            ms = -s+(ii-1);
            sm(ii, jj) = sqrt(s*(s+1) - ms*(ms-1));
        end

        S_s = (sp-sm)/2i;
    %%%%%%%%%%%%%% calculate Sz
    elseif strcmp(corrs(jj), 'z')
        S_s = diag([-s:1:s]);
    %%%%%%%%%%%%%% if it non of those, it is a unitary matrix    
    else % 'e'
        S_s = diag(ones(mult, 1));
    end
    totalout = kron(totalout, S_s);
end

