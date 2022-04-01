function [y, ax] = kv_hyscorefd(Sys, Exp, Opt, Shift)
% input:
% Sys.A/Q = [MHz], Exp.mwFreq = [GHz], Exp.Field = [mT], Exp.tau = [ns]
% Sys.Apa/Qpa = [degree]
tstart = cputime;
tau = safeget(Exp, 'tau', 0)*1e-9; %ns->s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Exp, 'phi')& isfield(Exp, 'theta'),
    nphi = length(Exp.phi);
    ntheta = length(Exp.theta);

    theta = reshape(Exp.theta(ones(nphi, 1), :), nphi*ntheta, 1)*pi/180;
    phi = reshape(Exp.phi(ones(ntheta, 1), :).', nphi*ntheta, 1)*pi/180;
    ww = ones(nphi*ntheta, 1);
elseif isfield(Opt, 'Symmetry'),
    Opt.Symmetry = safeget(Opt, 'Symmetry', 'Ci');
    Opt.nKnots = safeget(Opt, 'nKnots', 20);
    [phi, theta, ww] = sphgrid(Opt.Symmetry, Opt.nKnots);
    phi = phi.'; theta = theta.'; ww = ww.';
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

Exp.ExciteWidth = safeget(Exp, 'ExciteWidth', 0);
Exp.mwFreq = safeget(Exp, 'mwFreq', 0);
if isfield(Exp, 'gField'),
   if ~Exp.mwFreq & ~isfield(Exp, 'Field'), ...
           error('hycalc: there is no Exp.mwFreq to calculate B0');
   end
   Exp.Field = (Exp.mwFreq*1e+9*planck)/(Exp.gField*bmagn*1e-3);  % mT
else 
    Exp.Field = safeget(Exp, 'Field', 0);
end

if (Exp.ExciteWidth > 0) & (Exp.mwFreq>0)
    gcenter = fld2g(Exp.Field*1E-3, Exp.mwFreq*1E9);
    gw      = gcenter*Exp.ExciteWidth*1e6*(planck)/(2*bmagn*1e-3)/Exp.Field /sqrt(2*log(2));
    ak = exp(-2*((geff-gcenter)./gw).^2);
else
    ak = ones(ntheta*nphi, 1);
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
    
if max(ak) < 1e-3, disp('no resonance'); disperr(handles); return; end 
% Weights = orisel(Sys,Par,Ori)

if isfield(Sys, 'Nucs'), [Sys.I,Sys.gn] = nucdata(Sys.Nucs); end

Exp.tau = safeget(Exp, 'tau', 0);
Opt.Treshold = safeget(Opt, 'Treshold', 1e-3);

oriind = find(ak >1e-4);
Exp.phi = phi(oriind, 1);
Exp.theta = theta(oriind, 1);
Exp.Weights = ww(oriind, 1).*ak(oriind,1 )/max(ww(oriind, 1));

%  test sphere for orientation selection (doted sphere)
if safeget(Opt, 'ShowOri', 0)
    figure(125);
    [xx, yy, zz] = sphere(20);
    surf(xx, yy, zz, 'Linestyle', 'none'); %zz*0)
    colormap([0.5, 0.5, 1])
    shading interp;
    hold on;
    X =ang2vec(phi(oriind, 1), theta(oriind, 1));
    www = ww(oriind, 1).*ak(oriind,1 )/max(ww(:, 1));
    for ii = 1:length(X(1, :))
        plot3(X(1, ii), X(2, ii), X(3, ii), 'Marker', '.', 'LineStyle', 'none', 'Color', [www(ii), 0, 0]);
        hold on;
    end
  %  X1 =ang2vec(Exp.phi, Exp.theta);
  %  plot3(X1(1, :), X1(2, :), X1(3, :), '.r');
    hold off;
    
%    figure(124);
%     X =ang2vec(phi(:, 1), theta(:, 1));
%     tres = delaunay(X(:, 1), X(:, 2);
%     tres = voronoin(X);    
%     trimesh(tres, X(:, 1), X(:, 2), X(:, 1));
%         www = ww(:, 1).*ak(:,1 )/max(ww(:, 1));
%     for ii = 1:length(X(1, :))
%         plot3(X(1, ii), X(2, ii), X(3, ii), 'Marker', '*', 'LineStyle', 'none', 'Color', [www(ii), 0, 0]);
%         hold on;
%     end
  %  X1 =ang2vec(Exp.phi, Exp.theta);
  %  plot3(X1(1, :), X1(2, :), X1(3, :), '.r');
%     hold off;
    %figure(123); plot3(phi*180/pi, theta*180/pi, ak, '.')
    % xlabel('phi'); ylabel('theta');
end

% simulation for few nuclei
Parameters ={'A', 'Apa', 'Q', 'Qpa', 'Nucs', 'gn', 'I'};
% lets find number of nuclei (parameters can be common, it means size(???, 1) can be 1)
Sys.Apa = safeget(Sys, 'Apa', [0, 0, 0])*pi/180; % degree -> rad
Sys.Qpa = safeget(Sys, 'Qpa', [0, 0, 0])*pi/180; % degree -> rad

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
    for kk = 1:length(isfield_num)
        pk = isfield_num(kk); % kind of optimisation
        tPar = getfield(Sys, Parameters{pk});
        if size(tPar, 1)==1
            Sys = setfield(Sys, Parameters{pk}, tPar(ones(sim_size, 1), :));
        elseif size(tPar, 1)~=sim_size
            error({['Parameter Sys.', Parameters{pk}, ' has wrong number of rows.']; ['Correct number is ', num2str(sim_size)]});
        end
    end
end

if isfield(Sys, 'Aiso'),
    Sys.A = safeget(Sys, 'A', 0) + Sys.Aiso;
end

MatrA = zeros(safeget(Exp, 'nPoints', 500));
MatrB = MatrA;
% start calculation
AmpRation = safeget(Opt, 'AmpRatios', ones(sim_size, 1));
if (length(AmpRation)~=sim_size), error('length(AmpRatios)~=size(A, 1)'); end
for ci = 1:sim_size
    Sy = Sys;
    for kk = 1:length(isfield_num)
        pk = isfield_num(kk);
        tPar = getfield(Sy, Parameters{pk});
        Sy = setfield(Sy, Parameters{pk}, tPar(ci, :));
    end
    [HamA, HamB] = Local_spinham(Sy, Exp); 
    k = 0;
    for ang = 1:length(Exp.phi)
        [V_a(:, :, ang), tD_a] = eig(HamA(:, :, ang));
        D_a(:, 1, ang) = diag(tD_a);
        [V_b(:, :, ang), tD_b] = eig(HamB(:, :, ang));
        D_b(:, 1, ang) = diag(tD_b);
        if sum(imag(D_b(:, 1, ang)))>1e-9|sum(imag(D_a(:, 1, ang)))>1e-9
            k = max(k, max(sum(imag(D_b(:, 1, ang))),sum(imag(D_a(:, 1, ang))))); 
        end
    end
    if k % imag ~=0 means problems with hamiltonian, check Local_spinham
%         error(['Bad Hamiltonian, max(imag(D)):', num2str(k)]);
        disp(['Bad Hamiltonian, max(imag(D)):', num2str(k)]);
    end
    Str.V_a = V_a;
    Str.V_b = V_b;
    Str.E_a = D_a;
    Str.E_b = D_b;
    Str.vi = nmagn*Sy.gn*Exp.Field/planck;
    Str.I = Sy.I;
    Str.Weights = Exp.Weights;
    Str.Dim = safeget(Exp, 'nPoints', 500);
    Str.Board = safeget(Exp, 'MaxFreq', 10)*1e6;
    
    [tMatrA(:, :, ci), tMatrB(:, :, ci)] = Local_crosspeack_fd(Str, tau);
    if AmpRation(ci,1)~=1
        tMatrA(:, :, ci) = tMatrA(:, :, ci)*AmpRation(ci,1);
        tMatrB(:, :, ci) = tMatrB(:, :, ci)*AmpRation(ci,1);
    end
end
make_prod = safeget(Opt, 'ProdRule', 1);
kill_neg = safeget(Opt, 'KillNeg', 0);
dx = Str.Board*2/Str.Dim;
if sim_size>1
    [MatrA, MatrB] = make_fftcomb(tMatrA, tMatrB, kill_neg, make_prod, safeget(Sys, 'lw', 0)*1e6/dx);
else
    [MatrA, MatrB] = make_fftcomb(tMatrA, tMatrB, kill_neg, 0, safeget(Sys, 'lw', 0)*1e6/dx);
end
ax.x = linspace(-Str.Board, Str.Board, Str.Dim).';
ax.y = linspace(-Str.Board, Str.Board, Str.Dim).';

y = MatrA+MatrB; 

if safeget(Shift, 'ding', 0)
    x = [1:1000]/1.5;
    sound(min(Shift.ding, 1)*sin(x).*sin(x/2));
end
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

%----------------------------------------------------
function [MatrA, MatrB] = make_fftcomb(tMatrA, tMatrB, kill_neg, make_prod, lw)
% tMatrA, tMatrB  = [Exp.Dim, Exp.Dim, sim_size] (see main function)
% kill_neg = 1 or 0. If kill_neg ==1, negative part of the timedomain
% "spectra" will be == 0
% make_prod - amplitude of the combination freq. peaks. if 0, product rull
%   will be screwed
% lw  unitary linewidth IN POINTS, not in a MHZ (internal parameter)!!!

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
function [Ham_A, Ham_B] = Local_spinham(Sys, Exp)
% Generating spin hamiltonian for alpha/beta electron spin manifolds (S'=1/2)
% [Ham_A, Ham_B] = spinham(Sys, Exp)
% input units: mT, ns, MHz, rad
Wplanck = planck;
angl = length(Exp.phi);

Bo = Exp.Field*1e-3; % T
phi(1, 1, :) = Exp.phi;
theta(1, 1, :) = Exp.theta;

gxx = Sys.g(1);
gyy = Sys.g(2);
gzz = Sys.g(3);
Ix = sop([Sys.I], 'x');
Iy = sop([Sys.I], 'Y');
Iz = sop([Sys.I], 'z');

Ixx = sop([Sys.I], 'x')*sop([Sys.I], 'x');
Ixy = sop([Sys.I], 'x')*sop([Sys.I], 'y');
Ixz = sop([Sys.I], 'x')*sop([Sys.I], 'z');

Iyx = sop([Sys.I], 'y')*sop([Sys.I], 'x');
Iyy = sop([Sys.I], 'Y')*sop([Sys.I], 'y');
Iyz = sop([Sys.I], 'y')*sop([Sys.I], 'z');

Izx = sop([Sys.I], 'z')*sop([Sys.I], 'x');
Izy = sop([Sys.I], 'z')*sop([Sys.I], 'y');
Izz = sop([Sys.I], 'z')*sop([Sys.I], 'z');

E = sop([Sys.I], 'e');

if ~isfield(Sys, 'Apa'), Sys.Apa = [0, 0, 0]; end

A = diag(Sys.A)*1e6;
RM = erot(Sys.Apa);
A = RM*A*RM.';
calcQ = 0;
if isfield(Sys, 'Q')&~(Sys.I==1/2);
    Q = diag(Sys.Q)*1e6;
    if isfield(Sys, 'Qpa'), RMQ = erot(Sys.Qpa);
    else RMQ = diag(ones(3, 1)); 
    end
    Q = RMQ*Q*RMQ.';
    calcQ = 1;
end
sn = 2*Sys.I + 1;
snum = ones(1, sn);
angnum = ones(angl, 1);

geff = sqrt(gzz^2*cos(theta(snum, snum, :)).^2 + gxx^2*(sin(theta(snum, snum, :)).*cos(phi(snum, snum, :))).^2 ...
    + gyy^2*(sin(theta(snum, snum, :)).*sin(phi(snum, snum, :))).^2); 

we = (1/Wplanck)*bmagn*Bo*geff;
nu_n = 1/Wplanck*nmagn*Sys.gn*Bo;

Zeeman_A = -nu_n*(sin(theta(snum, snum, :)).*cos(phi(snum, snum, :)).*Ix(:, :, angnum)...
    + sin(theta(snum, snum, :)).*sin(phi(snum, snum, :)).*Iy(:, :, angnum)...
    + cos(theta(snum, snum, :)).*Iz(:, :, angnum));

HFI_A = (Ix(:, :, angnum)*A(1, 1) + Iy(:, :, angnum)*A(1, 2) + Iz(:, :, angnum)*A(1, 3)).*...
            (gxx./geff).*sin(theta(snum, snum, :)).*cos(phi(snum, snum, :))*1/2 + ...
        (Ix(:, :, angnum)*A(2, 1) + Iy(:, :, angnum)*A(2, 2) + Iz(:, :, angnum)*A(2, 3)).*...
            (gyy./geff).*sin(theta(snum, snum, :)).*sin(phi(snum, snum, :))*1/2 + ...
        (Ix(:, :, angnum)*A(3, 1) + Iy(:, :, angnum)*A(3, 2) + Iz(:, :, angnum)*A(3, 3)).*...
            (gzz./geff).*cos(theta(snum, snum, :))*1/2;

if calcQ
    Quad = Q(1,1)*Ixx(:, :, angnum) + Q(1, 2)*Ixy(:, :, angnum)+ Q(1, 3)*Ixz(:, :, angnum)+...
           Q(2,1)*Iyx(:, :, angnum) + Q(2, 2)*Iyy(:, :, angnum)+ Q(2, 3)*Iyz(:, :, angnum)+...
           Q(3,1)*Izx(:, :, angnum) + Q(3, 2)*Izy(:, :, angnum)+ Q(3, 3)*Izz(:, :, angnum);
    %disp(Quad);
else
    Quad = 0;
end

Ham_A = Zeeman_A + HFI_A + Quad;
Ham_B = Zeeman_A - HFI_A + Quad;

%--------------------------------------------------------
function [MatrA, MatrB] = Local_crosspeack_fd(Str, tau)
E_a = Str.E_a;
E_b = Str.E_b;
sn = size(E_a, 3);
I = Str.I;
nst = 2*I+1;
V_a = Str.V_a;
V_b = Str.V_b;

% generating all possible combination of 6 numbers numbers in range [1:nst]
idxrange = 1:nst;
b1(:, 1, 1, 1, 1, 1) = idxrange;
b2(1, :, 1, 1, 1, 1) = idxrange;
b3(1, 1, :, 1, 1, 1) = idxrange;
b4(1, 1, 1, :, 1, 1) = idxrange;
b5(1, 1, 1, 1, :, 1) = idxrange;
b6(1, 1, 1, 1, 1, :) = idxrange;

nstones = ones(nst, 1);
ii = reshape(b1(:, nstones, nstones, nstones, nstones, nstones), nst^6, 1);
kk = reshape(b2(nstones, :, nstones, nstones, nstones, nstones), nst^6, 1);
ll = reshape(b3(nstones, nstones, :, nstones, nstones, nstones), nst^6, 1);
nn = reshape(b4(nstones, nstones, nstones, :, nstones, nstones), nst^6, 1);
jj = reshape(b5(nstones, nstones, nstones, nstones, :, nstones), nst^6, 1);
mm = reshape(b6(nstones, nstones, nstones, nstones, nstones, :), nst^6, 1);

% other way to deal with this
% ii = [];
% kk = [];
% ll = [];
% nn = [];
% jj = [];
% mm = [];
% for cii = 1:(nst-1)
%     for ckk = (cii+1):nst
%         for cll = 1:(nst-1)
%             for cnn = (cll+1):nst
%                 for cjj = 1:nst
%                     for cmm = 1:nst
%                         ii(end+1) = cii;
%                         kk(end+1) = ckk;
%                         ll(end+1) = cll;
%                         nn(end+1) = cnn;
%                         jj(end+1) = cjj;
%                         mm(end+1) = cmm;                        
%                     end
%                 end
%             end
%         end
%     end
% end
%%%%%%%%% Frequencies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hz, Original construction:  [i, j, Nangles] 
wa = real(E_a(:, ones(size(E_a, 1), 1), :) - permute(E_a(:, ones(size(E_a, 1), 1), :), [2, 1, 3]));
wb = real(E_b(:, ones(size(E_b, 1), 1), :) - permute(E_b(:, ones(size(E_b, 1), 1), :), [2, 1, 3]));

%%%%%%%%% Transition moments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for aass = 1:sn
   M(:, :, aass) = V_a(:, :, aass)'*V_b(:, :, aass);
end
cM = conj(M);

Amp = [];
W = [];
freqind = 1;
CC = 1;

MatrA = zeros(Str.Dim, Str.Dim);
MatrB = MatrA;
ax1 = linspace(-Str.Board, Str.Board, Str.Dim);
AMax = 0;
BMax = 0;
%%%%%%%%% Calculation of crosspeaks frequencies and amplitudes
maxamp = 0;
nummax = 0;
for cc = 1:length(ii) % 1:nst^2*(cnt-1)^2 % 
     %if (ii(cc)== kk(cc))&(ll(cc)~= nn(cc))|(ii(cc)~= kk(cc))&(ll(cc)== nn(cc)), 
    if (ii(cc)== kk(cc))|(ll(cc)== nn(cc)),      % to kill w_ik=0 or w_ln=0
        continue;
    end
   % amppar = Mil*Mlj'*Mjn'*Mnk'*Mkm'*Mmi'
    ampar  = M(ii(cc), ll(cc), :).*cM(ll(cc), jj(cc), :).*M(jj(cc), nn(cc), :).*cM(nn(cc), kk(cc), :).*M(kk(cc), mm(cc), :).*cM(mm(cc), ii(cc), :);
    
%   if abs(sum(ampar)/sn)<1e-35, continue; end
%        if (ii(cc)~= kk(cc))&(ll(cc)~= nn(cc)), 
%             disp('Yes');
%         end

   CC = CC+1;
   % exppar = exp(-i*pi*(waij + wblm)*tau) + exp(i*pi*(wakj + wbnm)*tau)
   exppar = exp(-i*2*pi*(wa(ii(cc), jj(cc), :) + wb(ll(cc), mm(cc), :))*tau) + exp(i*2*pi*(wa(kk(cc), jj(cc), :) + wb(nn(cc), mm(cc), :))*tau);

   % Weights[Nphi] -> Weights[1, 1, Nphi], because ampar and exppar are [1, 1, Nphi]
   Aikln(:, 1) =  ampar.*exppar;
   Bikln(:, 1) =  conj(ampar).*exppar;
   %Amp(end+[1:sn], [1, 2]) = abs([real(Aikln) + imag(Aikln), real(Bikln) + imag(Bikln)]);
%     Amp(end+[1:sn], [1, 2]) = abs([Aikln, Bikln]);
%     W(end+[1:sn], [1, 2]) = [permute(wa(ii(cc), kk(cc), :), [3, 1, 2]), permute(wb(ll(cc), nn(cc), :), [3, 1, 2])];
   % "2+" is a correction of the mistake in the binning (~1 point shift)
   ind1 = floor(1+(wa(ii(cc), kk(cc), :)-ax1(1))/(abs(ax1(1)-ax1(2))));
   ind2 = floor(1+(wb(ll(cc), nn(cc), :)-ax1(1))/(abs(ax1(1)-ax1(2))));   
  
   find1 = find(ind1(:)>0 & ind1(:) <= Str.Dim);
   ind1 = ind1(find1);
   ind2 = ind2(find1);   
   find2 = find(ind2(:)>0 & ind2(:) <= Str.Dim);
    
%    for ss = 1:size(find2)
% %        if ind2(find1(ss))>0 & ind2(find1(ss))<=Str.Dim
% %           MatrA(ind1(1, 1, find1(ss)), ind2(1, 1, find1(ss))) = abs(real(Aikln(find1(ss))) + imag(Aikln(find1(ss))));
% %            MatrB(ind2(1, 1, find1(ss)), ind1(1, 1, find1(ss))) = abs(real(Bikln(find1(ss))) + imag(Bikln(find1(ss))));
%             MatrA(ind1(1, 1, find2(ss)), ind2(1, 1, find2(ss))) = (abs(real(Aikln(find2(ss)))) + abs(imag(Aikln(find2(ss)))))*Str.Weights(find2(ss), 1);
%             MatrB(ind2(1, 1, find2(ss)), ind1(1, 1, find2(ss))) = (abs(real(Bikln(find2(ss)))) + abs(imag(Bikln(find2(ss)))))*Str.Weights(find2(ss), 1);
%             if abs(real(Aikln(find2(ss)))) + abs(imag(Aikln(find2(ss))))>maxamp, maxamp = abs(real(Aikln(find2(ss)))) + abs(imag(Aikln(find2(ss))));
%                 nummax(1) = wa(ii(cc), kk(cc)); nummax(2) = wb(ll(cc), nn(cc)); nummax(3) = ss; disp(cc); end
%             if abs(real(Bikln(find2(ss)))) + abs(imag(Bikln(find2(ss))))>maxamp, maxamp = abs(real(Bikln(find2(ss)))) + abs(imag(Bikln(find2(ss)))); 
%                 nummax(1) = wa(ii(cc), kk(cc)); nummax(2) = wb(ll(cc), nn(cc)); nummax(3) = ss; end
% %        end
%    end
    if ~isempty(find2),
       ind1 = ind1(find2);
       ind2 = ind2(find2); 
       f_ind2 = find1(find2);
%        MatrA = bin2D(MatrA, ind1(:), ind2(:), (abs(real(Aikln(find2))) + abs(imag(Aikln(find2)))));
%        MatrB = bin2D(MatrB, ind2(:), ind1(:), (abs(real(Bikln(find2))) + abs(imag(Bikln(find2)))));   
       MatrA = bin2D(MatrA, ind1(:), ind2(:), (abs(Aikln(f_ind2)) + abs(Bikln(f_ind2))).*Str.Weights(f_ind2, 1));
       MatrB = bin2D(MatrB, ind2(:), ind1(:), (abs(Aikln(f_ind2)) + abs(Bikln(f_ind2))).*Str.Weights(f_ind2, 1));   
       MatrA = bin2D(MatrA, Str.Dim - ind1(:), Str.Dim - ind2(:), (abs(Aikln(f_ind2)) + abs(Bikln(f_ind2))).*Str.Weights(f_ind2, 1));
       MatrB = bin2D(MatrB, Str.Dim - ind2(:), Str.Dim - ind1(:), (abs(Aikln(f_ind2)) + abs(Bikln(f_ind2))).*Str.Weights(f_ind2, 1));   
   end
   clear Aikln Bikln;
end
   
function varargout = Local_gaga2(x, xo, y, yo, fwhh, varargin)
% out = gaga(x, xo, y, yo, fwhh)
%  
%     y = gaussian(x,x0, y, yo, fwhh)
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
