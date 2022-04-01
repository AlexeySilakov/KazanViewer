% Simple and fast routine for calculating ENDOR spectra
% of one I=1 nucleus coupled to S=1/2 electron.
% 
% [x,y] = sugar1(Sys,Exp,Opt)
% 
% Partially compatible by parameters with EasySpin by Stefan Stoll.
% Parameters list: 
%  Sys: g,gn,A,Apa,Q,Qpa,lwEndor
%  Exp: Field, Range, mwFreq
%  Opt: ThetaRange,PhiRange,nKnots,nPoints
% External dependencies:  safeget, binning

% Boris Epel, MPI for Bioinorganic Chemistry, 2002-04

%==========================================================================
function [Freqs, Spec] = sugar1b(Sys,Ex,Si)
%==========================================================================
planck = 6.6261e-034;
nmagn = 5.0508e-027;
pi = 3.141593;

Parameters ={'A', 'Apa', 'Q', 'Qpa', 'Nucs', 'gn', 'I', 'nNucs'};

% Sys.Apa = safeget(Sys, 'Apa', [0, 0, 0])*pi/180; % degree -> rad
% Sys.Qpa = safeget(Sys, 'Qpa', [0, 0, 0])*pi/180; % degree -> rad
% lets find number of simmulaitons (parameters can be common, it means
% size(XX, 1)==1)
if isfield(Sys, 'Nucs'), [Sys.I,Sys.gn] = nucdata(Sys.Nucs); end % EasySpin function
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

[phi,theta,sangle] = sphgrid('Ci',Si.nKnots);

ctheta = cos(theta);
stheta = sin(theta);

cphi = cos(phi);
sphi = sin(phi);

angle(1,:)= stheta(:).*cphi(:);
angle(2,:)= stheta(:).*sphi(:);
angle(3,:)= ctheta(:);

% angle(1,:)=reshape(stheta(:)*cphi(:).', 1, nphi*ntheta);
% angle(2,:)=reshape(stheta(:)*sphi(:).', 1, nphi*ntheta);
% angle(3,:) = reshape(ctheta(:)*ones(1,nphi), 1, nphi*ntheta);

% sangle = reshape(stheta(:)*ones(1,nphi), 1, nphi*ntheta);

Gx=angle(1,:)*Sys.g(1);
Gy=angle(2,:)*Sys.g(2);
Gz=angle(3,:)*Sys.g(3);

geff=sqrt(Gx.*Gx+Gy.*Gy+Gz.*Gz);

Spec = zeros(1,Ex.nPoints);


for nucs = 1:length(Sys.gn)
    Sy = Sys;
    Sy.A = Sys.A(nucs, :);
    Sy.Q = Sys.Q(nucs, :);
    Sy.Apa = Sys.Apa(nucs, :);
    Sy.Qpa = Sys.Qpa(nucs, :);
    Sy.I = Sys.I(nucs);
    Sy.gn = Sys.gn(nucs);

%     if Sy.I==1/2, disp('Warning: for case Sys.I=1/2 ''sugar'' may be used.');
%     elseif Sy.I>1
%         error('case Sys.I>1 is not supported.');
%     end
    
    R = erot(safeget(Sy, 'Apa', [0,0,0]));
    A = R * diag(Sy.A)*R.';
    Ix =local_sop(Sy.I, 'x');
    Iy =local_sop(Sy.I, 'y');
    Iz =local_sop(Sy.I, 'z');
    
    iix =local_sop(Sy.I, 'x');
    iiy =local_sop(Sy.I, 'y');
    iiz =local_sop(Sy.I, 'z');
      
    omega_l = Ex.Field/planck /1E9 * Sy.gn * nmagn;
    A = diag(Sy.A); % [MHz]
    RM = erot(Sy.Apa);
    A = RM*A*RM.';

    HFIx = (Ix*A(1, 1) + Iy*A(2, 1) + Iz*A(3, 1));
    HFIy = (Ix*A(1, 2) + Iy*A(2, 2) + Iz*A(3, 2));
    HFIz = (Ix*A(1, 3) + Iy*A(2, 3) + Iz*A(3, 3));
    
    
    hfx = (iix*A(1, 1) + iiy*A(2, 1) + iiz*A(3, 1));
    hfy = (iix*A(1, 2) + iiy*A(2, 2) + iiz*A(3, 2));
    hfz = (iix*A(1, 3) + iiy*A(2, 3) + iiz*A(3, 3));
    
%       Quadrupole coupling 
    Quad = 0;
    if isfield(Sy, 'Q')&(Sy.I>1/2);
        if sum(abs(Sy.Q)) ~= 0

            Q = diag(Sy.Q); %[MHz]
            if isfield(Sy, 'Qpa'),
                RMQ = erot(Sy.Qpa);
            else
                RMQ = diag(ones(3, 1));
            end
            Q = RMQ*Q*RMQ.';
            Quad = Q(1,1)*Ix*Ix + Q(1, 2)*Ix*Iy + Q(1, 3)*Ix*Iz +...
                       Q(2,1)*Iy*Ix + Q(2, 2)*Iy*Iy + Q(2, 3)*Iy*Iz +...
                       Q(3,1)*Iz*Ix + Q(3, 2)*Iz*Iy + Q(3, 3)*Iz*Iz;
        end
    end
    mult = 2*Sy.I+1;
    onn = ones(mult, 1);    
    if Ex.ExciteWidth > 0
        gcenter = fld2g(Ex.Field*1E-3, Ex.mwFreq*1E9);
        gw      = gcenter*Ex.ExciteWidth/28/Ex.Field /sqrt(2*log(2));
        ak = exp(-2*((geff-gcenter)./gw).^2);
    else
        ak = ones(1,ntheta*nphi);
    end
    
    Treshold = safeget(Si, 'Treshold', 0.01);
    MA = []; EA=MA; MB=MA; EB=MA; projA=[];
    
    for ii=1:length(phi)
        if ak(ii)<Treshold, continue; end
        HFI_A = HFIx*Gx(ii)./geff(ii)*1/2 + ...
            HFIy*Gy(ii)./geff(ii)*1/2 + ...
            HFIz*Gz(ii)./geff(ii)*1/2;

        fhhf = HFIx*Gx(ii)./geff(ii)*1/2 + ...
            HFIy*Gy(ii)./geff(ii)*1/2 + ...
            HFIz*Gz(ii)./geff(ii)*1/2 + 10000*(angle(1,ii)*iix + angle(2,ii)*iiy + angle(3,ii)*iiz);
        
        projA(end+1) = max(eig(fhhf))-10000;
        
        Zeeman_A = -omega_l*(angle(1,ii)*Ix + angle(2,ii)*Iy + angle(3,ii)*Iz);
        
        Ham_A = Zeeman_A - HFI_A + Quad;
        Ham_B = Zeeman_A + HFI_A + Quad;
        
        [VA, DA] = eig(Ham_A);
%         mm = real(VA'*Ix*VA).^2;
        dd = diag(DA);
        d1 = dd(:, onn)-dd(:, onn).';
        MA(end+1, :) = diag(d1, +1)*0+ak(ii);
        EA(end+1, :) = abs(real(diag(d1, +1)));
        
        [VB, DB] = eig(Ham_B);
%         mm = real(VB'*Ix*VB).^2;
        dd = diag(DB);
        d1 = dd(:, onn)-dd(:, onn).';
        MB(end+1, :) = diag(d1, +1)*0+ak(ii);
        EB(end+1, :) = abs(real(diag(d1, +1)));
    end
 
    Thole = safeget(Ex, 'THole', 0); % get times like tinv or tau [us]
    hole = safeget(Ex, 'Hole', 0);
    
    if hole | Thole
        if ~Thole
            Thole = 1/hole;
        end
        etha_s = abs(projA*Thole); 
        etype = safeget(Ex, 'Scheme', 'Davies');
        if strcmp(etype,'Davies')
            profile = sqrt(2)*etha_s./(etha_s.^2+0.5);
        else
             profile = 1/4*(1-cos(2*pi*etha_s));
        end
    else
        profile = projA*0+1;
    end
    oon = MA(1, :)*0+1;
    MA = MA.*profile(oon, :).';
    MB = MB.*profile(oon, :).';
    stepRF = (Ex.Range(2)-Ex.Range(1))/(Ex.nPoints-1);
    Freqs = [Ex.Range(1) : stepRF : Ex.Range(2)];
        
    tSpec1 = zeros(1,Ex.nPoints);
    tSpec2 = zeros(1,Ex.nPoints);
    
%     tSpec1 = binning(tSpec1, Freqs, w1(1,:), profile(1, :).*sangle.*ak);
%     tSpec2 = binning(tSpec2, Freqs, w1(2,:), profile(2, :).*sangle.*ak);
%     tSpec3 = binning(tSpec3, Freqs, w1(3,:), profile(3, :).*sangle.*ak);
%     tSpec4 = binning(tSpec4, Freqs, w1(4,:), profile(4, :).*sangle.*ak);

    tSpec1 = binning(tSpec1, Freqs, EA(:), MA(:));
    tSpec2 = binning(tSpec2, Freqs, EB(:), MB(:));
    tSpec = tSpec1+tSpec2;
    if Sys.lwEndor~=0
        mid = round(Ex.nPoints/2)+1;
        Line = lorentzian(Freqs,Freqs(mid),Sys.lwEndor);
        tdDecay = ifft(fftshift(Line*stepRF));
        td = ifft(tSpec,[],2).*tdDecay;
        tSpec = abs(fft(td,[],2))*Ex.nPoints;
    end
    hm = safeget(Ex, 'Harmonic', 0);
    if hm > 0,
        tSpec = pseumod(Freqs, tSpec, Sys.lwEndor/10, hm);
    end
    Spec = Spec+tSpec;
end

return

% input is (0,0)  (1,1)  (2,2)  (1,0)  (2,0)  (2,1)
%            1      2      3      4      5      6
% see Muha article for details
function w = p3eq(H)
sh = (H(1,:)+H(2,:)+H(3,:))/3;
H(1,:) = H(1,:) - sh;
H(2,:) = H(2,:) - sh;
H(3,:) = H(3,:) - sh;
W1=H(4,:).*conj(H(4,:));
W2=H(5,:).*conj(H(5,:));
W3=H(6,:).*conj(H(6,:));
Q=real(H(1,:)).*real(H(2,:)).*real(H(3,:))-...
    real(H(1,:)).*W3-real(H(2,:)).*W2-real(H(3,:)).*W1+...
    2*(real(H(4,:)).*(real(H(5,:)).*real(H(6,:))+imag(H(5,:).*imag(H(6,:)))))+...
    2*(imag(H(4,:)).*(imag(H(5,:)).*real(H(6,:))-real(H(5,:).*imag(H(6,:)))));
P=real(H(1,:)).*real(H(2,:))+real(H(1,:)).*real(H(3,:))+...
    real(H(2,:)).*real(H(3,:)) -W1-W2-W3;
ABS_P=abs(P);
C=(3./ABS_P).^1.5.*Q/2;
acosC=acos(C); 
P1=sqrt(4*ABS_P/3);
w=[P1.*cos(acosC/3)+sh;P1.*cos((acosC+4*pi)/3)+sh;P1.*cos((acosC+2*pi)/3)+sh];

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


