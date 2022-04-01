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
function [Freqs, Spec] = sugar1a(Sys,Ex,Si)
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

    if Sy.I==1/2, disp('Warning: for case Sys.I=1/2 ''sugar'' may be used.');
    elseif Sy.I>1
        error('case Sys.I>1 is not supported.');
    end
    
    R = erot(safeget(Sy, 'Apa', [0,0,0]));
    A = R * diag(Sy.A)*R.';
    
    % ntheta = Si.nKnots*2;
    % nphi = Si.nKnots*4;
    
    % PhiRange = safeget(Si, 'PhiRange', [0,2*pi]);
    % ThetaRange = safeget(Si, 'ThetaRange', [0, pi]);
    % steptheta = (ThetaRange(2) - ThetaRange(1))/(ntheta - 1);
    % stepphi = (PhiRange(2) - PhiRange(1))/(nphi - 1);
    % %
    % theta = ThetaRange(1):steptheta:ThetaRange(2);
    % phi   = PhiRange(1):stepphi:PhiRange(2);
    
    
    omega_l = Ex.Field/planck /1E9 * Sy.gn * nmagn;
    
    B = omega_l * angle;
    
    % G = g.*reshape([a_sthetacphi(:) a_sthetasphi(:) a_ctheta(:)], 3, ntheta*nphi);
    
    Axyz(1,:)=(Gx*A(1,1) + Gy*A(1,2) + Gz*A(1,3))./geff;
    Axyz(2,:)=(Gx*A(2,1) + Gy*A(2,2) + Gz*A(2,3))./geff;
    Axyz(3,:)=(Gx*A(3,1) + Gy*A(3,2) + Gz*A(3,3))./geff;
    
    sq2 = sqrt(2);
    
    Hzeeman(1,:) = -B(3,:);                        % (0,0)
    Hzeeman(2,:) = 0;                              % (1,1)
    Hzeeman(3,:) = B(3,:);                         % (2,2)
    Hzeeman(4,:) = (-B(1,:) - j*B(2,:))/2*sq2;     % (1,0)
    Hzeeman(5,:) = 0;                              % (2,0)
    Hzeeman(6,:) = (-B(1,:) - j*B(2,:))/2*sq2;     % (2,1)
    
    Hhpf(1,:) = Axyz(3,:);
    Hhpf(2,:) = 0;
    Hhpf(3,:) = -Axyz(3,:);
    Hhpf(4,:) = Axyz(1,:)/2*sq2 + j*Axyz(2,:)/2*sq2;
    Hhpf(5,:) = 0;
    Hhpf(6,:) = Axyz(1,:)/2*sq2 + j*Axyz(2,:)/2*sq2;
    
    Hhpf=Hhpf.*0.5;  %% e-spin
    
    % not tested yet
    if isfield(Sy,'Q')
        sz= size(angle,2);
        RQ = erot(safeget(Sy, 'Qpa', [0,0,0]));
        Q = RQ * diag(Sy.Q)*RQ.';
        HQuad(1,1:sz) = Q(3,3)/2;
        HQuad(2,1:sz) = -Q(3,3);
        HQuad(3,1:sz) = Q(3,3)/2;
        HQuad(4,1:sz) = (Q(1,3)+j*Q(2,3))*sq2/2;
        HQuad(5,1:sz) = (Q(1,1)-Q(2,2))/2 + j*Q(1,2);
        HQuad(6,1:sz) = -(Q(1,3)+j*Q(2,3))*sq2/2;
        H =  Hzeeman + Hhpf + HQuad;
        H1 = Hzeeman - Hhpf + HQuad;
    else
        H =  Hzeeman + Hhpf;
        H1 = Hzeeman - Hhpf;
    end
    
    e1 = p3eq(H);
    e2 = p3eq(H1);
    w(1,:) = abs(e1(1,:)-e1(2,:));
    w(2,:) = abs(e1(2,:)-e1(3,:));
    w(3,:) = abs(e2(1,:)-e2(2,:));
    w(4,:) = abs(e2(2,:)-e2(3,:));
    
    %%%%%%%%%%%%% alsi ... proper treatment of the blindspons
    Thole = safeget(Ex, 'THole', 0); % get times like tinv or tau [us]
    hole = safeget(Ex, 'Hole', 0);
    
    if hole | Thole
        
        
        H =  + Hhpf;
        H1 = - Hhpf;
        e1 = p3eq(H);
        e2 = p3eq(H1);
        w1(1,:) = real(e1(1,:)-e1(2,:));
        w1(2,:) = real(e1(2,:)-e1(3,:));
        w1(3,:) = real(e2(1,:)-e2(2,:));
        w1(4,:) = real(e2(2,:)-e2(3,:));
        
        if ~Thole
            Thole = 1/hole;
        end
        etha_s = abs(2*(abs(w1)-omega_l)*Thole); %w1 - is just a half of HF, so x2
        etype = safeget(Ex, 'Scheme', 'Davies');
        if strcmp(etype,'Davies')
            %lshape(Freqs, omega_l, hole, 'lorentz');
            profile = sqrt(2)*etha_s./(etha_s.^2+0.5);
        else
            
            profile = 1/4*(1-cos(2*pi*etha_s));
        end
        %     Spec = Spec.*profile;
    else
        profile = w*0+1;
    end
    
    
    stepRF = (Ex.Range(2)-Ex.Range(1))/(Ex.nPoints-1);
    Freqs = [Ex.Range(1) : stepRF : Ex.Range(2)];
    
    if Ex.ExciteWidth > 0
        gcenter = fld2g(Ex.Field*1E-3, Ex.mwFreq*1E9);
        gw      = gcenter*Ex.ExciteWidth/28/Ex.Field /sqrt(2*log(2));
        ak = exp(-2*((geff-gcenter)./gw).^2);
    else
        ak = ones(1,ntheta*nphi);
    end
    
    tSpec1 = zeros(1,Ex.nPoints);
    tSpec2 = zeros(1,Ex.nPoints);
    tSpec3 = zeros(1,Ex.nPoints);
    tSpec4 = zeros(1,Ex.nPoints);
%     tSpec1 = binning(tSpec1, Freqs, w1(1,:), profile(1, :).*sangle.*ak);
%     tSpec2 = binning(tSpec2, Freqs, w1(2,:), profile(2, :).*sangle.*ak);
%     tSpec3 = binning(tSpec3, Freqs, w1(3,:), profile(3, :).*sangle.*ak);
%     tSpec4 = binning(tSpec4, Freqs, w1(4,:), profile(4, :).*sangle.*ak);

    tSpec1 = binning(tSpec1, Freqs, w(1,:), profile(1, :).*sangle.*ak);
    tSpec2 = binning(tSpec2, Freqs, w(2,:), profile(2, :).*sangle.*ak);
    tSpec3 = binning(tSpec3, Freqs, w(3,:), profile(3, :).*sangle.*ak);
    tSpec4 = binning(tSpec4, Freqs, w(4,:), profile(4, :).*sangle.*ak);
    tSpec = tSpec1+tSpec2+tSpec3+tSpec4;
    if ~sum(tSpec), Spec = Spec+tSpec; continue; end
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
