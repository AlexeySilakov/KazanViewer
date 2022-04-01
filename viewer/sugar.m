% Simple and fast routine for calculating ENDOR spectra
% of one I=1/2 nucleus.
% 
% [x,y] = sugar(Sys,Exp,Opt)
% 
% Partially compatible with EasySpin TM
% Parameters list: 
%  Sys: g,gn,A,Apa,lwEndor
%  Exp: Field, Range, mwFreq
%  Opt: ThetaRange,PhiRange,nKnots,nPoints
%  External dependencies:  safeget, binning
%  not listed: Exp.asym - assymetery

% Boris Epel, MPI for Bioinorganic Chemistry, 2002-05

%==========================================================================
function [Freqs, Spec] = sugar(Sy,Ex,Si)
%==========================================================================
if nargin < 2
    error('Usage: sugar(Sys,Exp,Opt).');
end

% tm = clock;

if isfield(Sy, 'Nucs'), [Sy.I,Sy.gn] = nucdata(Sy.Nucs); end
if Sy.I==1, disp('sugar: for case Sys.I=1 ''sugar1'' must be used.');
elseif Sy.I>1
    disp('sugar: case Sys.I>1 is not supported.');
end
nucgfactor=Sy.gn;
planck = 6.6261e-034;
nmagn = 5.0508e-027;

Ex.Range=safeget(Ex, 'Range', [0 100]);
R = erot(Sy.Apa);
A = R * diag(Sy.A)*R.';

if isfield(Si, 'OriSelMatrix')
    phi = Si.OriSelMatrix(:, 1);
    theta = Si.OriSelMatrix(:, 2);
    angle(1,:)=sin(theta).*cos(phi);
    angle(2,:)=sin(theta).*sin(phi);
    angle(3,:) = cos(theta);
    sangle = Si.OriSelMatrix(:, 3).';
    Ex.ExciteWidth  = 0;
    ntheta = numel(sangle);
    nphi = 1;
else
    ntheta = Si.nKnots*2;
    nphi = Si.nKnots*4;
    
    PhiRange = safeget(Si, 'PhiRange', [0,2*pi]);
    ThetaRange = safeget(Si, 'ThetaRange', [0, pi]);
    steptheta = (ThetaRange(2) - ThetaRange(1))/(ntheta - 1);
    stepphi = (PhiRange(2) - PhiRange(1))/(nphi - 1);
    
    theta = ThetaRange(1):steptheta:ThetaRange(2);
    phi   = PhiRange(1):stepphi:PhiRange(2);
    
    ctheta = cos(theta);
    stheta = sin(theta);
    
    % ctheta = linspace(cos(ThetaRange(1)), cos(ThetaRange(2)), ntheta);
    % stheta = sqrt(1-ctheta.*ctheta);
    
    cphi = cos(phi);
    sphi = sin(phi);
    
    angle(1,:)=reshape(stheta(:)*cphi(:).', 1, nphi*ntheta);
    angle(2,:)=reshape(stheta(:)*sphi(:).', 1, nphi*ntheta);
    angle(3,:) = reshape(ctheta(:)*ones(1,nphi), 1, nphi*ntheta);
    
    sangle = (reshape(stheta(:)*ones(1,nphi), 1, nphi*ntheta));
end
% sangle = ones(1,nphi*ntheta);

omega_l = Ex.Field/planck /1E9 * nucgfactor * nmagn;

B = omega_l * angle;

Gx=angle(1,:)*Sy.g(1);
Gy=angle(2,:)*Sy.g(2);
Gz=angle(3,:)*Sy.g(3);

geff=sqrt(Gx.*Gx+Gy.*Gy+Gz.*Gz);
       
% G = g.*reshape([a_sthetacphi(:) a_sthetasphi(:) a_ctheta(:)], 3, ntheta*nphi);

Axyz(1,:)=(Gx*A(1,1) + Gy*A(1,2) + Gz*A(1,3))./geff;
Axyz(2,:)=(Gx*A(2,1) + Gy*A(2,2) + Gz*A(2,3))./geff;
Axyz(3,:)=(Gx*A(3,1) + Gy*A(3,2) + Gz*A(3,3))./geff;

Hzeeman(1,:) = -B(3,:);
Hzeeman(2,:) = -B(1,:) + j*B(2,:);

Hzeeman = Hzeeman.*0.5;

Hhpf(1,:)= Axyz(3,:);
Hhpf(2,:)= Axyz(1,:) - j* Axyz(2,:);

Hhpf=Hhpf.*0.25;

H =  Hzeeman + Hhpf;
H1 = Hzeeman - Hhpf;

w = 2 * sqrt(H(1,:).*H(1,:) + H(2,:).*conj(H(2,:)));
w1 = 2 * sqrt(H1(1,:).*H1(1,:) + H1(2,:).*conj(H1(2,:)));
% w = abs(2 * (H(1,:) + H(2,:)));
% w1 = abs(2 * (H1(1,:) + H1(2,:)));

stepRF = (Ex.Range(2)-Ex.Range(1))/(Ex.nPoints-1);
Freqs = [Ex.Range(1) : stepRF : Ex.Range(2)];

if Ex.ExciteWidth > 0
    gcenter = fld2g(Ex.Field*1E-3, Ex.mwFreq*1E9);
    gw      = gcenter*Ex.ExciteWidth/28/Ex.Field /sqrt(2*log(2));
    ak = exp(-2*((geff-gcenter)./gw).^2);
else
    ak = ones(1,ntheta*nphi);
end

as = safeget(Ex, 'asym', 0);

Spec = zeros(1,Ex.nPoints);
Spec = binning(Spec, Freqs, w, sangle.*ak*(1-as));
Spec = binning(Spec, Freqs, w1, sangle.*ak*(1+as));

if Sy.lwEndor
    mid = round(Ex.nPoints/2)+1;
    Line = lshape(Freqs,Freqs(mid),Sy.lwEndor);
    tdDecay = ifft(fftshift(Line*stepRF));
    td = ifft(Spec,[],2).*tdDecay;
    Spec = abs(fft(td,[],2))*Ex.nPoints;
end

hole = safeget(Ex, 'Hole', 0);
if hole > 0
    etype = safeget(Ex, 'Scheme', 'Davies');
    if strcmp(etype,'Davies')
       profile = 1-renorm(1./((Freqs-omega_l).^2+ hole.^2/4));
%       profile = lorentzian(Freqs, omega_l, hole, 'lorentz');    
%       profile = 1 - renorm(profile);
    else
      profile = sin(2*pi*(Freqs-omega_l)/hole).^2;  
    end
    Spec = Spec.*profile;
end

hm = safeget(Ex, 'Harmonic', 0);
if hm > 0, 
    Spec = kv_pseumod(Freqs, Spec, Sy.lwEndor/10, hm);
end

% tm = clock - tm;
% disp(sprintf('%f min %f sec', tm(5),tm(6)))

return
