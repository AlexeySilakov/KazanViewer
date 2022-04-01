% Simple and fast routine for calculating EPR spectra
% of one S = 1/2 and multiple I=1/2 nucleus.
% for the case Zeeman >> Hyperfine
% 
% [x,y] = epr05(Sys,Exp,Opt)
% 
% Partially compatible with EasySpin TM
% Parameters list: 
%  Sys: g,gn,A,Apa,lwEndor
%  Exp: Field, Range, mwFreq
%  Opt: ThetaRange,PhiRange,nKnots,nPoints
%  External dependencies:  safeget, binning
%  not listed: Exp.as - assymetery

% Boris Epel, MPI for Bioinorganic Chemistry, 2002-05

%==========================================================================
function [Field, Spec] = epr05(Sy,Ex,Si)
%==========================================================================
if nargin < 2
    error('Usage: sugar(Sys,Exp,Opt).');
end

% tm = clock;

if length(Sy.S)>1 | Sy.S~=0.5, disp('EPR05: Sys.S > 1 is not supported.'); return; end;
if isfield(Sy, 'Nucs'), [Sy.I,Sy.gn] = nucdata(Sy.Nucs); end
nucs = length(Sy.I);
if sum(Sy.I==0.5*ones(1,nucs))~=nucs, disp('EPR05: Sys.I > 1 is not supported.'); return; end;

planck = 6.6261e-034;
nmagn = 5.0508e-027;
Ex.Range=safeget(Ex, 'Range', [0 100]);
stepF = (Ex.Range(2)-Ex.Range(1))/(Ex.nPoints-1);
Field = linspace(Ex.Range(1), Ex.Range(2), Ex.nPoints);

for ii=1:nucs
  R = erot(Sy.Apa(ii,:));
  A(ii,:,:) = R * diag(Sy.A(ii,:))*R.';
end

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

cphi = cos(phi);
sphi = sin(phi);

angle(1,:)=reshape(stheta(:)*cphi(:).', 1, nphi*ntheta);
angle(2,:)=reshape(stheta(:)*sphi(:).', 1, nphi*ntheta);
angle(3,:) = reshape(ctheta(:)*ones(1,nphi), 1, nphi*ntheta);

sangle = (reshape(stheta(:)*ones(1,nphi), 1, nphi*ntheta));

Gx=angle(1,:)*Sy.g(1);
Gy=angle(2,:)*Sy.g(2);
Gz=angle(3,:)*Sy.g(3);

geff=sqrt(Gx.*Gx+Gy.*Gy+Gz.*Gz);
Fld_g1 = planck*Ex.mwFreq/bmagn*1E12;
       
% separate output 
soutp = strcmp(safeget(Si, 'Output','summed'), 'separate');

ww  = zeros(1,nphi*ntheta);
proj = allproj(nucs,[-1 1]*0.25);
nproj = size(proj,1);
if(soutp)
  Spec = zeros(nproj,Ex.nPoints);
else
  Spec = zeros(1,Ex.nPoints);
end
for ii=1:nucs
    omega_l = Fld_g1/2/planck/1E9 * Sy.gn(ii) * nmagn;
    B = omega_l * angle;

    Axyz(1,:)=(Gx*A(ii,1,1) + Gy*A(ii,1,2) + Gz*A(ii,1,3))./geff;
    Axyz(2,:)=(Gx*A(ii,2,1) + Gy*A(ii,2,2) + Gz*A(ii,2,3))./geff;
    Axyz(3,:)=(Gx*A(ii,3,1) + Gy*A(ii,3,2) + Gz*A(ii,3,3))./geff;
    
    Hhpf(1,:)= Axyz(3,:);
    Hhpf(2,:)= Axyz(1,:) - j* Axyz(2,:);
    
    ww(ii,:)  = 2 * sqrt(Hhpf(1,:).*Hhpf(1,:) + Hhpf(2,:).*conj(Hhpf(2,:)));
end

cc = 2*bmagn/planck*1E-9;
if(soutp)
    for ii=1:nproj
        ww1 = zeros(1,nphi*ntheta);
        for jj=1:nucs
            ww1(1,:)=ww1(1,:) + ww(jj,:)*proj(ii,jj)/cc;
        end
        Spec(ii,:) = binning(Spec(ii,:), Field, Fld_g1./geff+ww1, sangle);
    end
else
    for ii=1:nproj
        ww1 = zeros(1,nphi*ntheta);
        for jj=1:nucs
            ww1(1,:)=ww1(1,:) + ww(jj,:)*proj(ii,jj)/27.992;
        end
        Spec = binning(Spec, Field, Fld_g1./geff+ww1, sangle);
    end
end

if Sy.lw
    mid = round(Ex.nPoints/2)+1;
    Line = lshape(Field,Field(mid),Sy.lw);
    tdDecay = ifft(fftshift(Line*stepF));
    for ii=1:size(Spec,1)
        td = ifft(Spec(ii,:),[],2).*tdDecay;
        Spec(ii,:) = abs(fft(td,[],2))*Ex.nPoints;
    end
end

hm = safeget(Ex, 'Harmonic', 0);
if hm > 0, 
    Spec = pseumod(Field, Spec, Sy.lw/10, hm);
end

% tm = clock - tm;
% disp(sprintf('%f min %f sec', tm(5),tm(6)))
Field = Field*10;
return
