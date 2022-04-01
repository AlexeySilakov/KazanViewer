% Simple and fast routine for calculating ENDOR spectra
% of one I=1/2 nucleus with isotropic HF coupling (first order only).
% 
% [x,y] = cw_endor(Sys,Exp,Opt)
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
function [Freqs, Spec] = cw_endor(Sy,Ex,Si)
%==========================================================================
if nargin < 2
    error('Usage: sugar(Sys,Exp,Opt).');
end

if length(Sy.S)>1 | Sy.S~=0.5, disp('EPR05: Sys.S > 1 is not supported.'); return; end;
if isfield(Sy, 'Nucs'), [Sy.I,Sy.gn] = nucdata(Sy.Nucs); end
nucs = length(Sy.I);
nucgfactor=Sy.gn;

if(length(Sy.A)~=nucs) error('Nimber of Nucs not equal to number of A'); end

nmagn = 5.0508e-027;
planck = 6.6261e-034;
omega_l = Ex.Field/planck /1E9 * nucgfactor * nmagn;    

stepRF = (Ex.Range(2)-Ex.Range(1))/(Ex.nPoints-1);
Freqs = [Ex.Range(1) : stepRF : Ex.Range(2)];

Spec = zeros(1,Ex.nPoints);

% include <1> or not <0> HF enhancement
HFenhancement = 1;

if HFenhancement
    for ii=1:nucs
        Spec = binning(Spec, Freqs, omega_l(ii)+Sy.A(ii)/2, abs(omega_l(ii)+Sy.A(ii))/omega_l(ii));
        Spec = binning(Spec, Freqs, omega_l(ii)-Sy.A(ii)/2, abs(omega_l(ii)-Sy.A(ii))/omega_l(ii));
    end
else
    for ii=1:nucs
        Spec = binning(Spec, Freqs, omega_l(ii)+Sy.A(ii)/2, 1);
        Spec = binning(Spec, Freqs, omega_l(ii)-Sy.A(ii)/2, 1);
    end
end
Sy.lwEndor = safeget(Sy, 'lwEndor', 0.1);

if Sy.lwEndor
    mid = round(Ex.nPoints/2)+1;
    Line = lshape(Freqs,Freqs(mid),Sy.lwEndor, 0, safeget(Si, 'LineShape', 0));
    tdDecay = ifft(fftshift(Line*stepRF));
    td = ifft(Spec,[],2).*tdDecay;
    Spec = abs(fft(td,[],2))*Ex.nPoints;
end

hm = safeget(Ex, 'Harmonic', 0);
if hm > 0, 
    Spec = pseumod(Freqs, Spec, Sy.lwEndor/10, hm);
end

return
