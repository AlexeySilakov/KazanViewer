function [y, ax]= mimstau(varargin)

% Simple and fast routine for calculating tau dependent Mims ENDOR spectra
% of one I=1/2 nucleus.
% 
% [y, ax] = mimstau(Sys,Exp,Opt)
% y - data
% ax - struct containing axes (ax.x, ax.y, ax.xlabel, ax.ylabel, etc)
%
% Partially compatible with EasySpin TM
% Parameters list: 
%  Sys: g,gn,A,Apa,lwEndor
%  Exp: Field, Range, mwFreq, Tau [us], nTauPoints
%  Opt: ThetaRange,PhiRange,nKnots,nPoints
%  External dependencies:  safeget, binning
%  not listed: Exp.as - assymetery

% Based on the "sugar" script of Boris Epel
% Alexey Silakov, MPI for Bioinorganic Chemistry, 2009

%==========================================================================
% function [Freqs, Spec] = sugar(Sy,Ex,Si)
%==========================================================================

if nargin < 2
    error('Usage: mimstau(Sys,Exp), mimstau(Sys,Exp,Opt).');
end

Sy = varargin{1};
Ex = varargin{2};
Si = [];
if nargin>2,
    Si = varargin{3};
end

% tm = clock;

if isfield(Sy, 'Nucs'), [Sy.I,Sy.gn] = nucdata(Sy.Nucs); end
if Sy.I>1/2, disp('mimstau: Sys.I>=1 is not supported. Sys.I=1/2 will be calculated !'); end


nucgfactor=Sy.gn;
planck = 6.6261e-034;
nmagn = 5.0508e-027;

if isfield(Ex, 'gField'),
   if ~Ex.mwFreq & ~isfield(Ex, 'Field'), ...
           error('ESEEM: there is no Exp.mwFreq to calculate B0');
   end
   Ex.Field = (Ex.mwFreq*1e+9*planck)/(Ex.gField*bmagn*1e-3);  % mT
else 
   Ex.Field = safeget(Ex, 'Field', 0);
end

Ex.Range=safeget(Ex, 'Range', [0 100]);
R = erot(Sy.Apa);
A = R * diag(Sy.A)*R.';

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
% sangle = ones(1,nphi*ntheta);

omega_l = Ex.Field/planck /1E9 * nucgfactor * nmagn; % [MHz]

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

H_1 = Hhpf*2;
% H1_1 = -Hhpf;
w_1 = 2 * sqrt(H_1(1,:).*H_1(1,:) + H_1(2,:).*conj(H_1(2,:)));
% w1_1 = 2 * sqrt(H1_1(1,:).*H1_1(1,:) + H1_1(2,:).*conj(H1_1(2,:)));

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

mode = safeget(Si, 'MimsMode', 'fft');
nTauPoints = safeget(Si, 'nTauPoints', Ex.nPoints);
switch mode
    case 'time'
        Spec = zeros(nTauPoints,Ex.nPoints);
        if isfield(Ex, 'tau')
            if size(Ex.tau)>1
               tau_cur = linspace(Ex.tau(1), Ex.tau(2), nTauPoints);
            else
               tau_cur = linspace(0, Ex.tau(1), nTauPoints);
            end
        else
            % MaxFreq is explicitly specified
            Ex.MaxFreq = safeget(Ex, 'MaxFreq', Ex.Range(2));
            tau_cur = linspace(0, (nTauPoints-1)/Ex.MaxFreq/2, nTauPoints);
        end
        
        for ii = 1:nTauPoints
            profile = sin(2*pi*(w_1)*tau_cur(ii)).^2;

            Spec(ii, :) = binning(Spec(ii, :), Freqs, w, sangle.*ak*(1-as).*profile);
            Spec(ii, :) = binning(Spec(ii, :), Freqs, w1, sangle.*ak*(1+as).*profile);
            
            if Sy.lwEndor
                mid = round(Ex.nPoints/2)+1;
                Line = lshape(Freqs,Freqs(mid),Sy.lwEndor);
                tdDecay = ifft(fftshift(Line*stepRF));
                td = ifft(Spec,[],2).*tdDecay;
                Spec = abs(fft(td,[],2))*Ex.nPoints;
            end

        end
        
        y = Spec;
        ax.x = Freqs;
        ax.xlabel = 'ENDOR Frequency, MHz';
        ax.y = tau_cur;
        ax.ylabel = 'tau, us';
        
    case 'fft'
        Spec = zeros(nTauPoints,Ex.nPoints);
        
        ind11 = 1+floor( (w-Ex.Range(1))/(abs( (Ex.Range(2)-Ex.Range(1))/(Ex.nPoints-1))) );
        ind12 = 1+floor( (w1-Ex.Range(1))/(abs( (Ex.Range(2)-Ex.Range(1))/(Ex.nPoints-1))) );        
        
        if isfield(Ex, 'tau')
            % MaxFreq is based on the Nyquist frequency
            if size(Ex.tau)>1
                Ex.MaxFreq = (nTauPoints -1)/ (Ex.tau(2)-Ex.tau(1)) /2;
            else
                Ex.MaxFreq = (nTauPoints -1) / Ex.tau; % assume that it starts from tau=0
            end
        else
            % MaxFreq is explicitly specified
            Ex.MaxFreq = safeget(Ex, 'MaxFreq', Ex.Range(2));
        end
        
        ind2 = 1+floor( (w_1)/(abs(Ex.MaxFreq/(nTauPoints-1))) );
               
        AMP = sangle.*ak;
        
        find11 = find(ind11(:)>0 & ind11(:) <= Ex.nPoints);
        ind1_1 = ind11(find11);
        ind2_1 = ind2(find11);
        AMP_1 = AMP(find11);
                
        find12 = find(ind12(:)>0 & ind12(:) <= Ex.nPoints);
        ind1_2 = ind12(find12);
        ind2_2 = ind2(find12);
        AMP_2 = AMP(find12);
              
        find2_1 = find(ind2_1(:)>0 & ind2_1(:) <= nTauPoints);
        find2_2 = find(ind2_2(:)>0 & ind2_2(:) <= nTauPoints);
        
        if ~isempty(find2_1)
            Spec = bin2D(Spec, ind1_1(find2_1), ind2_1(find2_1), AMP_1(find2_1)*(1-as));
        end
        if ~isempty(find2_2)
            Spec = bin2D(Spec, ind1_2(find2_2), ind2_2(find2_2), AMP_2(find2_2)*(1+as));
        end
        
        if safeget(Sy, 'lwEndor', 0)
            
                mid = round(Ex.nPoints/2)+1;
                Line = lshape(Freqs,Freqs(mid),Sy.lwEndor);
                tdDecay = ifft(fftshift(Line*stepRF));
                td = ifft(Spec,[],2).*tdDecay(ones(nTauPoints, 1), :);
                Spec = abs(fft(td,[],2))*Ex.nPoints;
        end
        ax.y = Freqs.';
        ax.ylabel = 'ENDOR Frequency, MHz';
        ax.x = linspace(0, Ex.MaxFreq, nTauPoints).';
        ax.xlabel = '"Tau" Frequency, MHz';  
        
        if safeget(Sy, 'lw', 0)
            mid = round(nTauPoints/2)+1;
            Line = lshape(ax.x, ax.x(mid),Sy.lw);
            tdDecay = ifft(fftshift(Line* (ax.x(2)-ax.x(1)) ));
            td = ifft(Spec,[],1).*tdDecay(:, ones(Ex.nPoints, 1));
            Spec = abs(fft(td,[],1))*Ex.nPoints;
        end
        y = Spec;

end

hm = safeget(Ex, 'Harmonic', 0);
if hm > 0, 
    Spec = pseumod(Freqs, Spec, Sy.lwEndor/10, hm);
end



% tm = clock - tm;
% disp(sprintf('%f min %f sec', tm(5),tm(6)))

return
