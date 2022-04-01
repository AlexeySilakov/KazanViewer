% Simple and fast routine for calculating EPR and ENDOR spectra
% of spin polarised triplets including HF couplings.
% 
% [x,y] = smetana(Sys,Exp,Opt)
% 
% Partially compatible with EasySpin TM
% Parameters list: 
%  Sys: g,gn,A,Apa,lwEndor
%  Exp: Field, Range, mwFreq
%  Opt: ThetaRange,PhiRange,nKnots,nPoints
%  External dependencies:  safeget, binning
%  not listed: Exp.as - assymetery

% Alexey Silakov, MPI for Bioinorganic Chemistry, 2010-11
%     based on sugar by Borisd Epel
%==========================================================================
function [Freqs, Spec] = smetana(Sy,Ex,Op)
%==========================================================================
if nargin < 2
    error('Usage: smetana(Sys,Exp,Opt).');
end

% tm = clock;

if isfield(Sy, 'Nucs'), [Sy.I,Sy.gn] = nucdata(Sy.Nucs); end

if Sy.S > 1, error('SMETANA: Sys.S is larger than 1'); end

if ~(safeget(Op, 'EPR', 0))
    if Sy.I==1, disp('smetana: sorry, only I=1/2 is implemented');
    elseif Sy.I>1
        disp('moloko: case Sys.I>1 is not supported.');
    end
end

%%% simulation for a few nuclei
Parameters ={'A', 'Apa', 'Q', 'Qpa', 'Nucs', 'gn', 'I', 'nNucs'};

Sy.Apa = safeget(Sy, 'Apa', [0, 0, 0])*pi/180; % degree -> rad
Sy.Qpa = safeget(Sy, 'Qpa', [0, 0, 0])*pi/180; % degree -> rad
% lets find number of simmulaitons (parameters can be common, it means
% size(XX, 1)==1)
Sy.nNucs = safeget(Sy, 'nNucs', 1); % number of equivalent nuclei
sim_size = 1;
isfield_num = [];
for pk = 1:length(Parameters)
    if isfield(Sy, Parameters{pk})        
        tPar = getfield(Sy, Parameters{pk});
        isfield_num(end+1) = pk;        
        if size(tPar, 1)>1
            sim_size = size(tPar, 1);
        end
    end
end
if sim_size>1
    if size(Sy.I, 2) > 1,
        Sy.I = Sy.I.';
    end
    if size(Sy.gn, 2) > 1,
        Sy.gn = Sy.gn.';
    end
    for kk = 1:length(isfield_num)
        pk = isfield_num(kk); % kind of optimisation
        tPar = getfield(Sy, Parameters{pk});
        if size(tPar, 1)==1
            Sy = setfield(Sy, Parameters{pk}, tPar(ones(sim_size, 1), :));
        elseif size(tPar, 1)~=sim_size
%             error({['Parameter Sys.', Parameters{pk}, ' has wrong number of rows.']; ['Correct number is ', num2str(sim_size)]});
              error(['Parameter Sys.', Parameters{pk}, ' has wrong number of rows. Correct number is ', num2str(sim_size)]);
        end
    end
end

if isfield(Sy, 'Aiso'),
    Sy.A = safeget(Sys, 'A', 0) + Sy.Aiso;
end
% 


planck = 6.6261e-034;
nmagn = 5.0508e-027;

if safeget(Op, 'EPR', 0)
    Ex.Range=safeget(Ex, 'Range', [300 400]);
    [Freqs, Spec] = Local_orisel(Sy, Ex, Op);
    hm = safeget(Ex, 'Harmonic', 0);
    return;
else
    Ex.Range=safeget(Ex, 'Range', [0 100]);
    [phi, theta, wa] = Local_orisel(Sy, Ex, Op);
end

ctheta = cos(theta);
stheta = sin(theta);

    % ctheta = linspace(cos(ThetaRange(1)), cos(ThetaRange(2)), ntheta);
    % stheta = sqrt(1-ctheta.*ctheta);

    cphi = cos(phi);
    sphi = sin(phi);

    nangs = length(phi);
    
% NucsToCalc = safeget(Op, 'Nuclei', 1:size(Sy.A, 1));
if isfield(Op, 'Nuclei')
    NucsToCalc = Op.Nuclei;
else
    NucsToCalc = 1:size(Sy.A, 1);
end

Spec = zeros(1,Ex.nPoints);
stepRF = (Ex.Range(2)-Ex.Range(1))/(Ex.nPoints-1);
Freqs = [Ex.Range(1) : stepRF : Ex.Range(2)];

Sa = Sy.S;

nSysDim = 2*Sa+1;
MS = [-Sa:1:Sa];

for cA = 1:length(NucsToCalc)
   
    for cMS = 1:(nSysDim-1)
         ak = wa(:, cMS);
        nucidx = NucsToCalc(cA);
        nucgfactor=Sy.gn(nucidx);
        R = erot(Sy.Apa(nucidx, :));
        A = R * diag(Sy.A(nucidx, :))*R.';

        % angle(1,:)=reshape(stheta(:)*cphi(:).', 1, nphi*ntheta);
        % angle(2,:)=reshape(stheta(:)*sphi(:).', 1, nphi*ntheta);
        % angle(3,:) = reshape(ctheta(:)*ones(1,nphi), 1, nphi*ntheta);
        %
        % sangle = (reshape(stheta(:)*ones(1,nphi), 1, nphi*ntheta));
        angle(1, :) = stheta.*cphi;
        angle(2, :) = stheta.*sphi;
        angle(3, :) = ctheta;
        sangle = ones(1,nangs);

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

        Hhpf=Hhpf.*0.5;

        H =  Hzeeman + Hhpf*MS(cMS);
        H1 = Hzeeman + Hhpf*MS(cMS+1);

        w = 2 * sqrt(H(1,:).*H(1,:) + H(2,:).*conj(H(2,:)));
        w1 = 2 * sqrt(H1(1,:).*H1(1,:) + H1(2,:).*conj(H1(2,:)));
        % w = abs(2 * (H(1,:) + H(2,:)));
        % w1 = abs(2 * (H1(1,:) + H1(2,:)));



        % if Ex.ExciteWidth > 0
        %     gcenter = fld2g(Ex.Field*1E-3, Ex.mwFreq*1E9);
        %     gw      = gcenter*Ex.ExciteWidth/28/Ex.Field /sqrt(2*log(2));
        %     ak = exp(-2*((geff-gcenter)./gw).^2);
        % else
        %     ak = ones(1,ntheta*nphi);
        % end

        as = safeget(Ex, 'asym', 0);

        Spec = binning(Spec, Freqs, w, sangle.*ak.'*(1-as));
        Spec = binning(Spec, Freqs, w1, sangle.*ak.'*(1+as));

        
    end
end


if safeget(Sy, 'lwEndor', 0)
    mid = round(Ex.nPoints/2)+1;
    Line = lshape(Freqs,Freqs(mid),Sy.lwEndor);
    tdDecay = ifft(fftshift(Line*stepRF));
    td = ifft(Spec,[],2).*tdDecay;
    Spec = real(fft(td,[],2))*Ex.nPoints;
end

hole = safeget(Ex, 'Hole', 0);
if hole > 0
    etype = safeget(Ex, 'Scheme', 'Davies');
    if strcmp(etype,'Davies')
        profile = lshape(Freqs, omega_l, hole, 0, 0);
        profile = 1 - renorm(profile);
    else
        profile = sin(2*pi*(Freqs-omega_l)/hole).^2;
    end
    Spec = Spec.*profile;
end

hm = safeget(Ex, 'Harmonic', 0);
if hm > 0,
    Spec = pseumod(Freqs, Spec, Sy.lwEndor/10, hm);
end

% tm = clock - tm;
% disp(sprintf('%f min %f sec', tm(5),tm(6)))

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% calculate orientation selection and produce a set of angles
function varargout = Local_orisel(Sys, Exp, Opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
% function [Bx,Sa, Sb] = Local_orisel(Sys, xfreq, B0, Ex ,gt,D,E,krid, Blo,Bhi, npoint)
% Sys - structure
%   Sys.gn, Sys.I, Sys.nNucs, Sys.A [MHz], Sys.Apa [degree]
% xfreq - MW frequency, [GHz]
% D, E - zero-field spliting parameters either [1e-4*cm-1] or in [MHz]
% Opt.ZFunits -> D and E units: 'rcm' (default) or 'MHz'

%#########################################################################
%%%%%%%%%%%%%%%% get most of the parameters %%%%%%%%%%%%%%%%%%%%%%%%%
% gt - principles values of the g-tensor
% krid - nKnots
% Blo - lowest "recorded" field
% Bhi - highest field
% npoints - number of Field points
xfreq = Exp.mwFreq; %[GHz]

Ex = safeget(Exp, 'ExciteWidth', 10000); % [MHz] 
gt = safeget(Sys, 'g', [1, 1, 1]*2.0023);
DD = safeget(Sys, 'D', 0);

if Sys.S == 1/2, DD = [0, 0]; Sys.Type = 'norm'; end

if length(DD)==2
    D = DD(1);
    E = DD(2);    
else
    D = DD;
    E = safeget(Sys, 'E', 0);
end

krid = Opt.nKnots;

%#########################################################################
%############### generate a grid of orientations #########################
% 0 - rectangular, 1 - spiral (Ed, MR2), 2 - sphgrid (Stoll, Easyspin)
%%% NOTE: angles are arranged in rows !!!!!
gridtype = 1; 
if gridtype ==0
    ttthetas = linspace(0, 90, krid).'*pi/180;
    ttphis   = linspace(0, 90, krid).'*pi/180;
    thetas = reshape(ttthetas(ones(krid,1), :) , krid^2, 1);
    phis = reshape(ttphis(ones(krid,1), :).' , krid^2, 1);
    Weights = sin(thetas);
elseif gridtype ==1
    [phis, thetas]=Local_grid(krid);
%     phis = phis.';
%     thetas = thetas.';
    Weights = phis*0+1;
else
    [phis,thetas,Weights] = sphgrid('Ci',krid);
%%% NOTE: angles are arranged in rows !!!!!
% so....
    phis = phis.';
    thetas = thetas.';
    Weights = Weights.';
end

    
%#########################################################################
%############### get an EPR spectrum using a first order approach ########
sthetas2 = sin(thetas).^2;
cthetas2 = cos(thetas).^2;
sphis2 = sin(phis).^2;
cphis2 = cos(phis).^2;

%%%%%%%% geff^2 = n*g*g*n+ (n - unit vector along B0)
geff = sqrt(sthetas2.*cphis2*gt(1)^2  + sthetas2.*sphis2*gt(2)^2 + cthetas2*gt(3)^2);

% plank = 6.6261e-034; %Joule/Hz
% bohrmagn = 9.2740e-024; % Joule/Tesla.
% pbratio =  7.1448e-011;  % Tesla/Hz
pbratio =  7.1448e-02;  % mT/MHz
Bg = xfreq*1e3*pbratio./geff; % resonance without ZFS [mT]

%%%%%%%% ZFS 
%%%%%%%% wzfs = 3*n+*gt*Dt*gt+*n/2geff^2
%%%%%%%%%%%%%%%%%%%%%% CHECK %%%%%%%%%%%%%VVVV

rcmbohr = 4.66863 ;    % (bmagn/placnk/clight)   [1e-4cm-1/mT]
str = safeget(Opt, 'ZFunits', 'rcm');
if strcmp(str, 'rcm')
    convconst = rcmbohr; % 1/[1e-4*cm-1] -> 1/mT 
elseif strcmp(str, 'MHz')
    convconst = 1/pbratio; % 1/MHZ ->1/mT
else 
    error('Wrong ''Opt.ZFunits'' parameter (must be ''rcm'' for [1e-4cm-1] or ''MHz'' for [MHz])');
end
%%%%%%%%%%%%%%%%%%%% AAA

Ee = sthetas2.*(cphis2*gt(1)^2-sphis2*gt(2)^2)*E;
De = (-1/3*geff.^2 + gt(3)^2*cthetas2)*D;

gb = (geff.^3).*convconst; 
Dz = 1.5*(De+Ee)./gb; %  [mT]


Sa = Sys.S;
nSysDim = 2*Sa+1;

MS = [-Sa:1:Sa];

for cMS = 1:(nSysDim-1)

%     Ba = Bg - Dz;
%     Bb = Bg + Dz;
        dMS = MS(cMS) + MS(cMS+1);
    Br(:, cMS) = Bg + dMS*Dz;
    
    if strcmp(safeget(Sys, 'Type', 'ISC'), 'ISC')
        Wei(:, cMS) =  -dMS*Weights.*Dz/max(abs(Dz));
%         WeiB = -Weights.*Dz/max(abs(Dz));
    elseif strcmp(safeget(Sys, 'Type', 'RP'), 'RP')
        Wei(:, cMS) = -dMS*Weights;
%         WeiA = +Weights;
%         WeiB = -Weights;
    else
        Wei(:, cMS) = Weights;
%         WeiA = +Weights;
%         WeiB = +Weights;
    end
    %%%%%% alsi:: for my self :
    %%% planck*nu = g*bmagn*H -> Nu(exwidth)/Nu(mw) = H(exwidth)/H0
    if safeget(Opt, 'EPR', 0)
        B0 = mean(Bg);
        ExWidth = Ex/xfreq/1e3*B0; % [MHz] -> [mT]
    else
        B0 = Exp.Field; % [mT]
        ExWidth = Ex/xfreq/1e3*B0; % [MHz] -> [mT]
    end

    % include HF to the EPR spectrum (first order approach)
    useHFsel = safeget(Opt, 'HForisel', 0); % by default include this procedure to the calculation

    if useHFsel
        ctheta = cos(thetas);
        stheta = sin(thetas);
        cphi = cos(phis);
        sphi = sin(phis);
        angle(:, 1)=stheta.*cphi;
        angle(:, 2)=stheta.*sphi;
        angle(:, 3) = ctheta;
        nangles = length(angle);
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
                %%%% whf = nRGA
                LA(:, 1) = angle(:, 1) * At(1, 1) + angle(:, 2) * At(2, 1) + angle(:, 3) * At(3, 1);
                LA(:, 2) = angle(:, 1) * At(1, 2) + angle(:, 2) * At(2, 2) + angle(:, 3) * At(3, 2);
                LA(:, 3) = angle(:, 1) * At(1, 3) + angle(:, 2) * At(2, 3) + angle(:, 3) * At(3, 3);
                projA(:, k) = sqrt(LA(:, 1).^2 + LA(:, 2).^2 + LA(:, 3).^2)*1e6*planck/bmagn./geff*1e3; % [MHz] -> mT
                mult(k) = 2*Sys.I(ii)+1;

            end
        end
        nTrans = prod(mult);
        if length(mult)==1,
            MIs(:, 1) = (-(mult(1)-1)/2 : 1 : (mult(1)-1)/2 )';
        else
            prod2 = prod(mult(2:end));
            tA = (-(mult(1)-1)/2 : 1 : (mult(1)-1)/2)';
            MIs(:, 1) = kron(tA, ones(prod2, 1));

            for ii = 2:(length(mult)-1)
                tA = (-(mult(ii)-1)/2 : 1 : (mult(ii)-1)/2)';
                prod1 = prod(mult(1:(ii-1)));
                ttA = kron(ones(prod1, 1), tA);
                prod2 = prod(mult((ii+1):end));
                MIs(:, ii) = kron(ttA, ones(prod2, 1)); %  ???????????????????
            end
            tA = (-(mult(end)-1)/2 : 1 : (mult(end)-1)/2)';
            prod2 = prod(mult(1:(end-1)));
            MIs(:, length(mult)) = kron(ones(prod2, 1), tA);
        end
    end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if safeget(Opt, 'ShowEPRsel', 0)
% 
%         Blo = 0.99*min([min(Ba), min(Bb)]);%, B0-ExWidth
%         Bhi = 1.01*max([max(Ba), max(Bb)]);%, B0+ExWidth
%         npoint = 1000;
%         Bx = linspace(Blo,Bhi,1000).';
%         Sa = zeros(npoint, 1); Sb = Sa; S=Sa;
%     end

    %%%%%%%%%%%% just a trick to get EPR spectrum instead of ENDOR
    if safeget(Opt, 'EPR', 0)
        Blo = Exp.Range(1);%, B0-ExWidth
        Bhi = Exp.Range(2);%, B0+ExWidth
        npoint = Exp.nPoints;
        Bx = linspace(Blo,Bhi,npoint).';
        Sv(:, cMS) = zeros(npoint, 1);
        Ss(:, cMS) = zeros(npoint, 1);
        Opt.ShowEPRsel = 1;
    end

    if ~useHFsel
        ak(:, cMS) = exp(-2*((Br(:, cMS)-B0)./ExWidth).^2); % all in mT
%         akb = exp(-2*((Bb-B0)./ExWidth).^2); % all in mT
        if safeget(Opt, 'ShowEPRsel', 0)
            Sv(:, cMS) = binning(Sv(:, cMS), Bx, Br(:, cMS), Wei(:, cMS).*ak(:, cMS));
            Ss(:, cMS) = binning(Ss(:, cMS), Bx, Br(:, cMS), Wei(:, cMS));
%             Sb = binning(Sb, Bx, Bb, WeiB.*akb);
%             S = binning(S, Bx, Ba, WeiA);
%             S = binning(S, Bx, Bb, WeiB);
        end
    else
        % dE (epr) = nu_S + Sum(a*MIi)  (dMs = 1);
        % get all the projections MI1 ... MIn
        
        ak(:, cMS) = zeros(nangles, 1);
%         akb = zeros(nangles, 1);
        for ii = 1:nTrans
            ffre = zeros(nangles, 1);
            for jj = 1:length(mult)
                ffre = ffre + projA(:, jj)*MIs(ii, jj);
            end
            ak(:, cMS) = ak(:, cMS)+exp(-2*((Br(:, cMS)+ffre-B0)./ExWidth).^2); % all in mT
%             akb = akb+exp(-2*((Bb+ffre-B0)./ExWidth).^2); % all in mT
            if safeget(Opt, 'ShowEPRsel', 0)
                Sv(:, cMS) = binning(Sv(:, cMS), Bx, Br(:, cMS)+ffre, Wei(:, cMS).*ak(:, cMS));
                Ss(:, cMS) = binning(Ss(:, cMS), Bx, Br(:, cMS)+ffre, Wei(:, cMS));
%                 Sb = binning(Sb, Bx, Bb+ffre, WeiB.*akb);
%                 S = binning(S, Bx, Ba+ffre, WeiA);
%                 S = binning(S, Bx, Bb+ffre, WeiB);
            end
        end

    end



end
if safeget(Opt, 'ShowEPRsel', 0)
    if safeget(Opt, 'EPR', 0)
        step = Bx(end) - Bx(end-1);
        lw = safeget(Sys, 'lw', 1);
        wd = 0.5*lw/step; nn = round(4*wd);
        %%%%%% gaussian
        nrm = sqrt(2*pi)*wd;
        shape = exp(-0.5.*[-nn:nn].^2/wd^2)/nrm;
        
        S = sum(Ss, 2);
        
        S = local_convox(S.',shape);
        hm = safeget(Exp, 'Harmonic', 0);
        if hm > 0,
            S = kv_pseumod(Bx, S, Sys.lw, hm);
        end
        varargout{1} = Bx;
        varargout{2} = S;
        return;
        %%%%%%%%% out [X, Spec]
    else
        figure(Opt.ShowEPRsel); clf;
        step = Bx(end) - Bx(end-1);
        lw = safeget(Sys, 'lw', ExWidth/5);
        wd = 0.5*lw/step; nn = round(4*wd);
        %%%%%%% gaussian
        nrm = sqrt(2*pi)*wd;
        shape = exp(-0.5.*[-nn:nn].^2/wd^2)/nrm;
        
            S = sum(Ss, 2);
            
        S = local_convox(S.',shape);
        plot(Bx, S.', 'Color', [0, 0, 0]); hold on;
        
        for kk = 1:size(Sv, 2)
            Sv(:, kk) = local_convox(Sv(:, kk).',shape);
            plot(Bx, Sv(:, kk).', 'Color', [1, 0, 0]);
        end
%         Sb = local_convox(Sb.',shape);
%         plot(Bx, Sb.', 'Color', [0, 0, 1]);
    end
end

    w = ak.*Wei;
%     wb = akb.*WeiB;
    varargout{1} = phis;
    varargout{2} = thetas;
    varargout{3} = w;

    %%%%%%%%% out [phis, thetas, wa, wb]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phi, theta]=Local_grid(nKnots)
ncall = 0;
nselect= 0;

krid = nKnots;

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

    theta(end+1, 1) = thetap*pi/180; % rad
    phi(end+1, 1)   = phip*pi/180;   % rad
end
ww  = phi*0 + 1;


%============================================
% CONVOX
%
% convolute array x with window h
% 
% USAGE: y = convox(x,h)
%      
function y = local_convox(x,h)
n = length(x) ;
m = length(h) ; 
m2= (m-1)/2;
nn = n+2*m2;
hh=fliplr(h) ;
y = zeros(1,nn);
xx = [x(1)*ones(1,m2) x x(n)*ones(1,m2)];
%
for i = 1:nn-m+1
  y(i+m2) = xx(i:i+m-1)*hh' ;
end
y = y(m2+1:nn-m2);
return;

