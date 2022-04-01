function [Field, Spec] = EPRpolarized(Sys, Exp, Opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
% function [Field, Spec] = EPRpolarized(Sys, Exp, Opt)
% Sys - structure
%   Sys.gn, Sys.I, Sys.nNucs, Sys.A [MHz], Sys.Apa [degree]
% xfreq - MW frequency, [GHz]
% D, E - zero-field spliting parameters [1e-4*cm-1] or [MHz]
%
%#########################################################################
%%%%%%%%%%%%%%%% get most of the parameters %%%%%%%%%%%%%%%%%%%%%%%%%

xfreq = safeget(Exp, 'mwFreq', 9.3); %[GHz]
Brange = safeget(Exp, 'Range', [300, 400]); % [mT]
% B0 = Exp.Field; % [mT]
Ex = safeget(Exp, 'ExciteWidth', 0); % [MHz] 
gt = safeget(Sys, 'g', [1, 1, 1]*2.0023);
D = safeget(Sys, 'D', 0);
E = safeget(Sys, 'E', 0);

krid = safeget(Opt, 'nKnots', 50);

%#########################################################################
%############### generate a grid of orientations #########################
% 0 - ractangular, 1 - spiral (Ed, MR2), 2 - sphgrid (Stoll, Easyspin)
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
pbratio =  7.144773e-02;  %planck/bmagn*1e9 mT/MHz
rcmbohr = 4.668645 ;    % (bmagn/placnk/clight)   [1e-4cm-1/mT]

Bg = xfreq*1e3*pbratio./geff; % resonance without ZFS [mT]

%%%%%%%% ZFS 
%%%%%%%% wzfs = 3*n+*gt*Dt*gt+*n/2geff^2
%%%%%%%%%%%%%%%%%%%%%% CHECK %%%%%%%%%%%%%VVVV
str = safeget(Opt, 'ZFunits', 'rcm');
if strcmp(str, 'rcm')
    convconst = rcmbohr; % 1/[1e-4*cm-1] -> 1/mT 
elseif strcmp(str, 'MHz')
    convconst = 1/pbratio; % 1/MHZ ->1/mT
else 
    error('Wrong ''Opt.ZFunits'' parameter (must be ''rcm'' which means [1e-4cm-1] or ''MHz'' - [MHz])');
end

Ee = sthetas2.*(cphis2*gt(1)^2-sphis2*gt(2)^2)*E;
De = (-1/3*geff.^2 + gt(3)^2*cthetas2)*D;

gb = (geff.^3).*convconst; 
Dz = 1.5*(De+Ee)./gb; %  [mT]

Ba = Bg - Dz;
Bb = Bg + Dz;

%%%%%% alsi:: for my self :
%%% planck*nu = g*bmagn*H -> Nu(exwidth)/Nu(mw) = H(exwidth)/H0
% ExWidth = Ex/xfreq/1e3*B0; % [MHz] -> [mT]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% under construction !!!!!!!!!!!!!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%
% include HF to the EPR spectrum (first order approach)
useHFsel = safeget(Opt, 'HForisel', 0); % by default not include this procedure to the calculation

if useHFsel & isfield(Sys, 'A')
    Sys.nNucs = safeget(Sys, 'nNucs', size(Sys.A, 1));
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
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Blo = Brange(1); %0.99*min([min(Ba), min(Bb)]);%, B0-ExWidth
    Bhi = Brange(2); %1.01*max([max(Ba), max(Bb)]);%, B0+ExWidth
    npoint = safeget(Exp, 'nPoints', 1000);
    Bx = linspace(Blo,Bhi,npoint).';
    S = zeros(npoint, 1);% Sb = Sa; S=Sa;

if ~useHFsel
%     aka = exp(-2*((Ba-B0)./ExWidth).^2); % all in mT
%     akb = exp(-2*((Bb-B0)./ExWidth).^2); % all in mT
%     Sa = binning(Sa, Bx, Ba, Weights.*aka);
%     Sb = binning(Sb, Bx, Bb, -Weights.*akb);
    S = binning(S, Bx, Ba, Weights);
    S = binning(S, Bx, Bb, Weights); % !!!!!
else
    % dE (epr) = nu_S + Sum(a*MIi)  (dMs = 1);
    % get all the projections MI1 ... MIn
    nTrans = prod(mult);
    if length(mult)==1,
        MIs(:, 1) = (-(mult(1)-1)/2 : 1 : (mult(1)-1)/2 )';
    else
        prod2 = prod(mult(2:end));
        tA = (-(mult(1)-1)/2 : 1 : (mult(1)-1)/2)';
        MIs(:, 1) = kron(tA, ones(prod2, 1));

        for ii = 2:(length(mult)-1)
            tA = (-(mult(ii)-1)/2 : 1 : (mult(ii)-1)/2)';
            prod1 = prod(mult(1:ii));
            ttA = kron(ones(prod2, 1), tA);
            prod2 = prod(mult((ii+1):end));
            MIs(:, ii) = kron(ttA, ones(prod2, 1));
        end
        tA = (-(mult(end)-1)/2 : 1 : (mult(end)-1)/2)';
        prod2 = prod(mult(1:(end-1)));
        MIs(:, length(mult)) = kron(ones(prod2, 1), tA);
    end
%     aka = zeros(nangles, 1);
%     akb = zeros(nangles, 1);
    for ii = 1:nTrans
        ffre = zeros(nangles, 1);
        for jj = 1:length(mult)
            ffre = ffre + projA(:, jj)*MIs(ii, jj);
        end
%         aka = aka+exp(-2*((Ba+ffre-B0)./ExWidth).^2); % all in mT
%         akb = akb+exp(-2*((Bb+ffre-B0)./ExWidth).^2); % all in mT
        
%             Sa = binning(Sa, Bx, Ba+ffre, Weights.*aka);
%             Sb = binning(Sb, Bx, Bb+ffre, -Weights.*akb);
            S = binning(S, Bx, Ba+ffre, Weights);
            S = binning(S, Bx, Bb+ffre, -Weights);
            
    end

end

%     figure(Opt.ShowEPRsel); clf;
    step = Bx(end) - Bx(end-1);
    lw = safeget(Sys, 'lw', step*5);
    wd = 0.5*lw/step; nn = round(4*wd);
    %%%%%%% gaussian
    nrm = sqrt(2*pi)*wd;
    shape = exp(-0.5.*[-nn:nn].^2/wd^2)/nrm;
    Spec = local_convox(S.',shape);
%     Spec = S;
Field = Bx;
%     
%     S = local_convox(S.',shape);
%     plot(Bx, S.', 'Color', [0, 0, 0]); hold on;
%     Sa = local_convox(Sa.',shape);
%     plot(Bx, Sa.', 'Color', [1, 0, 0]);
%     Sb = local_convox(Sb.',shape);
%     plot(Bx, Sb.', 'Color', [0, 0, 1]); 


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

