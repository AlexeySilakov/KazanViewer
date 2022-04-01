%==========================================================================
function [Freqs, Spec] = mordor(Sys,Exp,Opt)
%==========================================================================
if nargin < 2
    error('Usage: mordor(Sys,Exp,Opt).');
end

% tm = clock;

if isfield(Sys, 'Nucs'), [Sys.I,Sys.gn] = nucdata(Sys.Nucs); end
if isfield(Exp, 'MaxFreq'), Exp.Range = [-Exp.MaxFreq, Exp.MaxFreq]; end
if ~isfield(Exp, 'HTA'), Exp.HTA = 1; end % High turn angle ... "1" means HTA is a normal pi pulse
nucgfactor=Sys.gn;
Lplanck = 6.6261e-034; % J/s
Lnmagn = 5.0508e-027;  % J/T
Lbmagn = 9.27401e-024; % J/T
mTtoMHz = Lplanck/Lbmagn*1e9;

Exp.Range=safeget(Exp, 'Range', [0 100]);

[phi, theta, weights] = Local_orisel(Sys, Exp, Opt);

if size(Sys.A, 1)==3
    A = Sys.A;
else
    R = erot(Sys.Apa);
    A = R * diag(Sys.A)*R.';
    if isfield(Sys, 'M')
        if size(Sys.M, 1)~=3
            M = diag(Sys.M);
        else
            M = Sys.M;
        end
        A = A*M; %%%%%%%%%% for this to work A(1, 3) is for SzIx
    end
end

R = erot(Sys.Qpa);
Q = R * diag(Sys.Q)*R.';

ctheta = cos(theta);
stheta = sin(theta);

cphi = cos(phi);
sphi = sin(phi);

angle(1,:)= stheta(:).*cphi(:);
angle(2,:)= stheta(:).*sphi(:);
angle(3,:)= ctheta(:);


Bo = Exp.Field/mTtoMHz;
Bon = Exp.Field*Lnmagn/Lplanck*1e-9;

Gx=angle(1,:)*Sys.g(1);
Gy=angle(2,:)*Sys.g(2);
Gz=angle(3,:)*Sys.g(3);

geff=sqrt(Gx.*Gx+Gy.*Gy+Gz.*Gz);
       
% G = g.*reshape([a_sthetacphi(:) a_sthetasphi(:) a_ctheta(:)], 3, ntheta*nphi);

% Axyz(1,:)=(Gx*A(1,1) + Gy*A(1,2) + Gz*A(1,3))./geff;
% Axyz(2,:)=(Gx*A(2,1) + Gy*A(2,2) + Gz*A(2,3))./geff;
% Axyz(3,:)=(Gx*A(3,1) + Gy*A(3,2) + Gz*A(3,3))./geff;

Sx = local_sop([Sys.S, Sys.I], 'xe');
Sy = local_sop([Sys.S, Sys.I], 'ye');
Sz = local_sop([Sys.S, Sys.I], 'ze');

Ix = local_sop([Sys.S, Sys.I], 'ex');
Iy = local_sop([Sys.S, Sys.I], 'ey');
Iz = local_sop([Sys.S, Sys.I], 'ez');

HF = Sx*A(1, 1)*Ix + Sx*A(2, 1)*Iy + Sx*A(3, 1)*Iz +...
     Sy*A(1, 2)*Ix + Sy*A(2, 2)*Iy + Sy*A(3, 2)*Iz +...
     Sz*A(1, 3)*Ix + Sz*A(2, 3)*Iy + Sz*A(3, 3)*Iz ;

 
HFIx = Ix*A(1, 1) +Iy*A(2, 1) +Iz*A(3, 1);
HFIy = Ix*A(1, 2) +Iy*A(2, 2) +Iz*A(3, 2);
HFIz = Ix*A(1, 3) +Iy*A(2, 3) +Iz*A(3, 3);
 
 if Sys.I>1/2
     HF1 = Ix*Q(1, 1)*Ix + Ix*Q(1, 2)*Iy + Ix*Q(1, 3)*Iz +...
          Iy*Q(2, 1)*Ix + Iy*Q(2, 2)*Iy + Iy*Q(2, 3)*Iz +...
          Iz*Q(3, 1)*Ix + Iz*Q(3, 2)*Iy + Iz*Q(3, 3)*Iz  ;
 else
     HF1 = Ix*0;
 end
nS = ones(length(Sx), 1);
nxI = 2*Sys.I+1; 
nxS = 2*Sys.S+1;

% %%%%%%%% idea: to get transition intensity, we need to know how strong is
% %%%%%%%% the eldor pulse ... Do not care about B1- exactly, as we can talk
% %%%%%%%% in terms of how many times we overrotate an allowed transiton. By
% %%%%%%%% knowing how long is the ordinary pi pulse, and knowing the
% %%%%%%%% attenuation between eldor and MW channels,we can just get that
% %%%%%%%% info....
% attEl = 0; %dB
% lenEl = 3000*10^(-attEl/10); %ns
% lenPi = 16; % ns
% length_eldor = 5; % in comparisoin within terms of allowed transition

Spec = zeros(Exp.nPoints, 1);
Freqs = linspace(Exp.Range(1), Exp.Range(2), Exp.nPoints)';
for ii = 1:length(phi)
    if weights(ii)<1e-3, continue; end
    lx = sin(theta(ii))*cos(phi(ii));
    ly = sin(theta(ii))*sin(phi(ii));
    lz = cos(theta(ii));
    lxdet = sin(theta(ii)+pi/2)*cos(phi(ii));
    lydet = sin(theta(ii)+pi/2)*sin(phi(ii));
    lzdet = cos(theta(ii)+pi/2);
    Ham = Sys.g(1)*Bo*lx*Sx + Sys.g(2)*Bo*ly*Sy + Sys.g(3)*Bo*lz*Sz +...
          Sys.gn*Bon*lx*Ix + Sys.gn*Bon*ly*Iy + Sys.gn*Bon*lz*Iz + HF+HF1;
    Det= lxdet*Sx +lydet*Sy + lzdet*Sz;
%     HH = HFIx*lx*Sys.g(1)+HFIy*ly*Sys.g(2)+HFIz*lz*Sys.g(3);
%     Ham = geff(ii)*Bo*Sz + Sys.gn*(Iz*lz+Ix*lx+Iy*ly)*Bon+HF1+ Sz*HH/geff(ii);
%     Det  =Sx;
    [V, D] = eig(Ham);
    
    M = sqrt(abs((V'*Det*V).^2))/2;
%     Meldor = 1-cos(2*pi*M*Exp.HTA).^2;
%     [ci, cj] = find(M>0.4); 
    
    eVal = diag(D);
    dd = real(eVal(:, nS) - eVal(:, nS)');
    
%     idx1 = find(ci<=nxI);
    vv = [];
    aa = [];
    
    for jj = 1:nxI
        for cc = [1:nxI]+nxI
            for kk = [1:nxI]+nxI
                if kk == cc, continue; end; % we hit the observer transition
                vv = [vv, dd( jj, kk ) - dd( jj, cc ) ];
                aa = [aa, (1-cos(pi*M( jj, cc )*Exp.HTA))*M( jj, kk )]; % probablilty to observe +probability to invert
%                 aa = [aa, M( jj, kk )*Meldor( jj, cc)]; % probablilty to observe +probability to invert
            
            end            
        end
    end
    for jj = [1:nxI]+nxI
        for cc = 1:nxI
            for kk = 1:nxI
                if kk == cc, continue; end; % we hit the observer transition
                vv = [vv, dd( jj, kk ) - dd( jj, cc ) ];
                aa = [aa, (1-cos(pi*M( jj, cc )*Exp.HTA))*M( jj, kk )];%
%                 aa = [aa, M( jj, kk )*Meldor( jj, cc)]; % probablilty to observe +probability to invert
            
            end            
        end
    end
%     for jj = 1:length(idx1)
%         for kk = 
%             if kk == cj(idx1(jj)), continue; end; % we hit allowed transition
%             vv = [vv, dd( ci(idx1(jj)), kk ) - dd( ci(idx1(jj)), cj(idx1(jj)) ) ];
%             aa = [aa, M( ci(idx1(jj)), kk )];
%         end
%     end
%     idx1 = find(ci>nxI);
%     for jj = 1:length(idx1)
%         for kk = [1:nxI]
%             if kk == cj(idx1(jj)), continue; end % we hit allowed transition
%             vv = [vv, dd( ci(idx1(jj)), kk ) - dd( ci(idx1(jj)), cj(idx1(jj)) ) ];
%             aa = [aa, M( ci(idx1(jj)), kk )];
%         end
%     end
    Spec = binning(Spec, Freqs, (vv), abs(aa)*weights(ii));
    Spec = binning(Spec, Freqs, -(vv), abs(aa)*weights(ii));
    
end
lwENDOR = safeget(Sys, 'lwEndor', 0);
if lwENDOR
    
    step = Freqs(2)-Freqs(1);
    
    wd = 0.5*lwENDOR/step; nn = round(4*wd);
    %%%%%%% gaussian
    nrm = sqrt(2*pi)*wd;
    aa = fftshift(gaussian(Freqs, (max(Freqs)+min(Freqs))/2, lwENDOR, safeget(Exp, 'Harmonic', 0)));
    iaa = ifft(aa/max(aa));
    iSpec = ifft(Spec);
    Spec = real(fft(iaa.*iSpec));
end
% if Sys.lwEndor
%     mid = round(Exp.nPoints/2)+1;
%     stepRF = Freqs(2)-Freqs(1);
%     Line = lshape(Freqs,Freqs(mid),Sys.lwEndor);
%     tdDecay = ifft(fftshift(Line*stepRF));
%     td = ifft(Spec,[],2).*tdDecay;
%     Spec = abs(fft(td,[],2))*Exp.nPoints;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% calculate orientation selection and produce a set of angles
function [phi, theta, weights] = Local_orisel(Sys, Exp, Opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
if isfield(Opt, 'OriSelFile') % load from an exsisting file intead of generating the grid
%     "phi", "theta" and "weights" 
    load(Opt.OriSelFile, '-mat'); % if there is no file, load will do error by itself, no need to protect that
    if exist('phi')&exist('theta')&exist('weights')
        if size(phi, 1)==1
            phi = phi.';
            theta = theta.';
            weights = weights.';
        end
        return;
    else
        filestruct = load(Opt.OriSelFile, '-mat');
        disp(filestruct);
        error(['file ', Opt.OriSelFile, ' does not contains valid grid parameters']);
    end
end
if isfield(Opt, 'OriSelMatrix') % load from an exsisting file intead of generating the grid
%     "phi", "theta" and "weights" 
%     load(Opt.OriSelFile, '-mat'); % if there is no file, load will do error by itself, no need to protect that
    if size(Opt.OriSelMatrix, 2)==3
%         if size(phi, 1)==1
            phi = Opt.OriSelMatrix(:, 1);
            theta = Opt.OriSelMatrix(:, 2);
            weights = Opt.OriSelMatrix(:, 3);
%         end
        return;
    else
        disp(size(Opt.OriSelMatrix));
        error(['Supplied OriSelMatrix does not contains valid grid parameters, [phi(:), theta(:), weights(:)] ']);
    end
end
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
%if isfield(Sys, 'Nucs'), [Sys.I,Sys.gn] = nucdata(Sys.Nucs); end % EasySpin function



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
    Sys.nNucs = safeget(Sys, 'nNucs', Sys.I*0+1);
    projA = zeros(nangles, sum(Sys.nNucs));
    k = 0;
    for ii = 1:length(Sys.nNucs)
        for jj = 1:Sys.nNucs(ii)
            k = k+1;
            if sum(Sys.Apa(ii, :))
                R = erot(Sys.Apa(ii, :));
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

