% function [varargout] = sugar(Sys,Exp,varargin)
% arbitrary number of nuclei

%==========================================================================
function [Freqs, Spec] = triplesugar1D(Si, Ex, Opt)
%==========================================================================
% Counting of the niclei
num_nucs = size(Si.A, 1);
CombNucs = nchoosek(1:num_nucs,2); % all possible combinations
Spec = zeros(Ex.nPoints,Ex.nPoints);
for ci = 1:size(CombNucs, 1)
    nuc1 = CombNucs(ci, 1);
    nuc2 = CombNucs(ci, 2);   
    [g1, g2]       = Local_getval(Si, 'g',   CombNucs(ci, :), num_nucs);
    [gNuc1, gNuc2] = Local_getval(Si, 'gn',  CombNucs(ci, :), num_nucs);
    [Apa1, Apa2]   = Local_getval(Si, 'Apa', CombNucs(ci, :), num_nucs);
    [A1, A2]       = Local_getval(Si, 'A',   CombNucs(ci, :), num_nucs);
    
    R1 = erot(Apa1);
    A1 = R1 * diag(A1)*R1.';
    
    R2 = erot(Apa2);
    A2 = R2 * diag(A2)*R2.';
    
    Opt.nKnots  = safeget(Opt, 'nKnots', 10);
    
    if isfield(Opt, 'ThetaRange')
        ntheta = Opt.nKnots*2;
        nphi = Opt.nKnots*4;
        PhiRange = safeget(Si, 'PhiRange', [0,2*pi]);
        ThetaRange = safeget(Si, 'ThetaRange', [0, pi]);
        steptheta = (ThetaRange(2) - ThetaRange(1))/(ntheta - 1);
        stepphi = (PhiRange(2) - PhiRange(1))/(nphi - 1);
        t_theta = ThetaRange(1):steptheta:ThetaRange(2);
        t_phi   = PhiRange(1):stepphi:PhiRange(2);
        theta = reshape(t_theta(ones(  nphi, 1), :)  , 1, nphi*ntheta);
        phi   = reshape(  t_phi(ones(ntheta, 1), :).', 1, nphi*ntheta);    
        sangle = sin(theta);
    else
        Opt.Symmetry = safeget(Opt, 'Symmetry', 'Ci');
        [phi,theta,sangle] = sphgrid(Opt.Symmetry, Opt.nKnots);
    end
    
    ctheta = cos(theta);
    stheta = sin(theta);
    
    cphi = cos(phi);
    sphi = sin(phi);
    
    angle(1,:)=stheta.*cphi;
    angle(2,:)=stheta.*sphi;
    angle(3,:)=ctheta;
    
    omega_l = Ex.Field/planck /1E9 * gNuc1 * nmagn;
    omega_2 = Ex.Field/planck /1E9 * gNuc2 * nmagn;
    
    B1 = omega_l * angle;
    B2 = omega_2 * angle;
    
    G1x=angle(1,:)*g1(1);
    G1y=angle(2,:)*g1(2);
    G1z=angle(3,:)*g1(3);
    
    G2x=angle(1,:)*g2(1);
    G2y=angle(2,:)*g2(2);
    G2z=angle(3,:)*g2(3);
    
    geff1=sqrt(G1x.*G1x+G1y.*G1y+G1z.*G1z);
    geff2=sqrt(G2x.*G2x+G1y.*G2y+G2z.*G2z);
    
    % G = g.*reshape([a_sthetacphi(:) a_sthetasphi(:) a_ctheta(:)], 3, ntheta*nphi);
    
    A1xyz(1,:)=(G1x*A1(1,1) + G1y*A1(1,2) + G1z*A1(1,3))./geff1;
    A1xyz(2,:)=(G1x*A1(2,1) + G1y*A1(2,2) + G1z*A1(2,3))./geff1;
    A1xyz(3,:)=(G1x*A1(3,1) + G1y*A1(3,2) + G1z*A1(3,3))./geff1;
    
    A2xyz(1,:)=(G2x*A2(1,1) + G2y*A2(1,2) + G2z*A2(1,3))./geff2;
    A2xyz(2,:)=(G2x*A2(2,1) + G2y*A2(2,2) + G2z*A2(2,3))./geff2;
    A2xyz(3,:)=(G2x*A2(3,1) + G2y*A2(3,2) + G2z*A2(3,3))./geff2;
    
    H1zeeman(1,:) = -B1(3,:);
    H1zeeman(2,:) = -B1(1,:) + j*B1(2,:);
    
    H2zeeman(1,:) = -B2(3,:);
    H2zeeman(2,:) = -B2(1,:) + j*B2(2,:);
    
    H1zeeman = H1zeeman.*0.5;
    H2zeeman = H2zeeman.*0.5;
    
    H1hpf(1,:)= A1xyz(3,:);
    H1hpf(2,:)= A1xyz(1,:) - j* A1xyz(2,:);
    
    H2hpf(1,:)= A2xyz(3,:);
    H2hpf(2,:)= A2xyz(1,:) - j* A2xyz(2,:);
    
    H1hpf=H1hpf.*0.25;
    H2hpf=H2hpf.*0.25;
    
    H1 =  H1zeeman + H1hpf;
    H1d = H1zeeman - H1hpf;
    
    H2 =  H2zeeman + H2hpf;
    H2d = H2zeeman - H2hpf;
    
    w1 = 2 * sqrt(H1(1,:).*H1(1,:) + H1(2,:).*conj(H1(2,:)));
    w1d = 2 * sqrt(H1d(1,:).*H1d(1,:) + H1d(2,:).*conj(H1d(2,:)));
    
    w2 = 2 * sqrt(H2(1,:).*H2(1,:) + H2(2,:).*conj(H2(2,:)));
    w2d = 2 * sqrt(H2d(1,:).*H2d(1,:) + H2d(2,:).*conj(H2d(2,:)));
    
    stepRF = (Ex.Range(2)-Ex.Range(1))/(Ex.nPoints-1);
    Freqs = [Ex.Range(1) : stepRF : Ex.Range(2)];
    
    k1= floor( (w1 - Ex.Range(1))/stepRF);
    k2= floor( (w2 - Ex.Range(1))/stepRF);
    k1d= floor( (w1d - Ex.Range(1))/stepRF);
    k2d= floor( (w2d - Ex.Range(1))/stepRF);
    
    if Ex.ExciteWidth > 0
        gcenter = fld2g(Ex.Field*1E-3, Ex.mwFreq*1E9);
        gw      = gcenter*Ex.ExciteWidth/28/Ex.Field /sqrt(2*log(2));
        ak = exp(-2*((geff1-gcenter)./gw).^2).*exp(-2*((geff2-gcenter)./gw).^2);
    else
        ak = ones(1,ntheta*nphi);
    end
    
    Spec = bin2D(Spec, k1, k2, sangle.*ak);
    Spec = bin2D(Spec, k2, k1, sangle.*ak);
    
    Spec = bin2D(Spec, k1d, k2d, sangle.*ak);
    Spec = bin2D(Spec, k2d, k1d, sangle.*ak);
    
    % diagonal elements
    Spec = bin2D(Spec, k1, k1, sangle.*ak);
    Spec = bin2D(Spec, k2, k2, sangle.*ak);
    Spec = bin2D(Spec, k1d, k1d, sangle.*ak);
    Spec = bin2D(Spec, k2d, k2d, sangle.*ak);
end
% linewidth
Line2D = gaga2(1:Ex.nPoints, (Ex.nPoints+1)/2,...
    1:Ex.nPoints, (Ex.nPoints+1)/2, [1, 1]*Si.lwEndor*Ex.nPoints/(Ex.Range(2)-Ex.Range(1)));

tdSpec = ifft2(fftshift(Spec));
tdDecay = ifft2(fftshift(Line2D)*stepRF);
td = tdSpec.*tdDecay/max(max(tdDecay));
Spec = fftshift(abs(fft2(td)))*Ex.nPoints;
if safeget(Opt, 'AddFig', 0)
    figure(safeget(Opt, 'AddFig', 0));
    plot(Freqs, max(Spec));
end
if Ex.TrFreq,
    pos = floor( (Ex.TrFreq - Ex.Range(1))/stepRF+0.5)+1; % "+1" because we start with index=1, not 0 like in C++
    ind = find(pos<size(Spec,2)& pos>0);
    if ~isempty(ind)
        Spec = Spec(:, pos(ind));
    else
        Spec = zeros(Ex.nPoints, 1);
    end
end
return

function varargout = gaga2(x, xo, y, yo, fwhh, varargin)
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
X = tx;
Y = ty;
Z = out;
if nargout > 1,
    varargout{1} = tx;
    varargout{2} = ty;
    varargout{3} = out;
elseif nargout == 1
    varargout{1} = out;
else
    disp('something is wrong');
end

function [Val1, Val2] = Local_getval(Str, Field, nums, max)
if isfield(Str, Field)
    Values = getfield(Str, Field);
else
    Values = 0;
end
if size(Values, 1)==max
    Val1 = Values(nums(1), :);
    Val2 = Values(nums(2), :);
elseif size(Values, 1)==1
    Val1 = Values(1, :);
    Val2 = Values(1, :);    
else
    error(['Parameter ',Str,'.',Field 'has wrong number of rows']);
end

    
