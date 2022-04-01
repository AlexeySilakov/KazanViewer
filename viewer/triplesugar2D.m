% function [Freq1, Freq2, Spec] = triplesugar2D(Sys,Exp,Opt)
% Attention! New Hamiltonian is used for any pair of nuclei.
% difference from EasySpin notation:
%  Sys.Amp= array of weights for different nuclei 
%   the weight of triple spectrum is calculated from 
%   the multiplication of weights of corresponding nuclei
%==========================================================================
function [Freqs, Freq2, Spec] = triplesugar2D(Sy,Ex,Si)
%==========================================================================

if nargin < 2
    error('Usage: triplesugar2D(Sys,Exp,Opt).');
end
% determine the number of spins
snumber=length(Sy.I);

if snumber > size(Sy.A,1) | snumber > size(Sy.Apa,1)
    error('size of I parameter has to be equal or smaller than first dimension of A and Apa!')
end

Spec = zeros(Ex.nPoints,Ex.nPoints);
Amp = safeget(Sy,'Amp',ones(1,snumber));
gn  = safeget(Sy,'gn',ones(1,snumber)*Sy.gn(1));

% create matrixes of angles
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
sangle = reshape(stheta(:)*ones(1,nphi), 1, nphi*ntheta);

stepRF = (Ex.Range(2)-Ex.Range(1))/(Ex.nPoints-1);
Freqs = [Ex.Range(1) : stepRF : Ex.Range(2)];
Freq2 = Freqs;

for sp1=1:(snumber-1)
    for sp2=(sp1+1):snumber
        R1 = erot(Sy.Apa(sp1,:));
        A1 = R1 * diag(Sy.A(sp1,:))*R1.';
        
        R2 = erot(Sy.Apa(sp2,:));
        A2 = R2 * diag(Sy.A(sp2,:))*R2.';
        
        G1x=angle(1,:)*Sy.g(sp1,1);
        G1y=angle(2,:)*Sy.g(sp1,2);
        G1z=angle(3,:)*Sy.g(sp1,3);
        
        G2x=angle(1,:)*Sy.g(sp2,1);
        G2y=angle(2,:)*Sy.g(sp2,2);
        G2z=angle(3,:)*Sy.g(sp2,3);
        
        geff1=sqrt(G1x.*G1x+G1y.*G1y+G1z.*G1z);
        geff2=sqrt(G2x.*G2x+G1y.*G2y+G2z.*G2z);
        
        omega_l = Ex.Field/planck /1E9 * gn(sp1) * nmagn;
        B = omega_l * angle;
        
        % G = g.*reshape([a_sthetacphi(:) a_sthetasphi(:) a_ctheta(:)], 3, ntheta*nphi);
        
        A1xyz(1,:)=(G1x*A1(1,1) + G1y*A1(1,2) + G1z*A1(1,3))./geff1;
        A1xyz(2,:)=(G1x*A1(2,1) + G1y*A1(2,2) + G1z*A1(2,3))./geff1;
        A1xyz(3,:)=(G1x*A1(3,1) + G1y*A1(3,2) + G1z*A1(3,3))./geff1;
        
        A2xyz(1,:)=(G2x*A2(1,1) + G2y*A2(1,2) + G2z*A2(1,3))./geff2;
        A2xyz(2,:)=(G2x*A2(2,1) + G2y*A2(2,2) + G2z*A2(2,3))./geff2;
        A2xyz(3,:)=(G2x*A2(3,1) + G2y*A2(3,2) + G2z*A2(3,3))./geff2;
        
        H1zeeman(1,:) = -B(3,:);
        H1zeeman(2,:) = -B(1,:) + j*B(2,:);
        
        H2zeeman(1,:) = -B(3,:);
        H2zeeman(2,:) = -B(1,:) + j*B(2,:);
        
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
        
        k1= floor( (w1 - Ex.Range(1))/stepRF);
        k2= floor( (w2 - Ex.Range(1))/stepRF);
        k1d= floor( (w1d - Ex.Range(1))/stepRF);
        k2d= floor( (w2d - Ex.Range(1))/stepRF);
        
        if Ex.ExciteWidth > 0
            % excitation profile
            gcenter = fld2g(Ex.Field*1E-3, Ex.mwFreq*1E9);
            gw      = gcenter*Ex.ExciteWidth/28/Ex.Field /sqrt(2*log(2));
            ak = exp(-2*((geff1-gcenter)./gw).^2).*exp(-2*((geff2-gcenter)./gw).^2);
        else
            ak = ones(1,ntheta*nphi);
        end
        % weight coefficients
        ak=ak*Amp(sp1)*Amp(sp2);

        bin2D(Spec, k1, k2, sangle.*ak);
        bin2D(Spec, k2, k1, sangle.*ak);
        
        bin2D(Spec, k1d, k2d, sangle.*ak);
        bin2D(Spec, k2d, k1d, sangle.*ak);
        
        % diagonal elements
       bin2D(Spec, k1, k1, sangle.*ak);
        bin2D(Spec, k2, k2, sangle.*ak);
        bin2D(Spec, k1d, k1d, sangle.*ak);
        bin2D(Spec, k2d, k2d, sangle.*ak);
  
        % linewidth
        mid = round(Ex.nPoints/2)+1;
        Line = lshape(Freqs,Freqs(mid),Sy.lwEndor,0, 1); %'lorentz'
        for jj=1:size(Freqs, 2); Line2D(jj, :)=Line.*Line(jj); end
        tdDecay = ifft2(fftshift(Line2D*stepRF));
        td = ifft2(Spec).*tdDecay;
        Spec = abs(fft2(td))*Ex.nPoints;
        
    end
end
return
