% EPR simulations with isotropic g factor, HPF
% and trigonal D. Second order perturbation theory.
% [xx,yy]=EPRgisoD(Sy,Ex,Op)
function [xx,yy]=EPRgisoD(Sy,Ex,Op)
tt=cputime;
%
S = safeget(Sy, 'S', 1);
I = safeget(Sy, 'I', 0);

% isotropic g factor
g=safeget(Sy,'g',2.0023*[1 1 1]);
g=sum(g)/length(g); 

% isotropic hpf
a=safeget(Sy,'A',0*[1 1 1]);
a=sum(a)/length(a); 

DD=safeget(Sy,'D',1000*[-1 -1 2]);
if sum(DD) > 10, error('Non traceless ZFS'); end
[D,idx] = max(abs(DD)); D = D*3/2;
eta = abs(diff(DD([1:3]~=idx)))/D/2;

freq = safeget(Ex, 'mwFreq', 34)*1E3;
nTheta = safeget(Ex,'nTheta',1001);
nPhi   = safeget(Ex,'nPhi',61);

cth    = 0:1/(nTheta-1):1;
cth2   = cth.^2; sth2   = 1 - cth2;
sth    = sqrt(sth2);
s2th   = 2*cth.*sth;
c2th   = cth2 - sth2;
c2phi  = cos(0:2*pi/(nPhi-1):2*pi);
s2phi  = sin(0:2*pi/(nPhi-1):2*pi);

A   = reshape(D*((3*cth2(:)-1)/2*ones(1,nPhi)+3*eta/2*sth2(:)*c2phi(:)'), 1, nPhi*nTheta);
Amp = ones(1,nPhi*nTheta);
Bp  = reshape(D/4*s2th(:)*(eta*c2phi(:)'-1), 1, nPhi*nTheta); 
pm  = reshape(D/4*i*2*eta*sth(:)*s2phi(:)', 1, nPhi*nTheta);
Bm  = Bp - pm; Bp  = Bp + pm;
Bpm = Bp.*Bm;
Cp  = reshape(D/4*(sth2(:)*ones(1,nPhi)+eta*(cth2(:)+1)*c2phi(:)'), 1, nPhi*nTheta); 
pm  = reshape(D/4*i*2*eta*cth(:)*s2phi(:)', 1, nPhi*nTheta);
Cm  = Cp - pm; Cp  = Cp + pm;
Cpm = Cp.*Cm;
SS1 = S*(S+1);

xx = linspace(Ex.Range(1),Ex.Range(2),Ex.nPoints);
yy = zeros(1,Ex.nPoints);
m = 0; m1 = 0; a = 0;
for M=-S:S-1
     Amp1 = SS1-M*(M+1);
%      for m = -I:I
%          for m = -I:I
            H = (freq - a*(m+m1)-A*(1+2*M))/g/13.9962;
            % second order corrections
            H = H - (2*Bpm*(8*((M+1)^3-M*M*M)+1-4*SS1)-...
                2*Cpm*(2*((M+1)^3-M*M*M)+1-2*SS1))/g/freq/13.9962;
            % hpf corrections
%             H = H - a*a/2*(I*(I+1)+(m+m1)*M*(M+1)-m1*m1-m*m)/g/freq/13.9962
            yy = binning(yy,xx,abs(H),Amp*Amp1);
%         end
%     end
end

yy = convspec(yy,xx(2)-xx(1),Sy.lw,safeget(Ex,'Harmonic',0));
disp(sprintf('Exec time: %5.1fs',cputime-tt))
