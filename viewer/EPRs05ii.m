% EPR simulation with s=1/2 and 2 nuclei with
% colinear HPf and quarupole tenzors.
% ATTENTION ! contains mistake(s) 
% [xx,yy]=ENDORs05ii(Sy,Ex,Op)
function [xx, yy] = EPRs05ii(Sy, Ex, Op)
tt=cputime;
%
S = safeget(Sy, 'S', 1/2);
if S~=1/2, error('Only S=1/2 supported by ENDORs05ii function'); end

if isfield(Sy, 'Nucs'), [Sy.I,Sy.gn] = nucdata(Sy.Nucs); end
II = safeget(Sy, 'I', 5/2);

% g factor
gg = safeget(Sy,'g',2.0023*[1 1 1]);

% a hpf
AA = safeget(Sy,'A',[1 1 1; 1 1 1]);

% q hpf
QQ = safeget(Sy,'Q',0*[-1 -1 2; -1 -1 2]/3);
if sum(sum(QQ)) > 0.001, error('Non-traceles Q'); end

freq = safeget(Ex, 'mwFreq', 34)*1E3;
B0 = safeget(Ex, 'Field', 1240);
lwENDOR    = safeget(Sy, 'lwEndor', .5);
lw         = safeget(Sy, 'lw', 5);
ExWidth = safeget(Ex, 'ExciteWidth', 5);

vi = -B0/1E9*Sy.gn*7622587;

% angles 
nTheta = safeget(Ex,'nTheta',300);
nPhi   = safeget(Ex,'nPhi',50);

on_s  = ones(1,nPhi*nTheta);
cth   = 0:1/(nTheta-1):1;
cth2  = cth.^2; sth2   = 1 - cth2;
cphi  = cos(0:pi/(nPhi-1):pi);
cphi2 = cphi.*cphi;
sphi2 = 1 - cphi2;
s2c2  = reshape(sth2'*cphi2,1,nTheta*nPhi);
s2s2  = reshape(sth2'*sphi2,1,nTheta*nPhi);
c2    = reshape(cth2(ones(1,nPhi),:)',1,nTheta*nPhi);
fAmp   = on_s([1,1],:);

g = sqrt(gg(1)*gg(1)*s2c2 + gg(2)*gg(2)*s2s2 + gg(3)*gg(3)*c2); 
a(1,:) = AA(1,1)*s2c2 + AA(1,2)*s2s2 + AA(1,3)*c2;
a(2,:) = AA(2,1)*s2c2 + AA(2,2)*s2s2 + AA(2,3)*c2;  
q(1,:) = QQ(1,1)*s2c2 + QQ(1,2)*s2s2 + QQ(1,3)*c2;
q(2,:) = QQ(2,1)*s2c2 + QQ(2,2)*s2s2 + QQ(2,3)*c2;  

nPoints  = safeget(Ex,'nPoints',1000);
xx = linspace(Ex.Range(1),Ex.Range(2),nPoints);
yy = zeros(1,nPoints);

S=1/2;
I=5/2;
SS1 = S*(S+1);
II1 = I*(I+1);

for M=-S:S-1
    % amplitude
    Amp1 = SS1-M*(M+1);
    %  third orded constants 
    Sm = SS1-M*(M-1); Sp = SS1-M*(M+1); Smp = Sm*Sp;
    Smm1 = Sm+2*M; Sm1 = Sm-2*M; Sm2 = Sm1-2*M; 
    Sp1 = SS1-2*(M+1); Sp2 = SS1-2*(M+1);
    Sm01 = Sm*Sm1; Sm12 = Sm1*Sm2; Smm10 = Smm1*Sm;
    Sp01 = Sp*Sp1; Sp12 = Sp1*Sp2; 
    Smp1 = Sm1 * Sp1;
    
    for m = -I:I
        for m1 = -I:I
            mm = [m;m1];
            Im = II1-mm.*(mm-1); Ip = II1-mm.*(mm+1);
%             H = (freq - a*(m+m1)-A*(1+2*M))/g/13.9962;
            H = (freq - sum(mm(:,on_s).*a,1))./g/13.9962;
            % second order corrections
%             H = H - (2*Bpm*(8*((M+1)^3-M*M*M)+1-4*SS1)-...
%                 2*Cpm*(2*((M+1)^3-M*M*M)+1-2*SS1))/g/freq/13.9962;
            %  third order
%             H = H - A.*(Bpm*(Sp1*(2*M+3)^3-Sp*(2*M+1)^3-Sm1*(2*M+1)^3+Sm*(2*M-1)^3)+...
%                 Cpm*((M+2)*Sp12-(M+1)*Sp01-M*Sm01+(M-1)*Smm10))/g/13.996/freq/freq;
%             H = H - real(Bp.*Bp.*Cm)*((2*M+3)*(2*M+5)*Sp12-(2*M+1)*(2*M+3)*Sp01+...
%                 (2*M+1)*(2*M-1)*Sm01-(2*M-1)*(2*M-3)*Smm10-...
%                 2*(2*M+1)*(2*M+3)*Smp1+2*(2*M-1)*(2*M+1)*Smp...
%                 )/g/13.996/freq/freq;
            % hpf second order corrections
            H = H - sum(a.*a/2.*(I*(I+1)+mm(:,on_s)*M*(M+1)-mm(:,on_s).*mm(:,on_s)),1)/g/freq/13.9962;
            % hpf third order corrections
            H = H - sum(a.*a.*a.*(((-M+m)*Sp1-(-M+m-1)*Sp).*Im(:,on_s)+...
                 ((M-mm(:,on_s))*Sm1-(M-mm(:,on_s)-1)*Sm).*Ip(:,on_s)),1)/4/g/freq/freq/13.9962;
            yy = binning(yy,xx,H,on_s*Amp1);
        end
    end
end

yy = convspec(yy,xx(2)-xx(1),lw,safeget(Ex,'Harmonic',0));

disp(sprintf('Exec time: %5.1fs',cputime-tt))