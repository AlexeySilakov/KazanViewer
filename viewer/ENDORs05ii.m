% ENDOR simulation with s=1/2 and 2 nuclei with
% colinear HPf and quarupole tenzors.
% ATTENTION ! contains mistake(s) 
% [xx,yy]=ENDORs05ii(Sy,Ex,Op)
function [xx, yy] = ENDORs05ii(Sy, Ex, Op)
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

fmin = planck*freq/bmagn/max(gg)*1E9-2*II*max(max(abs(AA)))*S/28-40;
fmax = planck*freq/bmagn/min(gg)*1E9+2*II*max(max(abs(AA)))*S/28+40;
fxx = linspace(fmin,fmax,1000);
fyy = zeros(1,1000);

for k1 = -II:II
    for k2 = -II:II
        m = [k1;k2];
        sec = sum(a.*a.*(35/4 - m.*m*on_s))/freq;
        
%       EPR  
        BB = (freq - sum(m(:,on_s).*a,1)-sec)./g/13.996246;
        fyy = binning(fyy,fxx,BB,on_s);

%       selection rules  
%         dB = abs(BB-B0);
%         idx = dB<(ExWidth*3/28);
%         Amp = exp(-(dB([1,1],idx)/ExWidth/28).^2);

        f = B0*13.996246*g + sum(m(:,on_s).*a,1)+sec;
        dF = abs(f-freq);
        idx = dF<(ExWidth*4);
        Amp = exp(-(dF([1,1],idx)/ExWidth).^2);

        mm = m(:,on_s);
        Qu = (2*mm([1,2],idx)+1).*q([1,2],idx);
        Qd = (-2*mm([1,2],idx)+1).*q([1,2],idx);
        secpU = -a([1,2],idx).*abs(a([1,2],idx)).*(2*mm([1,2],idx)+1+1)/(freq*4);
        secpD = -a([1,2],idx).*abs(a([1,2],idx)).*(2*mm([1,2],idx)-1+1)/(freq*4);
        secmU = -a([1,2],idx).*abs(a([1,2],idx)).*(2*mm([1,2],idx)+1-1)/(freq*4);
        secmD = -a([1,2],idx).*abs(a([1,2],idx)).*(2*mm([1,2],idx)-1-1)/(freq*4);

%       nuclear spin goes down
        midx = m~=-2.5;
%       electron spin goes down
        fdm = abs(-a(midx,idx)/2-vi+Qd(midx,:))+secmD(midx,:);
        yy = binning(yy,xx,fdm,Amp(midx,:));
%       electron spin goes up
        fdp = abs(a(midx,idx)/2-vi+Qd(midx,:))+secpD(midx,:);
        yy = binning(yy,xx,fdp,Amp(midx,:));
%       nuclear spin goes up 
        midx = m~=2.5;
%       electron spin goes down
        fum = abs(-a(midx,idx)/2+vi+Qu(midx,:))+secpU(midx,:);
        yy = binning(yy,xx,fum,Amp(midx,:));
%       electron spin is up
        fup = abs(a(midx,idx)/2+vi+Qu(midx,:))+secmU(midx,:);
        yy = binning(yy,xx,fup,Amp(midx,:));
    end
end

yy = convspec(yy,xx(2)-xx(1),lwENDOR,safeget(Ex,'Harmonic',0));

EPRfig = safeget(Op,'EPRfig',0);
if EPRfig > 0
    figure(EPRfig); 
    if Op.SimIdx == 1
        clf; hold on
        Harmonic = safeget(Ex,'HarmonicEPR',0);
        if isfield(Op, 'EPRexp')
%             [expEPR.ax,expEPR.y]=asciiread(Op.EPRexp,'\t',1);
            [expEPR.ax,expEPR.y]=brukerread(Op.EPRexp);
            plot(expEPR.ax.x/10, renorm(real(expEPR.y),1,0), 'b');
        end
        fyy = convspec(fyy,fxx(2)-fxx(1),lw,safeget(Ex,'HarmonicEPR',0));
        plot(fxx,renorm(fyy,1,0),'r'); 
        xlabel('Magnetic Field [mT]');
    end
    plot(B0*[1,1],[0,1],':')
    axis tight
end

disp(sprintf('Exec time: %5.1fs',cputime-tt))