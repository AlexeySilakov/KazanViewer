% [xx,yy]=pepman(Sys,Exp,Op)
% Exact evaluation of Hamiltonian 
% without any interpolation using EasySpin

function [xx,yy]=pepman(Sy,Ex,Op)
tt=cputime;

nTheta = safeget(Ex,'nTheta',41);
nPhi   = safeget(Ex,'nPhi',11);

th  = linspace(0,pi/2,nTheta);
phi = linspace(0,pi,nPhi);

t=reshape(th(ones(1,nPhi),:),1,nTheta*nPhi);
p=reshape(phi(ones(1,nTheta),:).',1,nTheta*nPhi);

xx = linspace(Ex.Range(1),Ex.Range(2),Ex.nPoints);
yy = zeros(1,Ex.nPoints);

applyI = 0;
if isfield(Sy,'I'), Sy=rmfield(Sy,'I'); applyI=1; end

for i=1:length(p)
    [Pos,Amp] = eigfields(Sy,Ex,[p(i);t(i)]);
    yy = binning(yy,xx,Pos,Amp*t(i));
end

% if applyI
%     yy1 = zeros(1,Ex.nPoints);
%     a=safeget(Sy,'A',0*[1 1 1]);
%     a=sum(a)/length(a); 
%     hpfl=[1,2,3,4,5,6,5,4,3,2,1];
%     for ii=1:size(hpfl)
%         yy1 = yy1+spline(xx,yy,xx+(i-5)*a/28); 
%     end
%     yy=yy1;
% end

yy = convspec(yy,xx(2)-xx(1),Sy.lw,safeget(Ex,'Harmonic',0));
disp(sprintf('Exec time: %5.1fs',cputime-tt))
