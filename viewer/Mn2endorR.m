%Mn at  Q  band
%------------------------------------------------------
function [x, y] = Mn2endor(Sys, Exp, Opt)

Npoints = Exp.nPoints;
width   = Exp.Range(2)-Exp.Range(1); %RF im MHz
step = width/Npoints;
pQ = Sys.Q*2/9; %quadrupole splittings in MHz
B0 = Exp.Field*10;
gx = Sys.g(3);
gy = Sys.g(2);
gz = Sys.g(1);
ax = Sys.Az;
ay = Sys.Ay;
az = Sys.Ax;
f0 = Exp.mwFreq*1E3;
w=Sys.lw*28;
tresW = 3*w;
vi = -B0*0.0011;
wN = Sys.lwEndor;
tresWN = 3*wN;
nTheta = Opt.nKnots*2;
nPhi = Opt.nKnots2;

dF = [1:Npoints]*step;
ENDOR = zeros(1,Npoints);


for m = 1:nTheta
  cT = m/nTheta;
  sT = sqrt(1-cT*cT);
  for q = 1:nPhi
    cP = cos(q*pi/2/nPhi);
    sP = sin(q*pi/2/nPhi);  
    
    a = ax*sT*sT*cP*cP + ay*sT*sT*sP*sP + az*cT*cT; %for all manganeses
    g = sqrt(gx*gx*sT*sT*cP*cP + gy*gy*sT*sT*sP*sP + gz*gz*cT*cT); 
   
    for k1 = 1:6
        m(1) = k1 - 3.5;
       for k2 = 1:6
          m(2) = k2 - 3.5;
              sec = sum(a.*a.*(35/4 - m.*m))/f0;
              f = B0*1.3996246*g + m*a'+sec;
              
              if abs(f-f0)<tresW
                weightE = exp(-(f-f0)^2/(w*w));
                
                for N = 1:2
                  Qu(N) = 3*pQ(N)*(2*m(N)+1)*(3*cT*cT-1)/2;
                  Qd(N) = 3*pQ(N)*(-2*m(N)+1)*(3*cT*cT-1)/2;
                  secpU(N) = -a(N)*abs(a(N))*(2*m(N)+1+1)/(f0*4);
                  secpD(N) = -a(N)*abs(a(N))*(2*m(N)-1+1)/(f0*4);
                  secmU(N) = -a(N)*abs(a(N))*(2*m(N)+1-1)/(f0*4);
                  secmD(N) = -a(N)*abs(a(N))*(2*m(N)-1-1)/(f0*4);
                end
               
                for p = 1:Npoints
                  for N = 1:2
                   if (m(N) ~= -2.5) %nuclear spin goes down
                      %electron spin goes down
                      fdm = abs(-a(N)/2-vi+Qd(N))+secmD(N);
                   if abs(fdm-dF(p))<tresWN ENDOR(p) = ENDOR(p) + weightE*exp(-((fdm-dF(p))/wN)^2); end
                      %electron spin goes up
                      fdp = abs(a(N)/2-vi+Qd(N))+secpD(N);
                   if abs(fdp-dF(p))<tresWN ENDOR(p) = ENDOR(p) + weightE*exp(-((fdp-dF(p))/wN)^2); end
                   end 
                  
                  if (m(N) ~= 2.5) %nuclear spin goes up 
                      %electron spin goes down
                      fum = abs(-a(N)/2+vi+Qu(N))+secpU(N);
                   if abs(fum-dF(p))<tresWN ENDOR(p) = ENDOR(p) + weightE*exp(-((fum-dF(p))/wN)^2); end
                      %electron spin is up
                      fup = abs(a(N)/2+vi+Qu(N))+secmU(N);
                   if abs(fup-dF(p))<tresWN ENDOR(p) = ENDOR(p) + weightE*exp(-((fup-dF(p))/wN)^2); end
                  end 
                 end 
                  
          end     
        end
     end
  end
 end
end    
    
ME = max(ENDOR);

x = dF;
y = real(ENDOR/ME);

EPRfig = safeget(Opt,'EPRfig',0);
if EPRfig > 0
    Exp.nPoints = Exp.nPointsEPR;
    [x1,y1] = Mn2EPR_R(Sys,Exp,Opt);
    figure(EPRfig);
    clf
    plot(x1/10,y1, 'r'); hold on;
    xlabel('Magnetic Field, mT');
    [expEPR.ax,expEPR.y]=asciiread(Opt.EPRexp,'\t',1);
    MS = max(expEPR.y);
    plot(expEPR.ax.x/10, expEPR.y/MS, 'b');
    axis tight
    for ii=1:length(Exp.Field)
        xpp = Exp.Field(ii);
        [xmin,xidx] = min(abs(expEPR.ax.x-xpp));
        plot(expEPR.ax.x(xidx),expEPR.y(xidx),'.')
    end
end