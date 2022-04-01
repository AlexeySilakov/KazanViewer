function [x, y] = Mn2EPR_R(Sys, Exp, Opt)


Npoints = Exp.nPoints;
width = 3000;
step = width/Npoints;
a1x = Sys.Ax(1);
a1y = Sys.Ay(1);
a1z = Sys.Az(1);
a2x = Sys.Ax(2);
a2y = Sys.Ay(2);
a2z = Sys.Az(2);
Bc = Exp.FieldC*10; %in Gauss
gx = Sys.g(1);
gy = Sys.g(2);
gz = Sys.g(3);
f0 = Exp.mwFreq*1e3;
w=Sys.lw*10;
tresW = 3*w;
nTheta = Opt.nKnots*2;
nPhi = Opt.nKnots2;

for p = 1:Npoints 
    field(p) = Bc + (p - Npoints/2)*step;
    sp(p) = 0;
end;

for m = 1:nTheta
  cT = (m-0.5)/nTheta;
  sT = sqrt(1-cT*cT);
  for q = 1:nPhi
    cP = cos(q*pi/2/nPhi);
    sP = sin(q*pi/2/nPhi);
    a1 = a1x*sT*sT*cP*cP + a1y*sT*sT*sP*sP + a1z*cT*cT;
    a2 = a2x*sT*sT*cP*cP + a2y*sT*sT*sP*sP + a2z*cT*cT;    
    g = sqrt(gx*gx*sT*sT*cP*cP + gy*gy*sT*sT*sP*sP + gz*gz*cT*cT);
    for k1 = 1:6
        for k2 = 1:6
              sec = a1*a1*(35/4-(k1-3.5)^2)+a2*a2*(35/4-(k2-3.5)^2);
              B0 = (f0 - ((k1-3.5)*a1 + (k2-3.5)*a2 )-sec/f0)/(g*1.3996246);
              lowP = max(1, round((B0 - tresW - Bc)/step + Npoints/2));
              highP = min(Npoints, round((B0 + tresW - Bc)/step + Npoints/2));
                for p = lowP:highP                       
                  if abs(B0-field(p))<tresW sp(p) = sp(p) + exp(-(B0-field(p))*(B0-field(p))/(w*w)); end
                end 
            
        end
    end
  end
end    

x = field;
ME = max(sp);
y = sp/ME;


