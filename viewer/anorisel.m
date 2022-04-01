function [phi, theta] = anorisel(gtensor, g, nKnots, Ori)
% [phi, theta] = anorisel(gtensor, g, nKnots, Ori)
% analytical solution for the orientation selection problem in the first
% order. Only the g-tensor is used.

%tg = sort(gtensor);
tg = gtensor;
a = tg(1); b = tg(2); c = tg(3);
r = g;

if r<min(tg)
  %  disp('g-value is out of range of the given g-tensor. Closest principal g-value has been taken'); 
    r = min(tg);

end
if  r>max(tg),   
 %   disp('g-value is out of range of the given g-tensor. Closest principal g-value has been taken'); 
    r = max(tg);
end   
if r == min(tg)
    theta = pi/2*ones(nKnots*2, 1);
    phi = zeros(nKnots*2, 1);
    return;
end
if r == max(tg)
    theta = zeros(nKnots*2, 1);
    phi = zeros(nKnots*2, 1);
    return;
end

% z bonds:
%-z^2*b^2+b^2*c^2+c^2*z^2-c^2*rr >0
% if r < b, zb1 == 0
% tz = z/c ->cos(theta);
zb1 = (r^2-b^2)/(c^2-b^2);
if zb1<=0, zb1 = 0; end
zb1 = sqrt(zb1);
%-z^2*a^2+c^2*z^2-c^2*rr+a^2*c^2<0
zb2 = (r^2-a^2)/(c^2-a^2);
if zb2<=0, zb2 = 0; end
zb2 = sqrt(zb2);
% zb1 < z < zb2
th = linspace(acos(zb1), acos(zb2), floor(nKnots));
for ii = 1:nKnots
    z = c*cos(th(ii));
    x = 1/(-b^2+a^2)*(-(-b^2+a^2)*(-z^2*b^2+b^2*c^2+c^2*z^2-c^2*r^2))^(1/2)*a/c;
    ph(ii) = acos(x/sqrt(1-cos(th(ii))^2)/a);
end
theta = [th, th(end:-1:1)]'; %, -th, -th(end:-1:1)];
phi = [ph, -ph(end:-1:1)]'; %, ph, -ph(end:-1:1)];
    % [X, Y, Z] = sphere(50);
    % surf(X*r, Y*r, Z*r); hold on;
    % [X, Y, Z] = ellipsoid(0,0,0,a,b,c,50);
    % surf(X, Y, Z);
    % axis equal;
    % hold on;
    % plot3(a*cos(ph).*sin(th), b*sin(ph).*sin(th), c*cos(th), 'g');
    %y = 1/(-b^2+a^2)*((-b^2+a^2)*(-z^2*a^2+c^2*z^2-c^2*rr+a^2*c^2))^(1/2)*b/c;