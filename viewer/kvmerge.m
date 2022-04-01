% function y = kvmerge(x,y1, x,y2, xs)
% merging of trace 1 + trace 2 atarting from xs
% data of second trace are extrapolated 

function [x,y] = kvmerge(x1,y1, x2,y2, xs)

xmin = min(x1);
xmax = max([max(x1),max(x2)]);
step = (max(x1)-xmin)/(length(x1)-1);
x = [xmin:step:xmax]';

idx1 = (x1 >= xmin) & (x1 <= xs);
idx2 = (x2 > xs) & (x2 <= xmax);

idxx1 = (x >= xmin) & (x <= xs);
idxx2 = (x > xs) & (x <= xmax);

y(idxx1,1) = y1(idx1);
y(idxx2,1) = spline(x2(idx2),y2(idx2),x(idxx2));

return

