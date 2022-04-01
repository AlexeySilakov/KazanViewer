function [newx, newy] = kv_spaceman(x, y, npoints)
% creates equally distributed X-axes and interpolate data correspondingly
% [newx, newy] = kv_spaceman(x, y, npoints)
% x,y - the original data.
% npoints - number of points for new axes.
[minx, min_ind] = min(x);
[maxx, max_ind] = max(x);
% correction alsi 17.05.06:
    [x, tx_ind] = sort(x);
    y = y(tx_ind);
        % check if x axes is uphill or downhill. if it is going down,
        % then reverse x and y. (problems with 'interp1q')
        %     if max_ind<min_ind, x = x(end:-1:1); y = y(end:-1:1); end
% correction.
newx = linspace(minx, maxx, npoints).';
% kk  =find((x(1:end-1) - x(2:end))~=0); % find dx=0 and KILL second points
% kk = [kk; length(y)]; % including last point
% txx = x(kk);
% tyy = y(kk);
% newy =  interp1q (txx, tyy, newx);
newy =  interp1q (x, y, newx);
% if there is downhill case, revese everything back
        %if max_ind<min_ind, newx = newx(end:-1:1); newy = newy(end:-1:1); end