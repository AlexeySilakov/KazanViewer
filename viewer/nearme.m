function varargout = nearme(X, xi)
% nearme(X, xi)
% k = nearme(X, xi)
% [k, v] = nearme(X, xi)
% 
%   return value (v) and index (k) 
%   of the nearest to the xi number of the array X 

[c, k] = min(abs(X - ones(size(X))*xi));
v = X(k);
switch nargout 
case 0
    varargout{1} = k;
case 1
    varargout{1} = k;
case 2
    varargout{1} = k;
    varargout{2} = v; 
end