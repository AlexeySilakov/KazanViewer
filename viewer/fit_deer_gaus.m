function y  = fit_deer_gaus(x, r0, dr, const, amp, varargin)
% function to fit deer time traces using gaussian distribution;
global kernel r lastx first_time;

if nargin>5
    fff = varargin{1};
else
    fff=0;
end
if isempty(first_time)
    fff = 1;
end
if ~fff
    if length(x)~=length(lastx)
        fff = 1;
    else
        res = find(x-lastx);
        fff = ~isempty(res);
    end
end
if fff
    first_time = 0;
    load('Pake_arrays.mat'); % Pake(r, t), t ns, r nm
    xmax = max(x); %
    xmin = min(x);
    [~, idx2] = min( abs(t-xmax)); idx2 = min(length(t), idx2+1);
    [~, idx1] = min( abs(t-xmin)); idx1 = max(1, idx1-1);
    kernel=Pake(:, idx1:idx2)-ones(size(Pake(:, idx1:idx2)));
    kernel = interp1(t(idx1:idx2), kernel.',  x);
    kernel = kernel.';
    lastx = x;
end

distr = exp(-(r-r0).^2/dr^2/4);
distr = 0.01*distr/sum(distr);
y = amp*exp(distr.'*kernel)+const;
y = y(:);
