% function [out] = renorm(in)
% function [out] = renorm(in, Amplitude)
% function [out] = renorm(in, Amplitude, 0) - preserve 0

function [out] = renorm(in, varargin)

if nargin>1,  ampl = varargin{1};
else, ampl = 1; end

maxx = max(in(:));
minn = min(in(:));

if nargin>2 & varargin{2}==0,  
    out = ampl * in./(max(abs([maxx,minn])));
else, 
    out = ampl * (in - minn)./(maxx-minn);
end

return