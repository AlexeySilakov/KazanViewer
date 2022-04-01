% function varargout=kv_asciiwrite(filename, src, varargin)

function varargout=kv_asciiwrite(filename, src, varargin)

opt = struct('is2d',false);
if nargin >= 2,
    if ~mod(nargin-2,2)
        for kk=1:2:nargin-3
            opt=setfield(opt, lower(varargin{kk}), varargin{kk+1});    
        end
    else, error('Wrong amount of arguments')
    end
else, error('Usage: kv_asciiwrite(filename, src). Options: is2D.');
end

if opt.is2d, 
    p = real(src.y);
    dlmwrite(filename,p,' ')
else
    y  = real(src.y);
    y1 = imag(src.y);
    x  = src.ax.x;
    if sum(y1)
        p = [x,y,y1];
    else
        p = [x,y];
    end
    save(filename, 'p', '-ascii')    
end
